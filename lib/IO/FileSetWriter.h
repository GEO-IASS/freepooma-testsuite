// -*- C++ -*-
//
// Copyright (C) 1998, 1999, 2000, 2002  Los Alamos National Laboratory,
// Copyright (C) 1998, 1999, 2000, 2002  CodeSourcery, LLC
//
// This file is part of FreePOOMA.
//
// FreePOOMA is free software; you can redistribute it and/or modify it
// under the terms of the Expat license.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Expat
// license for more details.
//
// You should have received a copy of the Expat license along with
// FreePOOMA; see the file LICENSE.
//

//-----------------------------------------------------------------------------
// This code is based on the DiscField implementation by Bill Humphrey.
//-----------------------------------------------------------------------------

/** @file
 * @ingroup IO
 * @brief
 * FileSetWriter<Dim> manages the writing of Arrays and Fields from
 * "DiscField" format files.
 */

#ifndef POOMA_IO_FILESETWRITER_H
#define POOMA_IO_FILESETWRITER_H

// Include files

#include "Array/Array.h"
#include "Engine/CompressibleBrick.h"
#include "Engine/RemoteEngine.h"
#include "Domain/Interval.h"
#include "Pooma/Pooma.h"
#include "Utilities/PAssert.h"

#include <fstream>
#include <string.h>  // memset
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////
//
// A class for writing DiskField fileSets. 
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
class FileSetWriter {
public:

  //---------------------------------------------------------------------------
  // Constructors.
  
  // This constructor is used when creating a .data and a .offset file
  // for writing. 
  
  FileSetWriter(const char *base, int fieldsPerRecord);

  //---------------------------------------------------------------------------
  // Destructor.
  
  ~FileSetWriter();

  //---------------------------------------------------------------------------
  // Modifiers.
  
  // Write the .data and .offset parts of the fileset.
  // Handles both Fields and Arrays.

  template<class Subject>
  void write(Subject &subject);

  //---------------------------------------------------------------------------
  // Typedefs

  // Public since we want to share it with DFOffsetData, and template
  // friends may not be supported everywhere. (This will go away when
  // we marshal the structure in a standard way.)

  typedef POOMA_INT64 Offset_t;
  
private:
    
  template<class T>
  struct DFOffsetData 
  {
    int        vnodedata_m[6 * Dim];
    union {
      bool       isCompressed_m;
      char       pad[8];
    } u;
    Offset_t   offset_m;
    T          compressedVal_m;
  };

  template<class Layout>
  void initializeRecord(const Layout &layout);
  
  template<class T>
  void writePatch(const Array<Dim, T, Remote<CompressibleBrick> > &a);

  template<class T>
  void writeValue(std::ofstream &fs, const T &val);

  template<class T>
  void writeValues(std::ofstream &fs, const T *vals, int n);

  void writeValue(std::ofstream &fs, const Interval<Dim> &val);
  
  void writeMetaFile();

  std::string baseFileName_m;
  std::vector<Interval<Dim> > domains_m;
  std::vector<int> numPatches_m;
  Interval<Dim> overallDomain_m;
  Offset_t currentOffset_m;
  int ioContext_m;
  int fieldsPerRecord_m, currentRecord_m, currentField_m;
  std::ofstream data_m, offset_m, layout_m;
};

////////////////////////////////////////////////////////////////////////////
//
// Sets up writing to this fileSet by intializing variables and opening
// files. 
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
FileSetWriter<Dim>::FileSetWriter(const char *base, int fieldsPerRecord)
: baseFileName_m(base), 
  currentOffset_m(0),
  fieldsPerRecord_m(fieldsPerRecord),
  currentRecord_m(0),
  currentField_m(0)
{
  // For now, we are always writing on context 0.
  
  ioContext_m = 0;

  // Only open the files if on the context doing I/O.
  
  if (Pooma::context() == ioContext_m)
    {
      std::string dataFileName = baseFileName_m + ".data";
      std::string offsetFileName = baseFileName_m + ".offset";
      std::string layoutFileName = baseFileName_m + ".layout";

      data_m.open(dataFileName.c_str(), std::ios::binary | std::ios::trunc);
      PInsist(data_m.is_open(), 
        "FileSetWriter::FileSetWriter - Couldn't create data file.");

      offset_m.open(offsetFileName.c_str(), std::ios::binary | std::ios::trunc);
      PInsist(offset_m.is_open(), 
        "FileSetWriter::FileSetWriter - Couldn't create offset file.");

      layout_m.open(layoutFileName.c_str(), std::ios::binary | std::ios::trunc);
      PInsist(layout_m.is_open(), 
        "FileSetWriter::FileSetWriter - Couldn't create layout file.");
    }
}

////////////////////////////////////////////////////////////////////////////
//
// The destructor closes our files if they are open.
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
FileSetWriter<Dim>::~FileSetWriter()
{
}

////////////////////////////////////////////////////////////////////////////
//
// User-callable write function: handles both Fields and Arrays.
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
template<class Subject>
void FileSetWriter<Dim>::write(Subject &subject)
{    
  // Make sure we have the correct dimension.
  
  PInsist(Dim == Subject::dimensions,
    "FileSetWriter::write - dimensions doesn't match DiskField.");

  // If this is the first record, we need to note the overal physical
  // domain.
  
  if (currentRecord_m == 0)
    overallDomain_m = subject.layout().innerDomain();
    
  // If this is the first field in the record, we need initialize some data. 

  if (currentField_m == 0)  
    initializeRecord(subject.layout());

  // Write the subject to the .data and .offset files.
  
  // Find out the number of materials and centering points.

  int nMaterials = numMaterials(subject);
  int nCentering = centeringSize(subject);
    
  // Make sure that this write will fit entirely inside a record.
  
  PInsist(currentField_m + nMaterials * nCentering <= fieldsPerRecord_m,
    "FileSetWriter::write - Too many fields in the record.");

  for (int m = 0; m < nMaterials; m++)
    {
      for (int c = 0; c < nCentering; c++)
        {
          // We need to write the current field# to the .offset file.

          writeValue(offset_m, currentField_m); 

          // Get a subfield view of the subject.
                   
          Subject s(subField(subject, m, c));
          
          // Now, we need to loop over all of the domains in the layout.
          
          for (int i = 0; i < domains_m.size(); i++)
            {
              // We have a slightly sticky problem in that the layout
              // holds vertex domains, which might not be corect for
              // some sub-fields (with, for example, cell centering).
              
              // We get around this by intersecting the domain from
              // the layout with the physical domain from the sub-field.
              
              // This is a no-op for arrays.
              
              Interval<Dim> d = intersect(domains_m[i], s.domain());
              
              // Create an array to receive the data. Make this
              // remote (owned by the I/O context) and compressible.
              
              Array<Dim, typename Subject::Element_t, Remote<CompressibleBrick> > a;
              a.engine() = 
                Engine<Dim, typename Subject::Element_t, Remote<CompressibleBrick> >
                  (ioContext_m, d);
                
              // Assign to the array to get the data. One might be able
              // to not do this if the subject already had no guard layers
              // and was on context 0 already. 
              
              a = s(d);
              
              // We need to make sure this assignment is done before
              // proceeding.
              
              Pooma::blockAndEvaluate();
              
              // Write this "vnode" to the .data and .offset files.
              
              // NOTE: this means that the "redundant" data in the
              // .offset and .layout files is not identical for non-vertex
              // centered fields.
              
              writePatch(a);
            }
                   
          currentField_m++;
        }
    }

  // Move on to the next record, if necessary. Write the .meta file here.
      
  if (currentField_m == fieldsPerRecord_m)
    {
      currentField_m = 0;
      currentRecord_m++;

      writeMetaFile();    
    }
}  

////////////////////////////////////////////////////////////////////////////
//
// Initializes the data for a record, writing out part of the .layout
// file in the process.
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
template<class Layout>
void FileSetWriter<Dim>::initializeRecord(const Layout &layout)
{
  // First, clear out the list of domains.
  
  domains_m.clear();
  
  // Now, go through the layout, getting out the nodes and intersecting their
  // owned domains, which include global guards that we don't want to write,
  // with the overall physical domain.
  
  typename Layout::const_iterator i = layout.beginGlobal();
  while (i != layout.endGlobal())
    {
      domains_m.push_back(intersect(i->domain(), overallDomain_m));
      ++i;
    }
   
  // The rest of this stuff is only done in the conect doing the I/O.
   
  if (Pooma::context() == ioContext_m)
    {
      // Store the number of patches (aka "vnodes") in this record.
      // We will need this for the .meta file.
       
      numPatches_m.push_back(domains_m.size());
       
      // Write the number of patches to the .layout file.
       
      writeValue<int>(layout_m, domains_m.size());
       
      // Write the domains.
       
      for (int i = 0; i < domains_m.size(); i++)
        writeValue(layout_m, domains_m[i]);
    }
} 
       
////////////////////////////////////////////////////////////////////////////
//
// Low level routine for writing patch data using the infamous DFOffsetData
// struct.
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
template<class T>
void FileSetWriter<Dim>::writePatch(const Array<Dim, T, Remote<CompressibleBrick> > &a)
{  
  // We only do this stuff if our context is doing I/O.
  
  if (Pooma::context() != ioContext_m)
    return;
    
  // We need to fill in this nutty structure.
    
  DFOffsetData<T> odata;

  memset(&odata, 0, sizeof(DFOffsetData<T>));  
  for (int i = 0; i < Dim; i++)
    {
      // NOTE: this ordering disagrees with the printed DiskField
      // documentation, but it DOES match the code. Bad Bill! Bad! :-)

      // Also, the offset is in terms of ELEMENTS, not bytes.

      odata.vnodedata_m[i * 6 + 1] = a.domain()[i].first();
      odata.vnodedata_m[i * 6 + 2] = 1;
      odata.vnodedata_m[i * 6 + 3] = a.domain()[i].size();
    }

  odata.u.isCompressed_m = compressed(a);
  
  if (odata.u.isCompressed_m)
    {
      odata.offset_m = 0;
      odata.compressedVal_m = a.engine().localEngine().compressedRead();
    }
  else
    {
      odata.offset_m = currentOffset_m;
      odata.compressedVal_m = 0;

      // Might as well just write the data here.
            
      writeValues(data_m, a.engine().localEngine().dataBlock().beginPointer(),
        a.domain().size());
      currentOffset_m += a.domain().size();
    }
    
  // Write to the .offset file.
  
  writeValue(offset_m, odata);
}

////////////////////////////////////////////////////////////////////////////
//
// Low level routines for writing a value or values to a file. Only
// executed on context doing I/O. 
//
// Also, performs error checking for write problems. 
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
template<class T>
void FileSetWriter<Dim>::
writeValue(std::ofstream &fs, const T &val)
{
  // We only want to write if we are on the I/O context.
  
  if (Pooma::context() == ioContext_m)
    {
      fs.write((char *) &val, sizeof(T));
    }
}

template<int Dim>
template<class T>
void FileSetWriter<Dim>::
writeValues(std::ofstream &fs, const T *vals, int n)
{
  // We only want to write if we are on the I/O context.
  
  if (Pooma::context() == ioContext_m)
    {
      fs.write((char *) vals, n * sizeof(T));
    }
}

template<int Dim>
void FileSetWriter<Dim>::
writeValue(std::ofstream &fs, const Interval<Dim> &val)
{
  // We only want to write if we are on the I/O context.
  
  if (Pooma::context() == ioContext_m)
    {
      for (int i = 0; i < Dim; i++)
        {
          // NOTE: this ordering disagrees with the printed DiskField
          // documentation, but it DOES match the code. Bad Bill! Bad! :-)
          
          writeValue<int>(fs, 0);
          writeValue<int>(fs, val[i].first());
          writeValue<int>(fs, 1);
          writeValue<int>(fs, val[i].size());
          writeValue<int>(fs, 0);
          writeValue<int>(fs, 0);
        }
    }
}

////////////////////////////////////////////////////////////////////////////
//
// Writes the .meta file. 
//
////////////////////////////////////////////////////////////////////////////

template<int Dim>
void FileSetWriter<Dim>::writeMetaFile()
{
  // We only want to write if we are on the I/O context.
  
  if (Pooma::context() != ioContext_m)
    return;

  std::string metaFileName(baseFileName_m);
  metaFileName += ".meta";
  
  std::ofstream fout(metaFileName.c_str(), std::ios::out | std::ios::trunc);
  PInsist(fout.is_open(), 
    "FileSetWriter - couldn't create .meta file.");
    
  fout << "Type           = unknown\n";
  fout << "Dim            = " << Dim << "\n";
  for (int i = 0; i < Dim; i++)
    {
      fout << "Domain         = "
        << overallDomain_m[i].first() << " "
        << overallDomain_m[i].last() << " "
        << "1\n";
    }
  fout << "Fields         = " << fieldsPerRecord_m << "\n";
  fout << "Records        = " << currentRecord_m << "\n";
  fout << "SMPs           = " << 1 << "\n";
  fout << "VnodesInRecord =";
  for (int j = 0; j < numPatches_m.size(); j++)
    {
      fout << " ";
      fout << numPatches_m[j];
    }
  fout << "\n";
  fout << "VnodeTally     =";
  int tally = 0;
  for (int k = 0; k < numPatches_m.size(); k++)
    {
      fout << " ";
      fout << tally;
      tally += numPatches_m[k];
    }
  fout << "\n";
}
    
#endif // POOMA_IO_FILESETWRITER_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: FileSetWriter.h,v $   $Author: richard $
// $Revision: 1.14 $   $Date: 2004/11/01 18:16:52 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
