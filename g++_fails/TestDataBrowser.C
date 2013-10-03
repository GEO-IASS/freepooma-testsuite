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

// ----------------------------------------------------------------------------
// TestDataBrowser.cpp , Tim Williams 10/20/1999
// 
// Tests some of the dbprint() functions, for Fields, Arrays, and
// DynamicArrays (DynamicArrays not yet put in --TJW
// 10/22/1999). Self-checking via comparison with hardcoded output checked by
// hand to be correct. Also includes some example function prototypes for
// calling print functions interactively from the debugger; must run under
// debugger and reset Inform object to one that outputs to the screen to test
// these.
// ----------------------------------------------------------------------------

// Include files:
#include "Pooma/Fields.h"
#include "Pooma/Particles.h"
#include "Utilities/Tester.h"
#include "DataBrowser/DataBrowser.h"

#include <fstream>

// Function prototypes for functions defined below main():
void hardCodedOutput(char* filename);
bool thediff(char* filename1, char* filename2);

// Particles class for testing DataBrowser printing functions on attributes:
const int PDim = 2;
class Specks : public Particles<MPDynamicUniform>
{
public:
  typedef Particles<MPDynamicUniform>   Base_t;
  typedef Base_t::AttributeEngineTag_t  AttributeEngineTag_t;
  typedef Base_t::ParticleLayout_t      ParticleLayout_t;
  typedef double                        AxisType_t;
  typedef Vector<PDim,AxisType_t>       PointType_t;
  // Constructor: set up layouts, register attributes
  Specks(const ParticleLayout_t &pl)
  : Particles<MPDynamicUniform>(pl)
  {
    addAttribute(pos);
    addAttribute(vel);
  }
  // Position and velocity attributes (as public members)
  DynamicArray< PointType_t, AttributeEngineTag_t >  pos;
  DynamicArray< PointType_t, AttributeEngineTag_t >  vel;
};


// Global typedefs; useful in making user-defined functions below:
// 1D
typedef UniformRectilinearMesh<1> Mesh1_t;
typedef Field<Mesh1_t, double> ScalarField1_t;
typedef Field<Mesh1_t, Vector<1> > VectorField1_t;
typedef Array<1, double, CompressibleBrick> ScalarArray1_t;
typedef Array<1, Vector<1>, CompressibleBrick> VectorArray1_t;
// 2D
typedef Array<2, double, CompressibleBrick> ScalarArray2_t;
typedef Array<2, Vector<2>, CompressibleBrick> VectorArray2_t;
typedef Array<2, Tensor<2,double,Antisymmetric>, Brick> TensorArray2_t;
// 3D
typedef Array<3, double, CompressibleBrick> ScalarArray3_t;
// 4D
typedef Array<4, double, CompressibleBrick> ScalarArray4_t;
// Particle attributes:
typedef DynamicArray<Specks::PointType_t, Specks::AttributeEngineTag_t> 
Attribute_t;

// User-defined nontemplate dbprint()-type functions:
void sfdbprint(const ScalarField1_t &f) { dbprint(f); }
void vfdbprint(const VectorField1_t &f) { dbprint(f); }
void sadbprint(const ScalarArray1_t &a) { dbprint(a); }
void vadbprint(const VectorArray1_t &a) { dbprint(a); }
void sa2dbprint(const ScalarArray2_t &a) { dbprint(a); }
void va2dbprint(const VectorArray2_t &a) { dbprint(a); }
void ta2dbprint(const TensorArray2_t &a) { dbprint(a); }
void pdbprint(const Attribute_t &pa) { dbprint(pa); }
// Subsetting functions:
// N.B.: these have to have separate names; some debuggers aren't smart enough
// to understand multiple prototypes of function with same name.
void esfdbprint(const ScalarField1_t &f, int i) { dbprint(f,i); }
void rsfdbprint(const ScalarField1_t &f, int ibase, int ibound,
                int istride) { dbprint(f,ibase,ibound,istride); }
void epdbprint(const Attribute_t &pa, int i) { dbprint(pa, i); }
void rpdbprint(const Attribute_t &pa, int base, int bound, int stride)
{ dbprint(pa, base,bound,stride); }

int main(int argc, char* argv[])
{
  // Initialize POOMA and Tester class.
  Pooma::initialize(argc, argv);
  Pooma::Tester tester(argc, argv);

  // Create an Inform object and attach to data browser
  Inform* foutp = 
    new Inform(NULL,"text.test.TestDataBrowser",Inform::out,0);
  Inform &fout = *foutp;
  dbSetInform(fout);
  dbSetCarReturn(3);
  dbSetDataPrecision(6);
  //  dbSetScientific(true);
  dbSetDataWidth(15);

  // 1D Vertex and cell domains:
  int nVerts = 9;
  int nCells = nVerts - 1;
  Interval<1> vertDomain(nVerts), cellDomain(nCells);

  // Create the 1D mesh; default origin and spacings:
  Mesh1_t mesh(vertDomain);
  
  // Create the 1D geometry:
  Centering<1> cell = canonicalCentering<1>(CellType, Continuous);
  DomainLayout<1> layout(vertDomain);
  
  fout << std::endl << "=========== 1D ============" << std::endl;

  // Make some 1D fields:
  ScalarField1_t s1(cell, layout, mesh);
  VectorField1_t v1(cell, layout, mesh);

  // Assign to spatially-varying values:
  s1.all() = positions(s1).comp(0);
  v1.all() = positions(v1);

  // Create some 1D Arrays:
  ScalarArray1_t sa1(cellDomain);
  VectorArray1_t va1(cellDomain);

  // Assign to spatially-varying values:
  sa1 = s1;
  va1 = v1;

  // Make sure assignments have completed before printing values
  Pooma::blockAndEvaluate();

  // Output the 1D Fields:
  dbprint(s1);
  dbprint(v1);

  // Output one of the 1D Arrays:
  Range<1> ss;
  for (int d = 0; d < 1; d++) {
    ss[d] = Range<1>(1, nCells - 2, 2);
  }
  dbprint(sa1(ss));
  dbprint(sa1, ss);
  dbprint(sa1, 0,3,2);
  dbprint(sa1,1);

  // 2D:
  fout << std::endl << "=========== 2D ============" << std::endl;

  ScalarArray2_t sa2(cellDomain, cellDomain);
  VectorArray2_t va2(cellDomain, cellDomain);
  TensorArray2_t ta2(cellDomain, cellDomain);
  for (int ix = 0; ix < nCells; ix++) {
    sa2(ix, cellDomain) = sa1;
    va2(ix, cellDomain).comp(0) = sa1;
    va2(ix, cellDomain).comp(1) = sa1;
    for (int i = 1; i < 2; i++) {
      for (int j = 0; j < i; j++) {
        ta2(ix, cellDomain).comp(i,j) = sa1;
      }
    }
  }

  // Make sure assignments have completed before printing values
  Pooma::blockAndEvaluate();

  dbprint(sa2);
  dbSetCarReturn(2);
  dbprint(va2);
  dbprint(sa2,1,1);
  dbprint(va2,1,1);
  dbSetDataWidth(5);
  dbprint(ta2);

  // 3D:
  fout << std::endl << "=========== 3D ============" << std::endl;

  ScalarArray3_t sa3(cellDomain, cellDomain, cellDomain);
  for (int i = 0; i < nCells; i++) {
    for (int j = 0; j < nCells; j++) {
      sa3(i, j, cellDomain) = sa1;
    }
  }

  // Make sure assignments have completed before printing values
  Pooma::blockAndEvaluate();

  dbSetDataWidth(10);
  dbprint(sa3,1,1,1);
//   dbprint(sa3);

  // 4D:
  fout << std::endl << "=========== 4D ============" << std::endl;

  ScalarArray4_t sa4(cellDomain, cellDomain, cellDomain, cellDomain);
  for (int k = 0; k < nCells; k++) {
    for (int j = 0; j < nCells; j++) {
      for (int i = 0; i < nCells; i++) {
        sa4(i, j, k, cellDomain) = sa1;
      }
    }
  }

  // Make sure assignments have completed before printing values
  Pooma::blockAndEvaluate();

  dbprint(sa4,1,1,1,1);
  Interval<4> ijkl;
  ijkl[0] = cellDomain; 
  ijkl[1] = cellDomain;
  ijkl[2] = Interval<1>(1,2); 
  ijkl[3] = Interval<1>(1,2);
  dbprint(sa4, ijkl);
  dbprint(sa4, 0, nCells-1, 1, 0, nCells-1, 1, 3, 4, 1, 1, 2, 1);
//   dbprint(sa3);

  // Particles (2D):
  fout << std::endl << "=========== 2D Particles============" << std::endl;

  Specks::ParticleLayout_t particleLayout(17); // 17 patches
  Specks specks(particleLayout);
  const int np = 50;
  specks.globalCreate(np);
  for (int p = 0; p < np; p++) {
    specks.pos(p) = Specks::PointType_t(1.0*p);
    specks.vel(p) = Specks::PointType_t(2.0*p);
  }

  dbprint(specks.pos);
  dbprint(specks.vel);
  dbprint(specks.vel, 0, 23, 3);

  // Write out "by hand" into another file what the previous field-printing
  // functions should have produced; this will be compared with what they
  // actually did produce:
  hardCodedOutput("text.correct.TestDataBrowser");

  // Compare the two files by mocking up the Unix "diff" command:
  delete foutp;
  tester.check(thediff("text.test.TestDataBrowser", 
                       "text.correct.TestDataBrowser"));

  int retval = tester.results("TestDataBrowser");
  Pooma::finalize();
  return retval;
}

//-----------------------------------------------------------------------------
// Mock up the Unix "diff" utility to compare two files:
//-----------------------------------------------------------------------------
bool thediff(char* filename1, char* filename2)
{
  bool same = true;
  char ch1, ch2;
  std::ifstream file1(filename1);
  std::ifstream file2(filename2);
  while (file1.get(ch1)) {          // Read file 1 char-by-char until eof
    if (file2.get(ch2)) {           // Read equivalent char from file 2
      if (ch1 != ch2) same = false; // If they're different,files are different
    } else {
      same = false;                 // If file 2 ends before file 1, different
    }
  }
  return(same);
}

//-----------------------------------------------------------------------------
void hardCodedOutput(char* filename)
{
  std::ofstream of(filename);
  of << "" << std::endl;
  of << "=========== 1D ============" << std::endl;
  of << "( -2:009:001) =            -1.5            -0.5             0.5" << std::endl;
  of << "                            1.5             2.5             3.5" << std::endl;
  of << "                            4.5             5.5             6.5" << std::endl;
  of << "                            7.5             8.5             9.5" << std::endl;
  of << "( -2:009:001) = (           -1.5) (           -0.5) (            0.5)" << std::endl;
  of << "                (            1.5) (            2.5) (            3.5)" << std::endl;
  of << "                (            4.5) (            5.5) (            6.5)" << std::endl;
  of << "                (            7.5) (            8.5) (            9.5)" << std::endl;
  of << "(000:002:001) =             1.5             3.5             5.5" << std::endl;
  of << "(001:005:002) =             1.5             3.5             5.5" << std::endl;
  of << "(000:002:002) =             0.5             2.5" << std::endl;
  of << "(001) =             1.5" << std::endl;
  of << "" << std::endl;
  of << "=========== 2D ============" << std::endl;
  of << "(000:007:001,000) =             0.5             0.5             0.5" << std::endl;
  of << "                                0.5             0.5             0.5" << std::endl;
  of << "                                0.5             0.5" << std::endl;
  of << "(000:007:001,001) =             1.5             1.5             1.5" << std::endl;
  of << "                                1.5             1.5             1.5" << std::endl;
  of << "                                1.5             1.5" << std::endl;
  of << "(000:007:001,002) =             2.5             2.5             2.5" << std::endl;
  of << "                                2.5             2.5             2.5" << std::endl;
  of << "                                2.5             2.5" << std::endl;
  of << "(000:007:001,003) =             3.5             3.5             3.5" << std::endl;
  of << "                                3.5             3.5             3.5" << std::endl;
  of << "                                3.5             3.5" << std::endl;
  of << "(000:007:001,004) =             4.5             4.5             4.5" << std::endl;
  of << "                                4.5             4.5             4.5" << std::endl;
  of << "                                4.5             4.5" << std::endl;
  of << "(000:007:001,005) =             5.5             5.5             5.5" << std::endl;
  of << "                                5.5             5.5             5.5" << std::endl;
  of << "                                5.5             5.5" << std::endl;
  of << "(000:007:001,006) =             6.5             6.5             6.5" << std::endl;
  of << "                                6.5             6.5             6.5" << std::endl;
  of << "                                6.5             6.5" << std::endl;
  of << "(000:007:001,007) =             7.5             7.5             7.5" << std::endl;
  of << "                                7.5             7.5             7.5" << std::endl;
  of << "                                7.5             7.5" << std::endl;
  of << "(000:007:001,000) = (            0.5,            0.5) (            0.5,            0.5)" << std::endl;
  of << "                    (            0.5,            0.5) (            0.5,            0.5)" << std::endl;
  of << "                    (            0.5,            0.5) (            0.5,            0.5)" << std::endl;
  of << "                    (            0.5,            0.5) (            0.5,            0.5)" << std::endl;
  of << "(000:007:001,001) = (            1.5,            1.5) (            1.5,            1.5)" << std::endl;
  of << "                    (            1.5,            1.5) (            1.5,            1.5)" << std::endl;
  of << "                    (            1.5,            1.5) (            1.5,            1.5)" << std::endl;
  of << "                    (            1.5,            1.5) (            1.5,            1.5)" << std::endl;
  of << "(000:007:001,002) = (            2.5,            2.5) (            2.5,            2.5)" << std::endl;
  of << "                    (            2.5,            2.5) (            2.5,            2.5)" << std::endl;
  of << "                    (            2.5,            2.5) (            2.5,            2.5)" << std::endl;
  of << "                    (            2.5,            2.5) (            2.5,            2.5)" << std::endl;
  of << "(000:007:001,003) = (            3.5,            3.5) (            3.5,            3.5)" << std::endl;
  of << "                    (            3.5,            3.5) (            3.5,            3.5)" << std::endl;
  of << "                    (            3.5,            3.5) (            3.5,            3.5)" << std::endl;
  of << "                    (            3.5,            3.5) (            3.5,            3.5)" << std::endl;
  of << "(000:007:001,004) = (            4.5,            4.5) (            4.5,            4.5)" << std::endl;
  of << "                    (            4.5,            4.5) (            4.5,            4.5)" << std::endl;
  of << "                    (            4.5,            4.5) (            4.5,            4.5)" << std::endl;
  of << "                    (            4.5,            4.5) (            4.5,            4.5)" << std::endl;
  of << "(000:007:001,005) = (            5.5,            5.5) (            5.5,            5.5)" << std::endl;
  of << "                    (            5.5,            5.5) (            5.5,            5.5)" << std::endl;
  of << "                    (            5.5,            5.5) (            5.5,            5.5)" << std::endl;
  of << "                    (            5.5,            5.5) (            5.5,            5.5)" << std::endl;
  of << "(000:007:001,006) = (            6.5,            6.5) (            6.5,            6.5)" << std::endl;
  of << "                    (            6.5,            6.5) (            6.5,            6.5)" << std::endl;
  of << "                    (            6.5,            6.5) (            6.5,            6.5)" << std::endl;
  of << "                    (            6.5,            6.5) (            6.5,            6.5)" << std::endl;
  of << "(000:007:001,007) = (            7.5,            7.5) (            7.5,            7.5)" << std::endl;
  of << "                    (            7.5,            7.5) (            7.5,            7.5)" << std::endl;
  of << "                    (            7.5,            7.5) (            7.5,            7.5)" << std::endl;
  of << "                    (            7.5,            7.5) (            7.5,            7.5)" << std::endl;
  of << "(001,001) =             1.5" << std::endl;
  of << "(001,001) = (            1.5,            1.5)" << std::endl;
  of << "(000:007:001,000) = ((    0  -0.5)(  0.5     0)) ((    0  -0.5)(  0.5     0))" << std::endl;
  of << "                    ((    0  -0.5)(  0.5     0)) ((    0  -0.5)(  0.5     0))" << std::endl;
  of << "                    ((    0  -0.5)(  0.5     0)) ((    0  -0.5)(  0.5     0))" << std::endl;
  of << "                    ((    0  -0.5)(  0.5     0)) ((    0  -0.5)(  0.5     0))" << std::endl;
  of << "(000:007:001,001) = ((    0  -1.5)(  1.5     0)) ((    0  -1.5)(  1.5     0))" << std::endl;
  of << "                    ((    0  -1.5)(  1.5     0)) ((    0  -1.5)(  1.5     0))" << std::endl;
  of << "                    ((    0  -1.5)(  1.5     0)) ((    0  -1.5)(  1.5     0))" << std::endl;
  of << "                    ((    0  -1.5)(  1.5     0)) ((    0  -1.5)(  1.5     0))" << std::endl;
  of << "(000:007:001,002) = ((    0  -2.5)(  2.5     0)) ((    0  -2.5)(  2.5     0))" << std::endl;
  of << "                    ((    0  -2.5)(  2.5     0)) ((    0  -2.5)(  2.5     0))" << std::endl;
  of << "                    ((    0  -2.5)(  2.5     0)) ((    0  -2.5)(  2.5     0))" << std::endl;
  of << "                    ((    0  -2.5)(  2.5     0)) ((    0  -2.5)(  2.5     0))" << std::endl;
  of << "(000:007:001,003) = ((    0  -3.5)(  3.5     0)) ((    0  -3.5)(  3.5     0))" << std::endl;
  of << "                    ((    0  -3.5)(  3.5     0)) ((    0  -3.5)(  3.5     0))" << std::endl;
  of << "                    ((    0  -3.5)(  3.5     0)) ((    0  -3.5)(  3.5     0))" << std::endl;
  of << "                    ((    0  -3.5)(  3.5     0)) ((    0  -3.5)(  3.5     0))" << std::endl;
  of << "(000:007:001,004) = ((    0  -4.5)(  4.5     0)) ((    0  -4.5)(  4.5     0))" << std::endl;
  of << "                    ((    0  -4.5)(  4.5     0)) ((    0  -4.5)(  4.5     0))" << std::endl;
  of << "                    ((    0  -4.5)(  4.5     0)) ((    0  -4.5)(  4.5     0))" << std::endl;
  of << "                    ((    0  -4.5)(  4.5     0)) ((    0  -4.5)(  4.5     0))" << std::endl;
  of << "(000:007:001,005) = ((    0  -5.5)(  5.5     0)) ((    0  -5.5)(  5.5     0))" << std::endl;
  of << "                    ((    0  -5.5)(  5.5     0)) ((    0  -5.5)(  5.5     0))" << std::endl;
  of << "                    ((    0  -5.5)(  5.5     0)) ((    0  -5.5)(  5.5     0))" << std::endl;
  of << "                    ((    0  -5.5)(  5.5     0)) ((    0  -5.5)(  5.5     0))" << std::endl;
  of << "(000:007:001,006) = ((    0  -6.5)(  6.5     0)) ((    0  -6.5)(  6.5     0))" << std::endl;
  of << "                    ((    0  -6.5)(  6.5     0)) ((    0  -6.5)(  6.5     0))" << std::endl;
  of << "                    ((    0  -6.5)(  6.5     0)) ((    0  -6.5)(  6.5     0))" << std::endl;
  of << "                    ((    0  -6.5)(  6.5     0)) ((    0  -6.5)(  6.5     0))" << std::endl;
  of << "(000:007:001,007) = ((    0  -7.5)(  7.5     0)) ((    0  -7.5)(  7.5     0))" << std::endl;
  of << "                    ((    0  -7.5)(  7.5     0)) ((    0  -7.5)(  7.5     0))" << std::endl;
  of << "                    ((    0  -7.5)(  7.5     0)) ((    0  -7.5)(  7.5     0))" << std::endl;
  of << "                    ((    0  -7.5)(  7.5     0)) ((    0  -7.5)(  7.5     0))" << std::endl;
  of << "" << std::endl;
  of << "=========== 3D ============" << std::endl;
  of << "(001,001,001) =        1.5" << std::endl;
  of << "" << std::endl;
  of << "=========== 4D ============" << std::endl;
  of << "(001,001,001,001) =        1.5" << std::endl;
  of << "" << std::endl;
  of << "~~~~~~~~~~~~~~ (0:7:1,0:7:1,1:2:1,1:2:1) ~~~~~~~~~~~~~~" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,1,1):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,001,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,002,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,003,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,004,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,005,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,006,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,007,001,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,2,1):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,001,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,002,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,003,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,004,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,005,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,006,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,007,002,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,1,2):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,001,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,002,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,003,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,004,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,005,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,006,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,007,001,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,2,2):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,001,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,002,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,003,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,004,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,005,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,006,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,007,002,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "" << std::endl;
  of << "~~~~~~~~~~~~~~ (0:7:1,0:7:1,3:4:1,1:2:1) ~~~~~~~~~~~~~~" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,3,1):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,001,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,002,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,003,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,004,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,005,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,006,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,007,003,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,4,1):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,001,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,002,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,003,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,004,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,005,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,006,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "(000:007:001,007,004,001) =        1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "                                   1.5        1.5" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,3,2):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,001,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,002,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,003,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,004,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,005,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,006,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,007,003,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "" << std::endl;
  of << "(0:7:1,0:7:1,4,2):" << std::endl;
  of << "----------------------------------------------------" << std::endl;
  of << "(000:007:001,000,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,001,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,002,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,003,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,004,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,005,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,006,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "(000:007:001,007,004,002) =        2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "                                   2.5        2.5" << std::endl;
  of << "" << std::endl;
  of << "=========== 2D Particles============" << std::endl;
  of << "(000:049:001) = (         0,         0) (         1,         1)" << std::endl;
  of << "                (         2,         2) (         3,         3)" << std::endl;
  of << "                (         4,         4) (         5,         5)" << std::endl;
  of << "                (         6,         6) (         7,         7)" << std::endl;
  of << "                (         8,         8) (         9,         9)" << std::endl;
  of << "                (        10,        10) (        11,        11)" << std::endl;
  of << "                (        12,        12) (        13,        13)" << std::endl;
  of << "                (        14,        14) (        15,        15)" << std::endl;
  of << "                (        16,        16) (        17,        17)" << std::endl;
  of << "                (        18,        18) (        19,        19)" << std::endl;
  of << "                (        20,        20) (        21,        21)" << std::endl;
  of << "                (        22,        22) (        23,        23)" << std::endl;
  of << "                (        24,        24) (        25,        25)" << std::endl;
  of << "                (        26,        26) (        27,        27)" << std::endl;
  of << "                (        28,        28) (        29,        29)" << std::endl;
  of << "                (        30,        30) (        31,        31)" << std::endl;
  of << "                (        32,        32) (        33,        33)" << std::endl;
  of << "                (        34,        34) (        35,        35)" << std::endl;
  of << "                (        36,        36) (        37,        37)" << std::endl;
  of << "                (        38,        38) (        39,        39)" << std::endl;
  of << "                (        40,        40) (        41,        41)" << std::endl;
  of << "                (        42,        42) (        43,        43)" << std::endl;
  of << "                (        44,        44) (        45,        45)" << std::endl;
  of << "                (        46,        46) (        47,        47)" << std::endl;
  of << "                (        48,        48) (        49,        49)" << std::endl;
  of << "(000:049:001) = (         0,         0) (         2,         2)" << std::endl;
  of << "                (         4,         4) (         6,         6)" << std::endl;
  of << "                (         8,         8) (        10,        10)" << std::endl;
  of << "                (        12,        12) (        14,        14)" << std::endl;
  of << "                (        16,        16) (        18,        18)" << std::endl;
  of << "                (        20,        20) (        22,        22)" << std::endl;
  of << "                (        24,        24) (        26,        26)" << std::endl;
  of << "                (        28,        28) (        30,        30)" << std::endl;
  of << "                (        32,        32) (        34,        34)" << std::endl;
  of << "                (        36,        36) (        38,        38)" << std::endl;
  of << "                (        40,        40) (        42,        42)" << std::endl;
  of << "                (        44,        44) (        46,        46)" << std::endl;
  of << "                (        48,        48) (        50,        50)" << std::endl;
  of << "                (        52,        52) (        54,        54)" << std::endl;
  of << "                (        56,        56) (        58,        58)" << std::endl;
  of << "                (        60,        60) (        62,        62)" << std::endl;
  of << "                (        64,        64) (        66,        66)" << std::endl;
  of << "                (        68,        68) (        70,        70)" << std::endl;
  of << "                (        72,        72) (        74,        74)" << std::endl;
  of << "                (        76,        76) (        78,        78)" << std::endl;
  of << "                (        80,        80) (        82,        82)" << std::endl;
  of << "                (        84,        84) (        86,        86)" << std::endl;
  of << "                (        88,        88) (        90,        90)" << std::endl;
  of << "                (        92,        92) (        94,        94)" << std::endl;
  of << "                (        96,        96) (        98,        98)" << std::endl;
  of << "(000:021:003) = (         0,         0) (         6,         6)" << std::endl;
  of << "                (        12,        12) (        18,        18)" << std::endl;
  of << "                (        24,        24) (        30,        30)" << std::endl;
  of << "                (        36,        36) (        42,        42)" << std::endl;
  of.close();
  return;
}
// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: TestDataBrowser.cpp,v $   $Author: richard $
// $Revision: 1.9 $   $Date: 2004/11/01 18:16:30 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
