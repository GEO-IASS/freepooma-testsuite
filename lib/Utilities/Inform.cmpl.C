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

//-----------------------------------------------------------------------------
// Inform non-template definitions.
//-----------------------------------------------------------------------------

// POOMA include files
#include "Utilities/Inform.h"
#include "Utilities/PAssert.h"

// string-based ostream subclass include files
#include <iomanip>
#include <fstream>

#if POOMA_NO_STRINGSTREAM
#include <strstream>
#else
#include <sstream>
#endif

// the size of the Inform internal character buffer,
// for those implementations which do not have ostringstream

const unsigned int Inform::bufSize = 32000;

// A mutex used to protect printing to the output stream, since that can
// be shared among many Inform's

Pooma::Mutex_t Inform::outputMutex_s;

// the current context ID and number of contexts for all Inform
// objects.  By default it looks like we're running on context 0
// of 1 total contexts.  But this can be changed by the underlying
// run-time system.

Inform::Context_t Inform::context_s = 0;
Inform::Context_t Inform::nContexts_s = 1;


//-----------------------------------------------------------------------------
// The InformStream class, defined here because it is only used internally
// by Inform.  InformStream stores information about a single stream
// connection.  The information is: ostream pointer, whether it needs to
// be closed by us, the destination context, and output level threshhold.
//-----------------------------------------------------------------------------

class InformStream
{
public:
  // Constructor which takes the ostream to use and the destination context
  InformStream(std::ostream *s, Inform::Context_t oc)
    : stream_m(s), close_m(false), outputContext_m(oc), level_m(0)
  {
    PAssert(s != 0);
  }

  // Constructor which takes a filename, and opens a file.
  InformStream(const char *fname, int mode, Inform::Context_t oc)
    : stream_m(0), close_m(true), outputContext_m(oc), level_m(0)
  {
    PAssert(fname != 0);
    PAssert(mode == Inform::out || mode == Inform::app);

    if (oc < 0 || oc == Inform::context()) {
      if (mode == Inform::out)
	stream_m = new std::ofstream(fname, std::ios::out);
      else
	stream_m = new std::ofstream(fname, std::ios::app);
    }
  }

  // Destructor: close the stream, if necessary.
  ~InformStream()
  {
    if (close_m && stream_m != 0)
      delete stream_m;
  }

  // get or set the output context
  Inform::Context_t outputContext() const
  {
    return outputContext_m;
  }

  void setOutputContext(Inform::Context_t val)
  { 
    outputContext_m = val;
  }

  // get or set the output threshold
  Inform::Level_t outputLevel() const
  {
    return level_m;
  }

  void setOutputLevel(Inform::Level_t val)
  {
    level_m = val;
  }

  // print out the given message to the output stream
  void print(Inform::Level_t l, const std::string &prefix, const char *msg)
  {
    if (shouldPrint(l))
      {
	// Print the prefix, if necessary.
        if (prefix.length() > 0)
          {
            // First print the base part of the prefix
            *stream_m << prefix;

            // If there are more than 1 contexts, insert the context ID
            // into the prefix:
            if (Inform::numContexts() > 1)
              {
                if ((outputContext_m == Inform::allContexts) ||
                    (outputContext_m == Inform::context()))
                  {
                    *stream_m << "{" << Inform::context() << "}";
                  }
              }

            // Add on the "> ":
            *stream_m << "> ";
          }

        // Then print the rest of the message and flush the stream
        *stream_m << msg << "\n";
        stream_m->flush();
      }
  }

private:
  // The stream to manage.
  std::ostream *stream_m;

  // Do we need to close (delete) this stream when done?
  bool close_m;

  // Which context should we write to.
  Inform::Context_t outputContext_m;

  // The output message threshold level
  Inform::Level_t level_m;

  // A boolean function that determines if we should print out the
  // current message based on:
  //   1. The output level
  //   2. The current context settings
  //   3. Do we have somewhere to print to
  bool shouldPrint(Inform::Level_t level)
  {
    PAssert(level >= 0);

    // Definitely do not print if there is no stream
    if (stream_m == 0)
      return false;

    // Do not print if the specified message level does not match
    // our output level requirements, or we're on the wrong context.
    return (level <= level_m &&
	    (outputContext_m == Inform::context() ||
	     outputContext_m == Inform::allContexts));
  }
};


//-----------------------------------------------------------------------------
// Inform endl manipulator.  When 'endl' is inserted into an Inform stream,
// it causes the current message to be printed, along with a newline.
//
// Inform flush manipulator.  Just calls 'flush' on the Inform.
//
// Inform lock and unlock manipulator.  Just calls the corresponding
// routines to lock access to the stream.
//-----------------------------------------------------------------------------

namespace std {

  Inform &endl(Inform &inf)
  {
    inf.flush();
    return inf;
  }

  Inform &flush(Inform &inf)
  {
    inf.flush();
    return inf;
  }

  Inform &lock(Inform &inf)
  {
    inf.lock();
    return inf;
  }

  Inform &unlock(Inform &inf)
  {
    inf.unlock();
    return inf;
  }

} // namespace std

//-----------------------------------------------------------------------------
// Create an Inform object which will print to just standard out with
// the given prefix and destination context (initially these are
// defaulted to 'no prefix' and 'just print on context 0').
// The initial output stream has an ID value of '0'.
//-----------------------------------------------------------------------------

Inform::Inform(const char *prefix, Context_t outputContext)
  : prefix_m(""), outputContext_m(outputContext), level_m(0),
    message_m(0), buffer_m(0), nextID_m(0)
{
  // Create a connection to cout
  open(outputContext);

  // Initialize our internal formatting buffer
  setup(prefix);
}


//-----------------------------------------------------------------------------
// Create an Inform object which will print to a file with the given
// name, opened either for overwrite or append operations.  The first
// and last arguments give the prefix and destination context as usual.
// If the destination context is 'allContexts', then a file will be
// created by all contexts.  If the destination context is just a single
// context, then only that context will have a file opened.
// The destination context for a file cannot be changed once it is set
// in the constructor or 'open' call.
// The initial output stream has an ID value of '0'.
//-----------------------------------------------------------------------------

Inform::Inform(const char *prefix, const char *fname, int writemode,
	       Context_t outputContext)
  : prefix_m(""), outputContext_m(outputContext), level_m(0),
    message_m(0), buffer_m(0), nextID_m(0)
{
  // Create a connection to the given file
  open(fname, writemode, outputContext);

  // Initialize our internal formatting buffer
  setup(prefix);
}


//-----------------------------------------------------------------------------
// Create an Inform object which will print to the given ostream,
// with a prefix and destination context.  The destination context
// for this case CAN be changed later.
// The initial output stream has an ID value of '0'.
//-----------------------------------------------------------------------------

Inform::Inform(const char *prefix, std::ostream &outstream, 
               Context_t outputContext)
  : prefix_m(""), outputContext_m(outputContext), level_m(0),
    message_m(0), buffer_m(0), nextID_m(0)
{
  // Create a connection to the given file
  open(outstream, outputContext);

  // Initialize our internal formatting buffer
  setup(prefix);
}


//-----------------------------------------------------------------------------
// The Inform destructor will flush all existing buffers, and
// close all files that it opened.
//-----------------------------------------------------------------------------

Inform::~Inform()
{
  // delete all the existing connections
  close();

  // delete the existing output buffer stream
  delete message_m;

  // delete the storage buffer, if necessary
  if (buffer_m != 0)
    delete [] buffer_m;
}


//-----------------------------------------------------------------------------
// Open a connection to a new stream, and return the ID for that
// stream.  The ID should be used in other manipulator calls, such
// as 'outputLevel(ID_t)'.  There are three forms for open, which
// correspond to the three types of constructors without prefixes:
//   1. output context: open new standard-out connection.
//   2. filename + output mode + output context: open new file.
//   3. ostream + output context: open new ostream connection.
// Upon successful completion, open returns the ID for the new connection.
// If an error occurs, open returns (-1).
//-----------------------------------------------------------------------------

Inform::ID_t Inform::open(Context_t oc)
{
  streams_m.insert(Value_t(nextID_m, new InformStream(&std::cout, oc)));
  return nextID_m++;
}

Inform::ID_t Inform::open(const char *fname, int mode, Context_t oc)
{
  streams_m.insert(Value_t(nextID_m, new InformStream(fname, mode, oc)));
  return nextID_m++;
}

Inform::ID_t Inform::open(std::ostream &outstream, Context_t oc)
{
  streams_m.insert(Value_t(nextID_m, new InformStream(&outstream, oc)));
  return nextID_m++;
}


//-----------------------------------------------------------------------------
// Close the specifed connection.
//-----------------------------------------------------------------------------

void Inform::close(ID_t id)
{
  // find the proper connection to close
  iterator s = streams_m.find(id);
  PAssert(s != streams_m.end());

  // delete the connection, and remove it from the map
  delete ((*s).second);
  streams_m.erase(s);
}


//-----------------------------------------------------------------------------
// Close all connections.
//-----------------------------------------------------------------------------

void Inform::close()
{
  // delete all the existing connections
  for (iterator a = streams_m.begin(); a != streams_m.end(); ++a)
    delete ((*a).second);

  // remove all the elements
  streams_m.erase(streams_m.begin(), streams_m.end());
}


//-----------------------------------------------------------------------------
// Return the current value for the output threshhold level, which is
// the highest level of message that will be printed.  If this is < 0,
// no messages will be printed.  You must specify which stream to return
// the settings for, which has a default value of '0'.
//-----------------------------------------------------------------------------

Inform::Level_t Inform::outputLevel(ID_t id) const
{
  // Find the proper connection
  InformStream *s = findStream(id);
  PAssert(s != 0);

  // Return the output level
  return s->outputLevel();
}


//-----------------------------------------------------------------------------
// Change the output threshhold level for the output stream specified
// in the second argument.  If the first argument is < 0, then this
// effectively turns off that stream, since a message level cannot be
// < 0.  The 'off' enumeration can be used for the first argument to
// indicate this.  If no ID is given, change the setting for all.
//-----------------------------------------------------------------------------

void Inform::setOutputLevel(Level_t newval, ID_t id)
{
  // Find the proper connection
  InformStream *s = findStream(id);
  PAssert(s != 0);

  // Change the output level
  s->setOutputLevel(newval);
}

void Inform::setOutputLevel(Level_t newval)
{
  for (iterator a = streams_m.begin(); a != streams_m.end(); ++a)
    (*a).second->setOutputLevel(newval);
}


//-----------------------------------------------------------------------------
// Return the current value for the destination context, for the specified
// ostream connection.  The value is either >= 0, indicating a single
// destination context, or it is < 0, indicating 'all contexts'.
//-----------------------------------------------------------------------------

Inform::Context_t Inform::outputContext(ID_t id) const
{
  // Find the proper connection
  InformStream *s = findStream(id);
  PAssert(s != 0);

  // Return the output context
  return s->outputContext();
}


//-----------------------------------------------------------------------------
// Change the destination context for the specified ostream connection.
// Note that for some destinations, the context cannot be changed.
// If ID is not given, change the context for all.
//-----------------------------------------------------------------------------

void Inform::setOutputContext(Context_t outputContext, ID_t id)
{
  // Find the proper connection
  InformStream *s = findStream(id);
  PAssert(s != 0);

  // Change the output context
  s->setOutputContext(outputContext);
}

void Inform::setOutputContext(Context_t outputContext)
{
  for (iterator a = streams_m.begin(); a != streams_m.end(); ++a)
    (*a).second->setOutputContext(outputContext);
}


//-----------------------------------------------------------------------------
// Print out the current message to the active streams.
//-----------------------------------------------------------------------------

void Inform::flush()
{
  // terminate the existing string
  *message_m << std::ends;

  outputMutex_s.lock();

  // get a pointer to the formatted string buffer
#if POOMA_NO_STRINGSTREAM
  std::string formatstr = buffer_m;
#else
  std::string formatstr = message_m->str();
#endif

  // create a buffer that we can modify
  char *outputbuf = new char[formatstr.length() + 2];

  // scan through the lines of the buffer, as separated by newlines,
  // and print out each line
  char *endbuf = outputbuf;
  char *begbuf = outputbuf;
  const char *formatbuf = formatstr.c_str();

  do {
    // advance endbuf to next set of endlines
    while (*formatbuf != '\n' && *formatbuf != '\0')
      *endbuf++ = *formatbuf++;

    // make sure there is a string terminator
    *endbuf = '\0';

    // if there was a newline, skip it in the formatbuf
    if (*formatbuf == '\n')
      ++formatbuf;

    // go through each connection, and print the prefix and message.
    for (iterator a = streams_m.begin(); a != streams_m.end(); ++a)
      (*a).second->print(level_m, prefix_m, begbuf);

    // skip begbuf to the beginning of the next line
    begbuf = endbuf;
  } while (*formatbuf != '\0');

  // reset the string buffer
#if POOMA_NO_STRINGSTREAM
  message_m->seekp(0, std::ios::beg);
#else
  message_m->str( std::string() );
#endif

  // delete the buffer copy and we're done
  delete [] outputbuf;

  outputMutex_s.unlock();
}


//-----------------------------------------------------------------------------
// Return a pointer to the InformStream with the given ID.  If it is
// not found, 0 is returned.
//-----------------------------------------------------------------------------

InformStream *Inform::findStream(ID_t id) const
{
  // try to find the ID in the map
  const_iterator s = streams_m.find(id);

  // return either the pointer, or 0 if not found
  if (s != streams_m.end())
    return (*s).second;
  else
    return 0;
}


//-----------------------------------------------------------------------------
// Perform setup information needed by each constructor
//-----------------------------------------------------------------------------

void Inform::setup(const char *prefix)
{
  // Initialize the prefix string
  setPrefix(prefix);

#if POOMA_NO_STRINGSTREAM
  // create an ostrstream, using a fixed-size character buffer for storage
  buffer_m = new char[Inform::bufSize];
  *buffer_m = '\0';
  message_m = new std::ostrstream(buffer_m, Inform::bufSize, std::ios::out);
#else
  // create an ostringstream, which uses its own string for storage
  buffer_m = 0;
  message_m = new std::ostringstream(std::ios::out);
#endif
}


//-----------------------------------------------------------------------------
// Change the prefix string to the given value
//-----------------------------------------------------------------------------

void Inform::setPrefix(const char *prefix)
{
  // Initialize the prefix string.  If prefix is empty or null, do not use a
  // prefix at all.  Otherwise, make it of the form "prefix{context}> "

  if (prefix == 0 || *prefix == '\0')
    prefix_m = "";
  else
    prefix_m = prefix;
}


// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Inform.cmpl.cpp,v $   $Author: richard $
// $Revision: 1.21 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
