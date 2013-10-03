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

#ifndef POOMA_UTILITIES_INFORM_H
#define POOMA_UTILITIES_INFORM_H

//-----------------------------------------------------------------------------
// Classes:
// Inform
//-----------------------------------------------------------------------------

/** @file
 * @ingroup Utilities
 * @brief
 * A general utility class which looks very much like an ostream,
 * which will format output to include an optional prefix string, and will
 * print out results to multiple other ostreams.
 *
 * When used in a parallel
 * environment, Inform will also print out the context number as part of
 * the prefix.
 */

//-----------------------------------------------------------------------------
// Include Files
//-----------------------------------------------------------------------------

#include "Pooma/Configuration.h"
#include "Utilities/PAssert.h"
#include "Threads/PoomaMutex.h"

#include <iostream>
#include <string>
#include <map>

#if POOMA_NO_STRINGSTREAM
# include <strstream>
#else
# include <sstream>
#endif

//-----------------------------------------------------------------------------
// Forward References
//-----------------------------------------------------------------------------

class InformStream;


///////////////////////////////////////////////////////////////////////////////
// namespace POOMA {

/**
 * A message is sent to an Inform object by treating it as an ostream,
 * then ending the message by sending the 'inform' manipulator.  In
 * fact, Inform works much like an ostream, although it may actually
 * just use stdio for I/O.
 *
 * Each message is assigned the current 'level of interest'; the lower
 * the level, the more important it is.  Each Inform object is also
 * set for a current level; messages with a level <= the current level
 * are displayed.  Level values >= 0 should be used to print values;
 * setting the output threshhold level to be < 0 will turn off printing
 * of all messages.
 *
 * By default, a new Inform object will only print out the message on
 * context 0.  You may change the node on which this prints with the
 * 'printContext(int)' method; if the argument is 'allContexts',
 * the message will be printed on ALL nodes, not just one.
 * Or, 'printAllContexts()' can be used instead.  The final
 * argument to the constructor may also be set to the context to print on.
 */

class Inform
{
public:
  //============================================================
  // Typedefs and enumerations
  //============================================================

  // The type of data used to specify a stream in accessor functions.  This
  // for now is just an int, used to look up the stream in the list.
  typedef int ID_t;

  // A typedef for the data used to indicate level values, and a typedef
  // for context values.
  typedef int Level_t;
  typedef int Context_t;

  // Enumeration listing the ways in which a file may be opened for writing,
  // the code used to indicate 'all contexts', and the code used to indicate
  // 'messages off' when setting the output level.
  enum { out, app };
  enum { allContexts = (-1) };
  enum { off = (-1), on = 0 };


  //============================================================
  // Constructors
  //============================================================

  // Create an Inform object which will print to just standard out with
  // the given prefix and destination context (initially these are
  // defaulted to 'no prefix' and 'just print on context 0').
  // The initial output stream has an ID value of '0'.
  Inform(const char *prefix = 0, Context_t outputContext = 0);

  // Create an Inform object which will print to a file with the given
  // name, opened either for overwrite or append operations.  The first
  // and last arguments give the prefix and destination context as usual.
  // If the destination context is 'allContexts', then a file will be
  // created by all contexts.  If the destination context is just a single
  // context, then only that context will have a file opened.
  // The destination context for a file cannot be changed once it is set
  // in the constructor or 'open' call.
  // The initial output stream has an ID value of '0'.
  Inform(const char *prefix, const char *fname, int writemode,
	 Context_t outputContext = 0);

  // Create an Inform object which will print to the given ostream,
  // with a prefix and destination context.  The destination context
  // for this case CAN be changed later.
  // The initial output stream has an ID value of '0'.
  Inform(const char *prefix, std::ostream &outstream,
	 Context_t outputContext = 0);


  //============================================================
  // Destructor
  //============================================================

  // The Inform destructor will flush all existing buffers, and
  // close all files that it opened.
  ~Inform();


  //============================================================
  // Prefix string access/modify methods
  //============================================================

  // Return the current prefix string:
  const std::string &prefix() const { return prefix_m; }

  // Change the prefix string to the given value, or empty if the
  // argument is null.

  void setPrefix(const char *prefix = 0);


  //============================================================
  // Output stream open/close methods
  //============================================================

  // Open a connection to a new stream, and return the ID for that
  // stream.  The ID should be used in other manipulator calls, such
  // as 'outputLevel(ID_t)'.  There are three forms for open, which
  // correspond to the three types of constructors without prefixes:
  //   1. output context: open new standard-out connection.
  //   2. filename + output mode + output context: open new file.
  //   3. ostream + output context: open new ostream connection.
  // Upon successful completion, open returns the ID for the new connection.
  // If an error occurs, open returns (-1).

  ID_t open(Context_t context = 0);
  ID_t open(const char *fname, int writemode, Context_t context = 0);
  ID_t open(std::ostream &outstream, Context_t context = 0);

  // Close the specifed connection.
  void close(ID_t);

  // Close all connections.
  void close();


  //============================================================
  // Inform message level methods
  //============================================================

  // Return the current value for the message level, which is the level
  // of the message that is currently being created.  This level will
  // be compared with the output threshold level for each active stream,
  // and if it is <= the threshold, the message will be printed.
  Level_t messageLevel() const { return level_m; }

  // Change the current value for the message level.
  Inform &setMessageLevel(Level_t newval) { level_m = newval; return *this; }


  //============================================================
  // Inform output threshold level methods
  //============================================================

  // Return the current value for the output threshhold level, which is
  // the highest level of message that will be printed.  If this is < 0,
  // no messages will be printed.  You must specify which stream to return
  // the settings for, which has a default value of '0'.
  Level_t outputLevel(ID_t id = 0) const;

  // Change the output threshhold level for the output stream specified
  // in the second argument.  If the first argument is < 0, then this
  // effectively turns off that stream, since a message level cannot be
  // < 0.  The 'off' enumeration can be used for the first argument to
  // indicate this.  If no ID value is given, change it for all connections.
  void setOutputLevel(Level_t newval, ID_t id);
  void setOutputLevel(Level_t newval);


  //============================================================
  // Inform destination context methods
  //============================================================

  // Return the current value for the destination context, for the specified
  // ostream connection.  The value is either >= 0, indicating a single
  // destination context, or it is < 0, indicating 'all contexts'.
  Context_t outputContext(ID_t id = 0) const;

  // Change the destination context for the specified ostream connection.
  // Note that for some destinations, the context cannot be changed.
  // If no ID value is given, change it for all connections.
  void setOutputContext(Context_t outputContext, ID_t id);
  void setOutputContext(Context_t outputContext);


  //============================================================
  // Inform creator context methods
  //============================================================

  // Return the current value for the creator's context, and the
  // total number of contexts that Inform objects believe exist.
  // These methods are static since they work with static data.
  inline static Context_t context() { return context_s; }
  inline static Context_t numContexts() { return nContexts_s; }

  // Set the current and total number of contexts for all Inform
  // objects.  These methods are static since they change static
  // data.  This is generally only done once, during initialization
  // of the runtime system of whatever program is using Inform objects.
  static inline void setContext(Context_t c) { context_s = c; }
  static inline void setNumContexts(Context_t n) { nContexts_s = n; }


  //============================================================
  // Inform basic operations
  //============================================================

  // Print out the current message to the active streams.
  void flush();
  void print() { flush(); }
  void output() { flush(); }
  
#if !POOMA_NO_IOSBASE_FMTFLAGS
  typedef std::ios_base::fmtflags FmtFlags_t;
#else
  typedef long FmtFlags_t;
#endif

  // return a reference to the internal ostream used to print messages
  std::ostream& stream() { return *message_m; }

  // functions used to change format state; used just as for iostreams
  FmtFlags_t 
    setf(FmtFlags_t setbits,FmtFlags_t field) 
    { return message_m->setf(setbits,field);}
  FmtFlags_t 
    setf(FmtFlags_t f) { return message_m->setf(f); }
  void /*long*/ unsetf(FmtFlags_t f) { message_m->unsetf(f); }
  long flags() const { return message_m->flags(); }
  long flags(FmtFlags_t f) { return message_m->flags(f); }
  int width() const { return message_m->width(); }
  int width(int w) { return message_m->width(w); }
  char fill() const { return message_m->fill(); }
  char fill(char c) { return message_m->fill(c); }
  int precision() const { return message_m->precision(); }
  int precision(int p) { return message_m->precision(p); }

  //============================================================
  // Mutex functions
  //============================================================

  void lock()   const { mutex_m.lock(); }
  void unlock() const { mutex_m.unlock(); }

private:
  //============================================================
  // Private typedefs and enumerations
  //============================================================

  // The type of storage used to hold the list of InformStream objects,
  // with a Size_t type and iterator types.
  typedef std::map<ID_t, InformStream *> StreamList_t;
  typedef StreamList_t::size_type        Size_t;
  typedef StreamList_t::value_type       Value_t;
  typedef StreamList_t::iterator         iterator;
  typedef StreamList_t::const_iterator   const_iterator;


  //============================================================
  // Private data members
  //============================================================

  // The name of this object; put at the start of each message.
  std::string prefix_m;

  // The context of the creator of this Inform object, set on construction:
  Context_t outputContext_m;

  // The current message level
  Level_t level_m;

  // The list of output destinations
  StreamList_t streams_m;

  // The ostringstream which sends the text to be printed out to a string.
  // We use an ostrstream if ostringstream is not available.

  //--------------------------------------------------------------------------
  // Gross hack alert!
  // Because of incompatibilities in various string stream classes,
  // we store a pointer to a type that depends on the platform.
  //--------------------------------------------------------------------------

#if POOMA_NO_STRINGSTREAM
  std::ostrstream *message_m;
#else
  std::ostringstream *message_m;
#endif

  // A character buffer used as storage by the ostrstream, if needed.
  char *buffer_m;

  // The fixed size of the character buffer.
  static const unsigned int bufSize;

  // The next ID value to use
  ID_t nextID_m;

  // A mutex for use in printing to this stream from multiple threads
  mutable Pooma::Mutex_t mutex_m;

  // A mutex used to protect printing to just the output streams.  This
  // is used in the "flush" method.
  static Pooma::Mutex_t outputMutex_s;

  // The local context number.  This is a static value that defaults to
  // zero; if you are using Inform in a parallel environment, then after
  // initializing the parallel env somebody should call "setContext(int)"
  // to tell this Inform (and all others) what this context number is.
  // We cannot determine it ourselves since it is runtime-system-specific.
  // The same is true for nContexts_s, the total number of contexts.
  static Context_t context_s;
  static Context_t nContexts_s;

  //============================================================
  // Private methods
  //============================================================

  // Return a pointer to the InformStream with the given ID.  If it is
  // not found, 0 is returned.
  InformStream *findStream(ID_t) const;

  // Perform setup information needed by each constructor
  void setup(const char *prefix);
};


//-----------------------------------------------------------------------------
// Inform manipulators
//-----------------------------------------------------------------------------

namespace std {

  // manipulator for signaling we want to send the message.
  extern Inform &endl(Inform &);
  extern Inform &flush(Inform &);

  // manipulator for signaling to lock/unlock Inform streams
  extern Inform &lock(Inform &);
  extern Inform &unlock(Inform &);

} // namespace std

// specialized version of operator<< to handle Inform-specific manipulators
inline Inform &operator<<(Inform &o, Inform &(*d)(Inform &))
{
  return d(o);
}

#if POOMA_NO_STD_IOSBASE
// specialized version of operator<< to handle ios manipulators
inline Inform &operator<<(Inform &o, ios &(*d)(ios &))
{
  d(o.stream());
  return o;
}
#else // !POOMA_NO_STD_IOSBASE
// specialized version of operator<< to handle ios_base manipulators
inline Inform &operator<<(Inform &o, std::ios_base &(*d)(std::ios_base &))
{
  d(o.stream());
  return o;
}
#endif // POOMA_NO_STD_IOSBASE


//-----------------------------------------------------------------------------
// Templated version of operator<< for Inform objects.  If you try to
// print an object to an Inform instance, and no operator<< is defined for
// that object that takes an Inform as the first argument, this version of
// operator<< will be used to just print the object to the internal ostream
// of the Inform.  Thus, for any class, you only need to define
//    ostream &operator<<(ostream &o, const classname &instance)
// which will then work with both Inform and general ostreams.
//
// If you are having problems with an ambiguity between another class's
// version of operator<< and this version, you can do the following:
//   1. Add a templated print(StreamType &) method to your class.
//   2. include <iosfwd> in your header
//   3. define operator<< taking an ostream &, and then call the print
//      method in the original class
//-----------------------------------------------------------------------------

template<class T>
inline Inform &operator<<(Inform &o, const T &val)
{
  o.stream() << val;
  return o;
}


//-----------------------------------------------------------------------------
// specialized version of operator<< to handle void * arguments
//-----------------------------------------------------------------------------

inline Inform &operator<<(Inform &o, const void *val)
{
  Inform::FmtFlags_t oldformat = 
    o.setf(std::ios::hex, std::ios::basefield);
  o.stream() << "0x" << reinterpret_cast<long>(val);
  o.setf(oldformat, std::ios::basefield);
  return o;
}


//-----------------------------------------------------------------------------
// specialized version of operator<< to handle long long type
//-----------------------------------------------------------------------------

#if defined(_LONGLONG)
inline Inform &operator<<(Inform &o, long long val)
{
  // cast to long double before sending to ostream
  o.stream() << static_cast<long double>(val);
  return o;
}
#endif // _LONGLONG


//-----------------------------------------------------------------------------
// specialized function for sending C strings to Inform object
//-----------------------------------------------------------------------------

inline Inform &operator<<(Inform &o, const char *s)
{
  o.stream() << s;
  return o;
}

//-----------------------------------------------------------------------------
// specialized function for sending strings to Inform object
//-----------------------------------------------------------------------------

inline Inform &operator<<(Inform &o, const std::string &s)
{
  o.stream() << s.c_str();
  return o;
}

//-----------------------------------------------------------------------------
// ostream_iterator for Inform output
//-----------------------------------------------------------------------------

#include <iterator>

template <class T>
class InformIterator 
{
public:
  typedef std::output_iterator_tag  iterator_category;
  typedef void   value_type;
  typedef void   difference_type;
  typedef void   pointer;
  typedef void   reference;

  InformIterator(Inform &s) : out_m(&s), delim_m(0) { }
  InformIterator(Inform &s, const char *d) : out_m(&s), delim_m(d) { }

  InformIterator &operator=(const T &value)
  {
    *out_m << value;
    if (delim_m != 0)
      *out_m << delim_m;
    return *this;
  }

  InformIterator &operator*()     { return *this; }
  InformIterator &operator++()    { return *this; }
  InformIterator &operator++(int) { return *this; }

private:
  Inform *out_m;
  const char *delim_m;
};

// } // namespace POOMA

//////////////////////////////////////////////////////////////////////

#endif // POOMA_UTILITIES_INFORM_H

// ACL:rcsinfo
// ----------------------------------------------------------------------
// $RCSfile: Inform.h,v $   $Author: richard $
// $Revision: 1.33 $   $Date: 2004/11/01 18:17:17 $
// ----------------------------------------------------------------------
// ACL:rcsinfo
