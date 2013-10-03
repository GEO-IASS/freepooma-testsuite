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
// ----------------------------------------------------------------------
// Delete test #1.
//-----------------------------------------------------------------------------

#include "Utilities/algorithms.h"
#include "Utilities/Tester.h"
#include "Pooma/Pooma.h"

#include <iostream>
#include <vector>
#include <ctime>
#include <set>
#include <algorithm>

#define SHIFTUP
#define BACKFILL

double timer();
int randN(int n);

typedef InformIterator<int> OutIter_t;

double timer()
{
  return double(std::clock()) / CLOCKS_PER_SEC;
}

void test(Pooma::Tester &t, 
          int num_elements, int num_trails, double delete_fraction, 
          bool check_results, bool print_results, bool delete_end);

int main(int argc, char *argv[])
{
  Pooma::initialize(argc, argv);

  Pooma::Tester tester(argc, argv);

  double start_main = timer();

  test(tester, 30,    1, .30, true,  true, false);
  test(tester, 30,    1, .30, true,  true, true);
  test(tester, 100,   1, .75, true,  true, false);
  test(tester, 100,   1, .75, true,  true, true);
  test(tester, 300, 100, .99, true, false, true);
  test(tester, 300, 100, .99, true, false, false);
  test(tester, 300, 100, .50, true, false, true);
  test(tester, 300, 100, .50, true, false, false);
  test(tester, 300, 100, .10, true, false, true);
  test(tester, 300, 100, .10, true, false, false);
  test(tester, 300, 100, .01, true, false, true);
  test(tester, 300, 100, .01, true, false, false);
  
  test(tester, 10000,    2,   .001, true, false, true);
  test(tester, 10000,    2,   .001, true, false, false);
  test(tester, 1000000,  1, .00001, true, false, true);
  test(tester, 1000000,  1, .00001, true, false, false);

  tester.out() << "Total run time = " << timer() - start_main 
               << " seconds." << std::endl;
  
  tester.out() << "All Done!" << std::endl;
  int res = tester.results("delete_test1 " );
  Pooma::finalize();

  return res;
}


typedef std::vector<int> KillList_t;
typedef std::vector<int> Array_t;

void delete_shiftup_orig (Array_t &data, const KillList_t &killList);
void delete_shiftup_test1(Array_t &data, const KillList_t &killList);
void delete_shiftup_test2(Array_t &data, const KillList_t &killList);
void delete_shiftup_test3(Array_t &data, const KillList_t &killList);

void delete_backfill_orig (Array_t &data, const KillList_t &killList);
void delete_backfill_test1(Array_t &data, const KillList_t &killList);
void delete_backfill_test2(Array_t &data, const KillList_t &killList);

// Generate a random number in the range [1,...,10]:

int randN(int n)
{
  return int(1 + n*(double(std::rand())/RAND_MAX));
}



void test(Pooma::Tester &t, 
          int num_elements, int num_trials, double delete_fraction, 
          bool check_results, bool print_results, bool delete_end)
{  
  const int  size     = num_elements;
  const int  trials   = num_trials;
  const int  inc_size = static_cast<int>(2.0 / delete_fraction);
  const bool verbose  = t.verbose() && print_results;
          
  using std::endl;
  
  KillList_t killList;
  killList.reserve(size);
  
  int i = randN(3) - 1; // First element in [0,1,2]
  
  while (i < size)
  {
    killList.push_back(i);
    i += randN(inc_size); // Increment by [1,...,inc_size]
  }
          
  if (delete_end)
    {
      // ensure that the endpoint is NOT deleted.
      if (killList[killList.size()-1] == size-1) 
        killList.resize(killList.size()-1);
    }
  else
    {
      // ensure that the endpoint IS deleted.
      if (killList[killList.size()-1] != size-1)
        killList.push_back(size-1);
    }
  
  t.out() << "Testing with " 
          << num_elements << " elements, " 
          << num_trials << " trials, "
          << "Deleting " << killList.size() << " elements." 
          << endl;

  OutIter_t outIter(t.out(), " ");

  if (verbose)
    {
      t.out() << "Kill List:" << endl;
      t.out() << "size   = " << killList.size() << endl;
      t.out() << "values = " ;
      std::copy(killList.begin(), killList.end(), outIter);
      t.out() << endl << endl;
    }
    
  Array_t data(size);
  
  for (int j = 0; j < size; ++j)
    data[j] = j;
    
  if (verbose)
    {
      t.out() << "Data:" << endl;
      t.out() << "size   = " << data.size() << endl;
      t.out() << "values = ";
      std::copy(data.begin(), data.end(), outIter);
      t.out() << endl << endl;
    }
  
  Array_t a = data;

  double start, finish, time;
  
  std::set<int> result_set;
  std::set<int> test_set;
  Array_t::iterator pos;
  Array_t::iterator end;
  typedef std::pair<std::set<int>::iterator,bool> Ret_t;

#ifdef BACKFILL
  
  time = 0;
  for (int j = 0; j < trials; ++j)
    {
      a = data;
      start = timer();
      delete_backfill_orig(a, killList);
      finish = timer();
      time += finish - start;
    }
    
  t.out() << "Time for original backfill = " << time 
          << " seconds." << endl;

  if (verbose)
    {
      t.out() << "After deleting with backfill:" << endl;
      t.out() << "size   = " << a.size() << endl;
      t.out() << "values = ";
      std::copy(a.begin(), a.end(), outIter);
      t.out() << endl << endl;
    }
    
  if (check_results)
    {
      pos = a.begin();
      end = a.end();
      while (pos != end)
        {
          Ret_t status = result_set.insert(static_cast<std::set<int>::key_type>(*pos++));
          t.check(status.second);
        }
    }
      
// Not try new algorithms:

  time = 0;
  for (int j = 0; j < trials; ++j)
    {
      a = data;
      start = timer();
      delete_backfill_test1(a, killList);
      finish = timer();
      time += finish - start;
    }

  t.out() << "Time for backfill test1 = " << time << " seconds." << endl;
  
  if (verbose)
    {
      t.out() << "After deleting with backfill:" << endl;
      t.out() << "size   = " << a.size() << endl;
      t.out() << "values = ";
      std::copy(a.begin(), a.end(), outIter);
      t.out() << endl << endl;
    }

  if (check_results)
    {
      pos = a.begin();
      end = a.end();
      while (pos != end) test_set.insert(*pos++);
      t.check(test_set == result_set);
      test_set.clear();
    }
      
  time = 0;
  for (int j = 0; j < trials; ++j)
    {
      a = data;
      start = timer();
      delete_backfill_test2(a, killList);
      finish = timer();
      time += finish - start;
    }

  t.out() << "Time for backfill test2 = " << time << " seconds." << endl;
  
  if (verbose)
    {
      t.out() << "After deleting with backfill:" << endl;
      t.out() << "size   = " << a.size() << endl;
      t.out() << "values = ";
      std::copy(a.begin(), a.end(), outIter);
      t.out() << endl << endl;
    }

  if (check_results)
    {
      pos = a.begin();
      end = a.end();
      while (pos != end) test_set.insert(*pos++);
      bool OK = test_set == result_set;
      if (!OK)
        {
          std::copy(test_set.begin(), test_set.end(), outIter);
          t.out() << endl;
          std::copy(result_set.begin(), result_set.end(), outIter);
          t.out() << "FAILED" << endl;
        }
      t.check(OK);
      test_set.clear();
    }
      
#endif

#ifdef SHIFTUP

  time = 0;
  for (int j = 0; j < trials; ++j)
    {
      a = data;
      start = timer();
      delete_shiftup_orig(a, killList);
      finish = timer();
      time += finish - start;
    }
    
  t.out() << "Time for original shiftup  = " << time
               << " seconds." << endl;

  if (verbose)
    {
      t.out() << "After deleting with shiftup:" << endl;
      t.out() << "size   = " << a.size() << endl;
      t.out() << "values = ";
      std::copy(a.begin(), a.end(), outIter);
      t.out() << endl << endl;
    }
    
  Array_t answer;
  if (check_results) answer = a;
  
// Not try new algorithms:

  time = 0;
  for (int j = 0; j < trials; ++j)
    {
      a = data;
      start = timer();
      delete_shiftup_test1(a, killList);
      finish = timer();
      time += finish - start;
    }

  t.out() << "Time for modified shiftup #1 = " << time 
               << " seconds." << endl;
  
  if (verbose)
    {
      t.out() << "After deleting with shiftup:" << endl;
      t.out() << "size   = " << a.size() << endl;
      t.out() << "values = ";
      std::copy(a.begin(), a.end(), outIter);
      t.out() << endl << endl;
    }

  if (check_results)
    t.check(a == answer);

  time = 0;
  for (int j = 0; j < trials; ++j)
    {
      a = data;
      start = timer();
      delete_shiftup_test2(a, killList);
      finish = timer();
      time += finish - start;
    }

  t.out() << "Time for modified shiftup #2 = " << time 
               << " seconds." << endl;
  
  if (verbose)
    {
      t.out() << "After deleting with shiftup:" << endl;
      t.out() << "size   = " << a.size() << endl;
      t.out() << "values = ";
      std::copy(a.begin(), a.end(), outIter);
      t.out() << endl << endl;
    }

  if (check_results)
    t.check(a == answer);
  
  time = 0;
  for (int j = 0; j < trials; ++j)
    {
      a = data;
      start = timer();
      delete_shiftup_test3(a, killList);
      finish = timer();
      time += finish - start;
    }

  t.out() << "Time for modified shiftup #3 = " << time 
               << " seconds." << endl;
  
  if (verbose)
    {
      t.out() << "After deleting with shiftup:" << endl;
      t.out() << "size   = " << a.size() << endl;
      t.out() << "values = ";
      std::copy(a.begin(), a.end(), outIter);
      t.out() << endl << endl;
    }

  if (check_results)
    t.check(a == answer);
  
#endif
} // test(...)

void delete_backfill_orig(Array_t &data, const KillList_t &killList)
{
  int killed = killList.size();
  int dom2 = data.size() - 1;
  int fill = dom2;
  int start = killList[0];
  int end   = killList[killed-1];

  // Now loop through the data and destroy the values

  int klIdx = killed - 1;
  while (start <= end)
    {
      if (end != fill)
        {
          data[end] = data[fill];
        }
      --fill;
      if (--klIdx < 0)
	break;
      end = killList[klIdx];
    }

  // Now resize the data block to the new size

  data.resize(data.size()-killed);
}

void delete_shiftup_orig(Array_t &data, const KillList_t &killList)
{
  int killed = killList.size();
  int dom2 = data.size() - 1;
  int start = killList[0];
  int end   = killList[killed-1];

  // Now loop through the data and destroy the values
  
  int insert = (-1);
  for (int i=0; i < killed; ++i)
    {
      int ks = killList[i];
      int kl = 1;
      while (i < (killed - 1) && (killList[i+1]) == (ks+kl))
	{
	  ++kl;
	  ++i;
	}

      int shift1 = ks + kl;
      if (shift1 <= dom2)
	{
	  int shift2 = dom2 + 1;
	  if (i < killed - 1)
	    shift2 = killList[i+1];

	  if (insert < 0)
	    insert = ks;
#if 1
	  while (shift1 < shift2)
	    data[insert++] = data[shift1++];
#else    
	  std::copy(data.begin()+shift1, data.begin()+shift2, 
	            data.begin()+insert);
	  insert += shift2-shift1;
	  shift1 += shift2-shift1;
#endif
	}
    }

  // Now resize the data block to the new domain size

  data.resize(data.size() - killed);
}

void delete_backfill_test1(Array_t &data, const KillList_t &killList)
{
  // Destroy parameters.  
  // All indices are offsets from the first element in the datablock.

  int killed = killList.size();
  
  int last        = data.size() - 1;  
  int firstKilled = killList[0];
  int lastKilled  = killList[killed-1]; 

  // Now loop through the data and destroy the desired elements.  
  // No data has to be moved if we're destoying values that are at
  // the end of the sequence. Thus we look for these first.
  
  // NOTE: This algorithm assumes that the killList is sorted.
  // (More specifically, it assumes that killList.first() is
  // the minimum element, that killList.last() is the maximal
  // element, and that if we do delete points at the end,
  // that all of these points are at the end of the killList
  // and sorted in order.)
  
  int killIndex = killed - 1;

  while (killIndex >= 0)
    {
    // The position in the underlying data block of the point 
    // to be deleted:
     
      int pos = killList[killIndex];
                  
    // If we're no longer deleting the last element, we break and
    // finish up below. 
        
      if (pos != last) break;

    // We're still on the last element. Deleting it is simply a matter
    // of moving our last pointer and continuing.
              
      --last;
      --killIndex;
    }

  // Now we're deleting non-last elements. This only executes if
  // the above loop was terminated early by the break.
  
  while (killIndex >= 0)
    {
    // Where we're deleting from:
        
      int pos = killList[killIndex];
          
    // For BackFill, we overwrite the deleted element with
    // the current last element and then decrement the "last"
    // index since that element has effectively been moved.
          
      data[pos] = data[last];
        
      --last;
      --killIndex;
    }

  data.resize(data.size() - killed);
}

void delete_shiftup_test1(Array_t &data, const KillList_t &killList)
{
  int killed = killList.size();
  int dom2 = data.size() - 1;
  int start = killList[0];
  int end   = killList[killed-1];

  // Now loop through the data and destroy the values
  
  int insert = killList[0];
  for (int i=0; i < killed; ++i)
    {
      int ks = killList[i];
      int kl = 1;
      while (i < (killed - 1) && (killList[i+1]) == (ks+kl))
	{
	  ++kl;
	  ++i;
	}

      int shift1 = ks + kl;
      if (shift1 <= dom2)
	{
	  int shift2 = dom2 + 1;
	  if (i < killed - 1)
	    shift2 = killList[i+1];

          int len = shift2 - shift1;
#if 1
	  std::copy(data.begin() + shift1, data.begin() + shift2, 
	            data.begin() + insert);
#else    
          std::memmove(data.begin()+insert, 
                       data.begin()+shift1, 
                       sizeof(int)*len);
#endif
	  insert += len;
	  shift1 += len;
	}
    }

  // Now resize the data block to the new domain size

  data.resize(data.size() - killed);
}


void delete_shiftup_test2(Array_t &data, const KillList_t &killList)
{
  const int killed = killList.size();
  const int size   = data.size();
  
  const std::vector<Array_t::value_type>::iterator data_begin = data.begin();
  
  // Now loop through the data and destroy the values
  
  int insert_index = killList[0];
  
      for (int i = 0; i < killed; ++i)
        {
          int copy_begin_index = killList[i] + 1;
          int inext = i + 1;
          while (inext < killed && copy_begin_index == killList[inext])
	    {
	      ++copy_begin_index;
	      ++inext;
	    }
	  i = inext - 1;

          if (copy_begin_index < size)
	    {
  	      int copy_end_index = size;
	      if (inext < killed)
	        copy_end_index = killList[inext];
	        
	      std::copy(data_begin + copy_begin_index, 
			data_begin + copy_end_index,
			data_begin + insert_index);
	      insert_index += copy_end_index - copy_begin_index;
	    }
        }
      
  // Now resize the data block to the new domain size

  data.resize(data.size() - killed);
}



void delete_shiftup_test3(Array_t &data, const KillList_t &killList)
{  
  // Use our generic algorithm to delete the elements:
  using Pooma::Algorithms::delete_shiftup;
  
  int killed = delete_shiftup(data.begin(), data.end(), 
			      killList.begin(), killList.end());

  // Now resize the data block to the new domain size:

  data.resize(data.size() - killed);
}

void delete_backfill_test2(Array_t &data, const KillList_t &killList)
{  
  // Use our generic algorithm to delete the elements:
  using Pooma::Algorithms::delete_backfill;

  int killed = delete_backfill(data.begin(), data.end(), 
			       killList.begin(), killList.end());

  // Now resize the data block to the new domain size:

  data.resize(data.size() - killed);
}
