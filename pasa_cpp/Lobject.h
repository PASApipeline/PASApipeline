#ifndef __Lobject__
#define __Lobject__

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

extern bool DEBUG;

using namespace std;

class Lobject {
 public: 
  Lobject(int index, int num_alignments);
    
  int index;
  int num_alignments;
  
  vector<bool> contained_cdna_indices; 
  
  int LscoreF;
  int LscoreR;
  int combined_score;  /* LscoreF + LscoreR - num_contained_indices */
  Lobject* toLptr;
  Lobject* fromLptr;
  string toString();

  void setTraceIndices (vector<int>);
  vector<int> getTraceIndices();

  void setContainedIndices(vector<int>);
  int num_contained_indices; //includes self.

  int num_unique_contained(Lobject& other); /* returns number alignments contained in 
					       this, not in other */
 private:

  vector<int> traceIndices;
  
};

#endif
