
#include "common_structs.h"
#include "common_subs.h"
#include <string>

using namespace std;

#ifndef __ALIGNMENT_SEGMENT__
#define __ALIGNMENT_SEGMENT__

class Alignment_segment {
 public:
  
  Alignment_segment (int genomic_lend, int genomic_rend); //constructor
  Alignment_segment (struct coordset&); // copies coordset info over to local new coordset.

  struct coordset& get_coords();
  void set_coords(int lend, int rend);
  
  void set_type (string); 
  string get_type (); 
  
  void set_left_splice_junction(bool);
  bool get_left_splice_junction();
  
  void set_right_splice_junction(bool);
  bool get_right_splice_junction();

  string toString();
  string type; // first|last|internal|single|?
  
 private:
  void init();
  
  struct coordset coords; //genome coords of segment.
  bool has_left_splice_junction;
  bool has_right_splice_junction;
  
};

#endif

