#include <iostream>
#include <string>
#include <vector>
#include "alignment_segment.h"
#include "common_structs.h"

#ifndef __CDNA_ALIGNMENT__
#define __CDNA_ALIGNMENT__

using namespace std;

class CDNA_alignment {

 public: 
  CDNA_alignment(vector<Alignment_segment> seglist, char orient); //constructor
  
  
  void set_orientation (char);
  char get_orientation ();
  
  void set_title (string);
  string get_title ();
  
  void set_coords (int lend, int rend);
  struct coordset& get_coords ();
  
  void add_alignment_segment(Alignment_segment);
  vector<Alignment_segment>& get_alignment_segments();
  void delete_all_segments();

  string toString();
  string toAlignIllustration(int subtract, int rel_max, int lineLength);
  
  int num_segments;


 private:
  string title; // comname, header, accession, whatever you want to call it.
  char orient;
  vector<Alignment_segment> alignment_segs;
    
  struct coordset coords;
    
  void init(CDNA_alignment*);
  void refineAlignment();
  int coordConverter(int coord, int subtract, int rel_max, int linelength);
};

#endif
