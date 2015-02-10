#include "cdna_alignment.h"
#include <iostream>
#include <sstream>
#include <algorithm>

// constructor
CDNA_alignment::CDNA_alignment (vector<Alignment_segment> incomingSegs, char orient) {
  int num_segments = incomingSegs.size();
  if (num_segments <= 0) {
    cerr << "Empty list of alignment segments. Fatal." << endl;
    exit(1);
  }
  
  // constructor 
  init(this);
  this->orient = orient;
  this->alignment_segs = incomingSegs;
  this->refineAlignment();

}


void CDNA_alignment::init(CDNA_alignment* alignment) {
  this->title = "";
  this->coords.lend = 0;
  this->coords.rend = 0;

}


bool segmentSortCriteria (Alignment_segment a1, Alignment_segment a2) {
  if (a1.get_coords().lend <= a2.get_coords().lend) {
    return (true);
  } else {
    return (false);
  }
}


void CDNA_alignment::refineAlignment() {
  vector<Alignment_segment>& seglist = this->get_alignment_segments();
  
  int size = seglist.size();
  if (size > 1) {
    // sort them according to lend position:
    sort (seglist.begin(), seglist.end(), segmentSortCriteria);
  }
  
  int min_lend = seglist[0].get_coords().lend;
  int max_rend = seglist[size-1].get_coords().rend;
    
  this->set_coords(min_lend,max_rend);

  num_segments = seglist.size();
  // classify each segment according to type.
  if (num_segments == 1) {
    seglist[0].type = "single";
  } else {
    for (int i=0; i < num_segments; i++) {
      Alignment_segment& seg = seglist[i];
      if (i == 0) {
	seg.type = "first";
	seg.set_right_splice_junction(true);
      } else if (i == num_segments-1) {
	seg.type = "last";
	seg.set_left_splice_junction(true);
      } else {
	// must be internal
	seg.type = "internal";
	seg.set_left_splice_junction(true);
	seg.set_right_splice_junction(true);
      }
    }
  }
}

void CDNA_alignment::set_orientation (char orient) {
  this->orient = orient;
}

char CDNA_alignment::get_orientation () {
  return (this->orient);
}

void CDNA_alignment::set_title (string t) {
  this->title = t;
}

string CDNA_alignment::get_title () {
  return (this->title);
}

void CDNA_alignment::set_coords (int l, int r) {
  this->coords.lend = l;
  this->coords.rend = r;
}

struct coordset& CDNA_alignment::get_coords() {
  return (this->coords);
}

void CDNA_alignment::add_alignment_segment (Alignment_segment as) {
  this->alignment_segs.push_back(as);
}

vector<Alignment_segment>& CDNA_alignment::get_alignment_segments() {
  return (this->alignment_segs);
}

void CDNA_alignment::delete_all_segments () {
  this->alignment_segs.clear();
}

string CDNA_alignment::toString() {
  ostringstream os;
  os << get_title() << "," << orient << ",";
  for (int i = 0; i < alignment_segs.size(); i++) {
    os << alignment_segs[i].toString();
    if (i != alignment_segs.size() - 1) {
      os << ",";
    }
  }
  string output = os.str();
  return (output);
}

string CDNA_alignment::toAlignIllustration (int subtract, int rel_max, int lineLength) {
  char* token = new char[lineLength+1];
  for (int i=0; i < lineLength; i++) {
    token[i] = ' '; // initialize to whitespace
  }
  token[lineLength] = '\0'; // null char, terminator.
  char orient = get_orientation();
  int max_r_rel = 0;
  for (int i=0; i < alignment_segs.size(); i++) {
    Alignment_segment& seg = alignment_segs[i];
    struct coordset& coords = seg.get_coords();
    int lend = coords.lend;
    int rend = coords.rend;
    int l_rel = coordConverter(lend, subtract, rel_max, lineLength);
    int r_rel = coordConverter(rend, subtract, rel_max, lineLength);
    if (max_r_rel < r_rel && r_rel <= lineLength) {
      max_r_rel = r_rel;
    } 
    // init seg representation.
    for (int j=l_rel; j<= r_rel; j++) {
      token[j] = '-';
    }

    // add left splice indicator
    if (seg.get_left_splice_junction()) {
      token[l_rel] = '<';
    } else if (seg.type != "first" && seg.type != "single") {
      token[l_rel] = '|';
    }
    
    // add right splice indicator
    if (seg.get_right_splice_junction()) {
      token[r_rel] = '>';
    } else if (seg.type != "last" && seg.type != "single") {
      token[r_rel] = '|';
    }
  }
  token[max_r_rel+1] = '\0'; // truncate string after last char.
  
  string outline (token);
  char orientation [2];
  orientation[0]=orient;
  orientation[1]='\0';
  string myOrient (orientation);
  outline += " (" + myOrient + ") " + get_title();
  return (outline);
}

int CDNA_alignment::coordConverter (int coord, int subtract, int rel_max, int linelength) {
  return ((int) ( (coord-subtract)/(float)rel_max * linelength + 0.5));
}
    
