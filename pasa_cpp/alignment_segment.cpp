#include "alignment_segment.h"
#include <sstream>

Alignment_segment::Alignment_segment (int genomic_lend, int genomic_rend) {
  if (genomic_lend > genomic_rend) {
    swap_ints(genomic_lend, genomic_rend);
  }
  this->init();
  this->coords.lend = genomic_lend;
  this->coords.rend = genomic_rend;
  
}


Alignment_segment::Alignment_segment(struct coordset& coords) {
  init();
  set_coords(coords.lend, coords.rend);
}


void Alignment_segment::init() {
  this->type = "?";
  this->has_left_splice_junction = false;
  this->has_right_splice_junction = false;
}


struct coordset& Alignment_segment::get_coords() {
  return (this->coords);
}
  
void Alignment_segment::set_coords (int lend, int rend) {
  coords.lend = lend;
  coords.rend = rend;
}

void Alignment_segment::set_type (string t) {
  type = t;
}

string Alignment_segment::get_type () {
  return (type);
}

void Alignment_segment::set_left_splice_junction(bool val) {
  has_left_splice_junction = val;
}

bool Alignment_segment::get_left_splice_junction() {
  return (has_left_splice_junction);
}

void Alignment_segment::set_right_splice_junction(bool val) {
  has_right_splice_junction = val;
}

bool Alignment_segment::get_right_splice_junction() {
  return (has_right_splice_junction);
}

string Alignment_segment::toString() {
  ostringstream os;
  os <<  coords.lend << "-" << coords.rend;
  string output = os.str();
  return (output);
}

