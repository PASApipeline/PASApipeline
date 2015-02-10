#include "Lobject.h"

Lobject::Lobject (int index, int num_alignments) {
  this->index = index;
  this->num_alignments = num_alignments;
  LscoreF = 0; 
  LscoreR = 0;
  toLptr = 0;
  fromLptr = 0;

  // init the containment array:
  contained_cdna_indices.resize(num_alignments);
  for (int i = 0; i < num_alignments; i++) {
    contained_cdna_indices[i] = false;
  }
  
}

void Lobject::setContainedIndices(vector<int> indices) {
  num_contained_indices = 0;
  
  for (int i=0; i < indices.size(); i++) {
    contained_cdna_indices[indices[i]] = true;
    
    LscoreF++;
    LscoreR++;
    num_contained_indices++;
  }
  
}



int Lobject::num_unique_contained (Lobject& other) {
  int num = 0;
  
  for (int i=0; i < num_alignments; i++) {
    if (contained_cdna_indices[i] && ! other.contained_cdna_indices[i]) {
      num++;
    }
  }
  
  return (num);
}


string Lobject::toString () {
  ostringstream os;
  
  os << "Lobject index: [" << index << "] has LscoreF: " << LscoreF 
     << ", LscoreR: " << LscoreR 
     << ", combinedScore: " << combined_score 
     << endl << "contains the following alignment indices:" << endl;
  
  for (int i = 0; i < num_alignments; i++) {
    if (contained_cdna_indices[i]) {
      os << "\tindex: " << index << " contains " << i << endl;
    }
  }
  return (os.str());
}

void Lobject::setTraceIndices(vector<int> traces) {
  traceIndices = traces;
}

vector<int> Lobject::getTraceIndices() {
  return (traceIndices);
}

