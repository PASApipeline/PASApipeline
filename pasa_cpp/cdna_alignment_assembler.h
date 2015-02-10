#ifndef __CDNA_ALIGNMENT_ASSEMBLER__
#define __CDNA_ALIGNMENT_ASSEMBLER__

#include "Lobject.h"
#include "cdna_alignment.h"
#include "alignment_segment.h"
#include <string>
#include <vector>

using namespace std;

extern bool DEBUG;

class CDNA_alignment_assembler {
 public:
  CDNA_alignment_assembler(vector<CDNA_alignment>& incomingAlignments); /* constructor
									    the CDNA_alignment pointers
									    are copied over to a 	
									    local vector */

  ~CDNA_alignment_assembler(); /* Destructor. */  
  
  
  
  void assembleAlignments();
  
  vector<CDNA_alignment> get_assemblies ();  

  void set_fuzzlength(int);
  
  string toAlignIllustration(int lineLength); //default 100


 private:
  vector<CDNA_alignment>& alignments;  //after sorted by lend, must stay in sorted positions.
  vector<CDNA_alignment> assemblies; //assemblies added in order. 

  vector<vector<int> > assembly_containment_list; // holds list of alignment indices contained in assemblies.
  int fuzzlength; // default set to 20 bp
  
  bool canMerge(CDNA_alignment&, CDNA_alignment&);
  CDNA_alignment mergeAlignments (CDNA_alignment& A, CDNA_alignment& B); 
  
  vector<Lobject> Lobjects;

  void do_full_Fscan();
  void do_full_Rscan();

  bool encapsulates(CDNA_alignment& A, CDNA_alignment& B); //returns true if alignment A encapuslates the span of alignment B
  void determine_compatibilities_and_encapsulations();
  vector<int> forwardTrace (int); // index to begin trace.
  vector<int> backTrace(int); 
  CDNA_alignment create_assembly(vector<int>); 
  vector<int> get_top_scoring_alignment();
  vector<int> get_alignment_assembly_nucleating_at_alignment_index(int index);
  vector<int> unique_entries(vector<vector<int> >);
  void populateLobjects();
  
  bool** compatibilities;
  bool**  encapsulations;
  int num_alignments;
  
  Lobject* get_max_missing_Lobj(vector<Lobject*>&, map<int,bool>&);
  
};

#endif


