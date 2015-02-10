#include "cdna_alignment_assembler.h"
#include <algorithm>
#include <map>
#include <iostream>
#include <sstream>

bool sort_CDNA_alignments (CDNA_alignment a, CDNA_alignment b) {
  
  struct coordset& acoords = a.get_coords();
  struct coordset& bcoords = b.get_coords();
  int a_lend = acoords.lend;
  int b_lend = bcoords.lend;
  
  if (a_lend < b_lend) {
    return (true);
  } else {
    return (false);
  }
}


CDNA_alignment_assembler::CDNA_alignment_assembler (vector<CDNA_alignment>& incomingAlignments) : alignments (incomingAlignments) {
  
  // sort alignments by lend position:
  if (DEBUG) {
    cout << "-sorting alignments by lend position." << endl;
  }
  
  sort (alignments.begin(), alignments.end(), sort_CDNA_alignments);
  if (DEBUG) {
    cout << "-done sorting alignments." << endl;
    cout << "Here's the alignments orderd by lend pos." << endl;
    for (int i=0; i < alignments.size(); i++) {
      cout << i << " " << alignments[i].toString() << endl;
    }
  }
  
  fuzzlength = 20; //default setting.
  
  num_alignments = incomingAlignments.size();
  
  // initialize pointers.
  compatibilities = 0;
  encapsulations = 0; 
  
  
}

// Destructor.
CDNA_alignment_assembler::~CDNA_alignment_assembler () {
  
  if (compatibilities != 0) {
    free2Darray(compatibilities);
  }
  if (encapsulations != 0) {
    free2Darray(encapsulations);
  }
}


bool sort_Lobjects_via_combined_Lscore_R_F (Lobject* a, Lobject* b) {
  
  if (a->combined_score < b->combined_score) {
    return(true);
  } else {
    return (false);
  }
}



void CDNA_alignment_assembler::assembleAlignments() {
  if (DEBUG) { cout << "instantiating required compatibility and encapsulation 2D-arrays." << endl; }
  compatibilities = twoDarray(num_alignments, num_alignments, false);
  encapsulations = twoDarray(num_alignments, num_alignments, false);
  
  if (DEBUG) {
    cout << "Assembling alignments." << endl;
    cout << "Determining compatibilities and encapsulations." << endl;
  }

  determine_compatibilities_and_encapsulations();
  
  if (DEBUG) { cout << "Populating Lobjects." << endl; }
  
  populateLobjects();
  
  if (DEBUG) { cout << "Doing full Fscan." << endl; }
  
  do_full_Fscan();
  
  if (DEBUG) { cout << "Getting maximum scoring alignment assembly." << endl; }
  
  vector<int> topAssemblyIndices = get_top_scoring_alignment();
  CDNA_alignment assembly = create_assembly(topAssemblyIndices);
  
  ostringstream assemblyTitle;
  assemblyTitle << "assembly_" << assemblies.size();
  assembly.set_title(assemblyTitle.str());
  assemblyTitle.str(""); // clear it for later.
  assemblies.push_back(assembly);
  assembly_containment_list.push_back(topAssemblyIndices);
  
  if (topAssemblyIndices.size() == num_alignments) {
    return; // all alignments assembled.
  }
  
  // find all maximal assemblies for alignments not included
  // in the most maximal assembly
  // In pasa2, we do this from scratch, ignoring our best assembly.  We'll refind it and better deal with tie situations.
  
  /* reinit */
  assemblies.clear();
  assembly_containment_list.clear();
  

  // Do Rscan so the Fscan,Rscan trace will provide maximal assemblies from specified starting indices.
  if (DEBUG) { cout << "Doing full Rscan." << endl; }
  do_full_Rscan();
  
  // get list of missing alignments:
  map<int,bool> accountedFor;
  // initialize
  for (int i=0; i < num_alignments; i++) {
    accountedFor[i] =false;
  }
  
  /* pasa2 don't do it
  // track already accounted for alignments.
  for (int i=0; i < topAssemblyIndices.size(); i++) {
  accountedFor[topAssemblyIndices[i]] = true;
  }
  */
  
  // gather alignments not accounted for
  vector<Lobject*> untraversedLobjs;
  for (int i=0; i < num_alignments; i++) {
    if (!accountedFor[i]) {
      // in pasa2, this is all of 'em.
      untraversedLobjs.push_back(&Lobjects[i]);
    }
  }
  
    
  // assign combined scores
  for (int i = untraversedLobjs.size() - 1; i >= 0; i--) {
    Lobject& nextBestLobj = *(untraversedLobjs[i]);
    int nucleatingIndex = nextBestLobj.index;
    int combined_score = nextBestLobj.LscoreF + nextBestLobj.LscoreR - nextBestLobj.num_contained_indices;
  
    if (DEBUG) { 
      cout << "Combined Lscore for index: " << nucleatingIndex << " is " 
           << combined_score 
           <<" (LscoreF: " << nextBestLobj.LscoreF 
           << ", LscoreR: " << nextBestLobj.LscoreR 
           << ")" <<  endl; 
    }
    
    nextBestLobj.combined_score = combined_score;

    vector<int> alignmentIndices = get_alignment_assembly_nucleating_at_alignment_index(nucleatingIndex);
    nextBestLobj.setTraceIndices(alignmentIndices);
    

  }
  
  
  // sorted in order of increasing combined scores.
  // get alignments in order of decreasing combined score until
  // all alignments are accounted for.

  sort (untraversedLobjs.begin(), untraversedLobjs.end(), sort_Lobjects_via_combined_Lscore_R_F);

  
  // first, lets bin those assemblies that have the same number of elements and order them by 
  // the number of missing alignments, desc.
  
  vector< vector<Lobject*> > lobjsSameSizeVec;
  int max_index = untraversedLobjs.size() - 1;
  int curr_score = untraversedLobjs[ max_index ]->combined_score;
  vector<Lobject*> curr_bin;
  
  // remember that untraversedLobjs in sorted asc
  if (DEBUG) { cout << "\n\n\n**** Binning lobjs with same combined scores.\n";}
  
  for (int i = untraversedLobjs.size() - 1; i >= 0; i--) {
    Lobject* lobj = untraversedLobjs[i];
    if (DEBUG) { cout << "-adding lobj[" << i << "] to bin. " << endl; }
    int score = lobj->combined_score;
    
    if (DEBUG) { cout << "curr_score: " << curr_score << ", lobj score: " << score << endl; }
    
    if (score == curr_score) {
      if (DEBUG) { cout << "\tscores same, adding to current bin." <<endl; }
      curr_bin.push_back(lobj);
    }
    else {
      if (DEBUG) { cout << "\tscores diff, archiving current bin, creating new bin with this lobj." << endl; }
      lobjsSameSizeVec.push_back(curr_bin);
      curr_bin.clear();
      curr_bin.push_back(lobj);
      curr_score = score;
    }
  }
  
  if (curr_bin.size() > 0) {
    if (DEBUG) { cout << "-pushing last lobj bin onto the stack." << endl; }
    lobjsSameSizeVec.push_back(curr_bin);
    curr_bin.clear();
  }
  
  
  if (DEBUG) {
    cout << "\n-describing binned Lobjs:\n\n";
    cout << "There are " << lobjsSameSizeVec.size() << " bins of Lobjs.\n\n";
    for (int i = 0; i < lobjsSameSizeVec.size(); i++) {
      vector<Lobject*> lobj_bin = lobjsSameSizeVec[i];
      cout << "\tbin: " << i << " of size: " << lobj_bin.size() << endl; 
      for (int j=0; j < lobj_bin.size(); j++) {
        cout << "\t\tLobj: " << lobj_bin[j]->index << endl; }
    }
  }
  
  // now lobjsSameSizeVec is sorted desc
  for (int i = 0; i < lobjsSameSizeVec.size(); i++) {
    vector<Lobject*> lobj_bin = lobjsSameSizeVec[i];
    
    if (DEBUG) { cout << "-**** Analyzing lobj_bin at bin pos: " << i << endl; }
        
        
    while (lobj_bin.size() > 0) {
      Lobject* max_missing_Lobj = get_max_missing_Lobj(lobj_bin, accountedFor);
      
      if (max_missing_Lobj == NULL) {
        // no more missing alignments in this bin.  go to next bin.
        if (DEBUG) { cout << "no more missing alignments in bin " << i << ". Trying next bin." << endl << endl; }
        break;
      }
      
      int nucleatingIndex = max_missing_Lobj->index;
      vector<int> alignmentIndices = max_missing_Lobj->getTraceIndices();
      // see if there are any unconsumed alignments in this assembly
      bool hasUnconsumed = false;
      for (int j=0; j < alignmentIndices.size(); j++) {
        int index = alignmentIndices[j];
        if (! accountedFor[index]) {
          hasUnconsumed = true;
          break;
        }
      }
      if (hasUnconsumed) {
        CDNA_alignment newAssembly = create_assembly(alignmentIndices);
        assemblyTitle << "assembly_" << assemblies.size();
        newAssembly.set_title(assemblyTitle.str());
        assemblyTitle.str("");
        assemblies.push_back(newAssembly);
        assembly_containment_list.push_back(alignmentIndices);
        
        // track newly consumed alignments
        for (int j=0; j<alignmentIndices.size(); j++) {
          int index = alignmentIndices[j];
          accountedFor[index] = true;
        }
        
        bool all_accounted_for = true;
        for (int j=0; j < num_alignments; j++) {
          if (! accountedFor[j]) {
            all_accounted_for = false;
            break;
          }
        }
        if (all_accounted_for) {
          // done.
          return;
        }
      } else {
        // lacks unconsumed alignment
        if (DEBUG) { cout << "Apparently, no unconsumed alignments in assembly nucleating at index: " << nucleatingIndex << endl; }
        // so none of the others in this bin will have an unconsumed alignment either
        break; // break inner
      }

      if (DEBUG) { 
        cout << "Now removing max Lobj from current bin\n";
        cout << "Before bin size: " << lobj_bin.size() << endl;
      }
      
      vector<Lobject*> new_lobj_bin; // store the rest that haven't been examined yet
      for (int i=0; i < lobj_bin.size(); i++) {
        Lobject* lobj = lobj_bin[i];
        if (lobj != max_missing_Lobj) {
          new_lobj_bin.push_back(lobj);
        }
      }
      
      if (DEBUG) { cout << "doing replacment.  new bin size: " << new_lobj_bin.size()  << ", old bin size: " << lobj_bin.size() << endl; }
      lobj_bin = new_lobj_bin;
      if (DEBUG) { cout << "after replacement, bin size: " << lobj_bin.size() << endl; }
            
    }  // end while

  } // end foreach bin
  
  // if you've gotten this far, you must not have found assemblies for 
  // all alignments.  Death ensues.
  cerr << "Not all alignments were accounted for by alignment assemblies." << endl;
  exit(5);
}


vector<CDNA_alignment> CDNA_alignment_assembler::get_assemblies () {
  return (assemblies);
}

void CDNA_alignment_assembler:: set_fuzzlength(int length) {
  fuzzlength = length;
}

bool CDNA_alignment_assembler::canMerge(CDNA_alignment& a1, CDNA_alignment& a2) {
  // check for overlap between alignments:
  if (! overlap(a1.get_coords(), a2.get_coords())) {
    if (DEBUG) { cout << "-can't merge: alignment coordsets don't overlap." << endl; }
    return (false);
  }
  
  // make sure orientation is equivalent:
  if (a1.get_orientation() != a2.get_orientation()) {
    if (DEBUG) { cout << "-can't merge: diff orientations." << endl; }
    return (false);
  }
  
  // check overlapping segments to ensure non-conflicting segments
  vector<Alignment_segment>& a1_segments = a1.get_alignment_segments();
  vector<Alignment_segment>& a2_segments = a2.get_alignment_segments();
  
  // align segment orders between a1 and a2
  int starting_a1 = -1;
  int starting_a2 = -1;
  for (int i=0; i < a1_segments.size(); i++) {
    Alignment_segment& a1_seg = a1_segments[i];
    for (int j=0; j < a2_segments.size(); j++) {
      Alignment_segment& a2_seg = a2_segments[j];
      if (overlap(a1_seg.get_coords(), a2_seg.get_coords())) {
        starting_a1 = i;
        starting_a2 = j;
        break;
      }
    }
    if (starting_a1 != -1 && starting_a2 != -1) {
      break;
    }
  }
  
  if (starting_a1 == -1 || starting_a2 == -1) {
    // couldn't align two segments of alignments a1 and a2
    if (DEBUG) { cout << "can't merge: couldn't align two segments of the alignments." << endl; }
    return(false);
  }
  
  if (! (starting_a1 == 0 || starting_a2 == 0)) {
    // couldn't map first segment of either to the other
    if (DEBUG) { cout << "can't merge: couldn't map first segment of either to the other." << endl; }
    return (false);
  }
  
  // check for compatible introns/exons within overlapping segments of alignments
  while (starting_a1 < a1_segments.size() && starting_a2 < a2_segments.size()) {
    Alignment_segment& a1_seg = a1_segments[starting_a1];
    Alignment_segment& a2_seg = a2_segments[starting_a2];
    struct coordset& a1_seg_coords = a1_seg.get_coords();
    struct coordset& a2_seg_coords = a2_seg.get_coords();
    int a1_lend = a1_seg_coords.lend;
    int a1_rend = a1_seg_coords.rend;
    int a2_lend = a2_seg_coords.lend;
    int a2_rend = a2_seg_coords.rend;
    
    if (overlap(a1_seg_coords, a2_seg_coords)) {
      
      //Analyze left splice junction
      
      // see if have identical splice sites
      if (a1_seg.get_left_splice_junction() || a2_seg.get_left_splice_junction()) {
        if (a1_seg.get_left_splice_junction() && a2_seg.get_left_splice_junction() && a1_lend != a2_lend) {
          // have different splice sites
          if (DEBUG) { cout << "can't merge: diff left splice sites." << endl; }
          return (false);
        } else if (a1_seg.get_left_splice_junction() && (a2_lend + fuzzlength < a1_lend)) {
          // not within fuzzdist
          if (DEBUG) { cout << "can't merge: left splice analysis, not within fuzz distance." << endl; }
          return (false);
        } else if (a2_seg.get_left_splice_junction() && (a1_lend + fuzzlength < a2_lend)) {
          // not within fuzzdist
          if (DEBUG) { cout << "can't merge: left splice analysis, not within fuzz distance." << endl;}
          return (false);
        }
      }
      
      // Analyze right splice junction
      if (a1_seg.get_right_splice_junction() || a2_seg.get_right_splice_junction()) {
        
        // see if identical splice sites:
        if (a1_seg.get_right_splice_junction() && a2_seg.get_right_splice_junction() && a1_rend != a2_rend) {
          // diff splice sites
          if (DEBUG) { cout << "can't merge: diff right splice sites." << endl; }
          return (false);
        } else if (a1_seg.get_right_splice_junction() && (a2_rend - fuzzlength > a1_rend)) {
          if (DEBUG) { cout << "can't merge: right splice analysis, not within fuzz distance." << endl; }
          return (false);
          // not within fuzzlength
        } else if (a2_seg.get_right_splice_junction() && (a1_rend - fuzzlength > a2_rend)) {
          // not within fuzzlength
          if (DEBUG) { cout << "can't merge: right splice analysis, not within fuzz distance." << endl; }
          return (false);
        }
      }
    } else { // no overlap
      // two ordered segments do not overlap each other.
      if (DEBUG) { cout << "can't merge: Two ordered segments do not overlap each other." << endl; }
      return (false);
    }
    starting_a1++;
    starting_a2++;
  }
  
  // passed all tests
  if (DEBUG) { cout << "-Merge possible. Passed all tests." << endl; }
  return (true);
}



CDNA_alignment CDNA_alignment_assembler::mergeAlignments(CDNA_alignment& A, CDNA_alignment& B) {
  map<int,bool> leftsplicecoords;
  map<int,bool> rightsplicecoords;
  
  char orientation = A.get_orientation();  // A and B should have identical orientations.
  
  // examine a1 coordinates
  vector<Alignment_segment>& a1_segments = A.get_alignment_segments();
  for (int i=0; i < a1_segments.size(); i++) {
    Alignment_segment& a1_seg = a1_segments[i]; 
    struct coordset& a1_coordset = a1_seg.get_coords();
    int lend = a1_coordset.lend;
    int rend = a1_coordset.rend;
    if (a1_seg.get_left_splice_junction()) {
      leftsplicecoords[lend] = true;
    }
    if (a1_seg.get_right_splice_junction()) {
      rightsplicecoords[rend] = true;
    }
  }
  
  // examine a2 coordinates
  vector<Alignment_segment>& a2_segments = B.get_alignment_segments();
  for (int i=0; i < a2_segments.size(); i++) {
    Alignment_segment& a2_seg = a2_segments[i]; 
    struct coordset& a2_coordset = a2_seg.get_coords();
    int lend = a2_coordset.lend;
    int rend = a2_coordset.rend;
    if (a2_seg.get_left_splice_junction()) {
      leftsplicecoords[lend] = true;
    }
    if (a2_seg.get_right_splice_junction()) {
      rightsplicecoords[rend] = true;
    }
  }
  
  vector<struct coordset> merged_coords;
  
  for (int i=0; i < a1_segments.size(); i++) {
    Alignment_segment& a1_seg = a1_segments[i];
    struct coordset& a1_coordset = a1_seg.get_coords();
    int a1_lend = a1_coordset.lend;
    int a1_rend = a1_coordset.rend;
    int merged_lend = -1;
    int merged_rend = -1;
    for (int j=0; j < a2_segments.size(); j++) {
      Alignment_segment& a2_seg = a2_segments[j];
      struct coordset& a2_coordset = a2_seg.get_coords();
      int a2_lend = a2_coordset.lend;
      int a2_rend = a2_coordset.rend;
      if (overlap(a1_coordset, a2_coordset)) {
        
        // determine merged_lend
        if (leftsplicecoords[a1_lend]) {
          merged_lend = a1_lend;
        } else if (leftsplicecoords[a2_lend]) {
          merged_lend = a2_lend;
        } else {
          merged_lend = min(a1_lend, a2_lend);
        }
        
        // determine merged_rend
        if (rightsplicecoords[a1_rend]) {
          merged_rend = a1_rend;
        } else if (rightsplicecoords[a2_rend]) {
          merged_rend = a2_rend;
        } else {
          merged_rend = max(a1_rend, a2_rend);
        }
        break;
      }
    }
    struct coordset merged_coordset;
    if (merged_lend != -1 && merged_rend != -1) {
      // adding overlapped coords
      merged_coordset.lend = merged_lend;
      merged_coordset.rend = merged_rend;
    } else {
      // must not have been any overlap.  Keep the a1 coords.
      merged_coordset.lend = a1_lend;
      merged_coordset.rend = a1_rend;
    }
    merged_coords.push_back(merged_coordset);
  }
  
  // add the unconsumed a2 coordsets
  for (int i=0; i < a2_segments.size(); i++) {
    Alignment_segment& a2_seg = a2_segments[i];
    struct coordset& a2_coords = a2_seg.get_coords();
    int a2_lend = a2_coords.lend;
    int a2_rend = a2_coords.rend;
    bool overlapFlag = false;
    for (int j=0; j < merged_coords.size(); j++) {
      struct coordset& m_coords = merged_coords[j];
      if (overlap(a2_coords, m_coords)) {
        overlapFlag = true;
        break;
      }
    }
    if (! overlapFlag) {
      // consuming a2 non-overlapping coordset
      struct coordset m_coords;
      m_coords.lend = a2_lend;
      m_coords.rend = a2_rend;
      merged_coords.push_back(m_coords);
    }
  }
  
  vector<Alignment_segment> new_seg_list;
  for (int i=0; i < merged_coords.size(); i++) {
    struct coordset& coords = merged_coords[i];
    Alignment_segment new_seg (coords);
    new_seg_list.push_back(new_seg);
  }
  
  CDNA_alignment merged_alignment (new_seg_list, orientation);
  
  return (merged_alignment);
  
}


void CDNA_alignment_assembler::do_full_Fscan() {
  
  for (int i=1; i < num_alignments; i++) {
    //must compare to previous alignments
    Lobject& Lobj = Lobjects[i];
    int top_score = 0;
    int top_scoring_index = -1;
    CDNA_alignment& i_alignment = alignments[i];
    
    for (int j = i-1; j >= 0; j--) {
      
      if (DEBUG) { cout << "FSCAN: comparing alignment " << i << " to " << j << endl; }
      
//      CDNA_alignment& j_alignment = alignments[j];
      Lobject& prevLobj = Lobjects[j];
      
      bool compatible = compatibilities[i][j];
      
      bool containment =  encapsulations[i][j] || encapsulations[j][i];
      
      if (DEBUG) { cout << "\tcompat: " << compatible << ", contained: " << containment << endl; }
      
      if (compatible && !containment) {
//        int curr_Lscore = prevLobj.LscoreF;
//        int Cscore = Lobj.num_unique_contained(prevLobj);
//        int curr_total_score = curr_Lscore + Cscore;
        int curr_total_score = prevLobj.LscoreF + Lobj.num_unique_contained(prevLobj);
        if (DEBUG) { cout << "FSCAN(" << i << "," << j << ") \tcurr total score: " << curr_total_score << endl;}
        
        if (curr_total_score > top_score) {
          top_scoring_index = j;
          top_score = curr_total_score;
        }
      }
    }
    if (top_scoring_index > -1) {
      Lobject& topLobj = Lobjects[top_scoring_index];
      Lobj.fromLptr = &topLobj;
      Lobj.LscoreF = top_score;
      if (DEBUG) { cout << "-Assigning fromLptr from index " << i << " -> " << top_scoring_index << endl; }
    }
  }
}



void CDNA_alignment_assembler::do_full_Rscan() {
  for (int i= num_alignments - 2; i >= 0; i--) {
    //must compare to previous alignments
    Lobject& Lobj = Lobjects[i];
    int top_score = 0;
    int top_scoring_index = -1;
    CDNA_alignment& i_alignment = alignments[i];
    
    for (int j = i+1; j < num_alignments; j++) {
      if (DEBUG) { cout << "RSCAN: comparing alignment " << i << " to " << j << endl; }
//      CDNA_alignment& j_alignment = alignments[j];
      Lobject& nextLobj = Lobjects[j];
      
      bool compatible = compatibilities[i][j];
      
      bool containment =  encapsulations[i][j] || encapsulations[j][i];
      if (DEBUG) { cout << "\tcompat: " << compatible << ", contained: " << containment << endl; }
      
      if (compatible && !containment) {
//        int curr_Lscore = nextLobj.LscoreR;
//        int Cscore = Lobj.num_unique_contained(nextLobj);
//        int curr_total_score = curr_Lscore + Cscore;
        int curr_total_score = nextLobj.LscoreR + Lobj.num_unique_contained(nextLobj);
        if (DEBUG) { cout << "RSCAN(" << i << "," << j << ") \tcurr total score: " << curr_total_score << endl;}
        
        if (curr_total_score > top_score) {
          top_scoring_index = j;
          top_score = curr_total_score;
        }
      }
    }
    if (top_scoring_index > -1) {
      Lobject& topLobj = Lobjects[top_scoring_index];
      Lobj.toLptr = &topLobj;
      Lobj.LscoreR = top_score;
      if (DEBUG) { cout << "-Assigning toLptr from index " << i << " -> " << top_scoring_index << " with top score " << top_score << endl; }
    }
  }
}

bool CDNA_alignment_assembler::encapsulates (CDNA_alignment& A, CDNA_alignment& B) {
  
  // checks to see if B span is contained within A span
  
  struct coordset& Acoords = A.get_coords();
  int a_lend = Acoords.lend;
  int a_rend = Acoords.rend;
  
  struct coordset& Bcoords = B.get_coords();
  int b_lend = Bcoords.lend;
  int b_rend = Bcoords.rend;
  
  if (b_lend >= a_lend && b_rend <= a_rend) {
    return (true);
  } else {
    return (false);
  }
  
}

void CDNA_alignment_assembler::determine_compatibilities_and_encapsulations() {
  
  // All vs. All comparison of alignments
  for (int i=0; i < num_alignments; i++) {
        
    for (int j=i+1; j < num_alignments; j++) {
    
      if (DEBUG) { cout << "can merge " << i << " to " << j << " ?" << endl; }
      bool mergeable = false;
      if (canMerge(alignments[i], alignments[j])) {
        compatibilities[i][j] = true;
        compatibilities[j][i] = true;
        mergeable = true;
      }
      if (mergeable) { 
        // check for encapsulation
        if (encapsulates(alignments[i],alignments[j])) {
          if (DEBUG) { cout << "alignment " << i << " encapsulates " << j << endl; }
          encapsulations[i][j] = true;
        }
        if (encapsulates(alignments[j],alignments[i])) {
          if (DEBUG) { cout << "alignment " << j << " encapsulates " << i << endl; }
          encapsulations[j][i] = true;
        }
        
      }
    }
  }
  
}


void CDNA_alignment_assembler::populateLobjects () {
  for (int i=0; i < num_alignments; i++) {
    Lobject L (i, num_alignments);
    vector<int> contained;
    contained.push_back(i); // include itself in the containment list.
    for (int j=0; j < num_alignments; j++) {
      if (encapsulations[i][j]) {
        if (DEBUG) { cout << "Populating Lobj: align " << i << " contains " << j << endl; }
        contained.push_back(j);
      }
    }
    L.setContainedIndices(contained);
    Lobjects.push_back(L);
  }
}



vector<int> CDNA_alignment_assembler::forwardTrace (int startIndex) {
  if (DEBUG) { cout << "Beginning forwardTrace, starting at index: " << startIndex << endl; }
  map<int,bool> tracker;
  Lobject* Lobj = & Lobjects[startIndex];
  while (Lobj != 0) {
    if (DEBUG) { cout << "trace index: " << Lobj->index << endl; }
    if (DEBUG) { cout << Lobj->toString(); }
    
    for (int i = 0; i < num_alignments; i++) {
      if (Lobj->contained_cdna_indices[i]) {
        tracker[i] = true;
      }
    }
    Lobj = Lobj->toLptr;
  }
  
  vector<int> unique;
  map<int,bool>::iterator pos;
  for (pos = tracker.begin(); pos != tracker.end(); pos++) {
    int index = pos->first;
    bool is_contained = pos->second;
    if (is_contained) {
      unique.push_back(index);
    }
  }
  return (unique);
}

vector<int> CDNA_alignment_assembler::backTrace(int startIndex) {
  if (DEBUG) { cout << "Beginning backTrace, starting at index: " << startIndex << endl; }
  map<int,bool> tracker;
  Lobject* Lobj = & Lobjects[startIndex];
  while (Lobj != 0) {
    if (DEBUG) { cout << "trace index: " << Lobj->index << endl;
                 cout << Lobj->toString(); }
    
    for (int i = 0; i < num_alignments; i++) {
      if (Lobj->contained_cdna_indices[i]) {
        tracker[i] = true;
      }
    }
    Lobj = Lobj->fromLptr;
  }
  
  vector<int> unique;
  map<int,bool>::iterator pos;
  for (pos = tracker.begin(); pos != tracker.end(); pos++) {
    int index = pos->first;
    bool is_contained = pos->second;
    if (is_contained) {
      unique.push_back(index);
    }
  }
  return (unique);
}

CDNA_alignment CDNA_alignment_assembler::create_assembly(vector<int> Alignment_index_listing) {
  // assemble alignments in order from left to right.
  if (Alignment_index_listing.empty()) {
    cerr << "empty list of indices, can't create assembly." << endl;
    exit(6);
  }
  
  sort(Alignment_index_listing.begin(), Alignment_index_listing.end());
  
  int alignment_index = Alignment_index_listing[0];
  
  CDNA_alignment assembly = alignments[alignment_index]; // copy the first alignment
  for (int i = 1; i < Alignment_index_listing.size(); i++) {
    alignment_index = Alignment_index_listing[i];
    CDNA_alignment& nextAlignment = alignments[alignment_index];
    CDNA_alignment newAssembly = mergeAlignments(assembly, nextAlignment);
    assembly = newAssembly;
  }
  
  return (assembly);
}


vector<int> CDNA_alignment_assembler::get_top_scoring_alignment() {
  // find the highest Lscore from Forward scan
  int top_score = 0;
  int top_scoring_index = -1;
  for (int i=0; i < Lobjects.size(); i++) {
    Lobject* L = & Lobjects[i];
    int Lscore = L->LscoreF;
    if (DEBUG) { cout << "LscoreF of alignment " << i << " is " << Lscore << endl; }
    if (Lscore > top_score) {
      top_scoring_index = i;
      top_score = Lscore;
    }
  }
  if (DEBUG) { cout << "Top score is alignment " << top_scoring_index << " with LscoreF of " << top_score << endl; }
  vector<int> indexListing = backTrace(top_scoring_index);
  return (indexListing);
}

vector<int> CDNA_alignment_assembler::get_alignment_assembly_nucleating_at_alignment_index(int index) {
  vector<vector<int> > vecvec;
  vecvec.push_back(backTrace(index));
  vecvec.push_back(forwardTrace(index));
  return (unique_entries(vecvec));
}

vector<int> CDNA_alignment_assembler::unique_entries(vector<vector<int> > vecvec) {
  // requires implementation
  map<int,bool> uniqueMap;
  for (int j=0; j < vecvec.size(); j++) {
    vector<int> myvec = vecvec[j];
    for (int i=0; i <myvec.size(); i++) {
      int entry = myvec[i];
      uniqueMap[entry] = true;
    }
  }
  
  vector<int> uniqueEntries;
  map<int,bool>::iterator pos;
  for (pos = uniqueMap.begin(); pos != uniqueMap.end(); pos++) {
    int entry = pos->first;
    uniqueEntries.push_back(entry);
  }
  
  return(uniqueEntries);
}

string CDNA_alignment_assembler::toAlignIllustration (int lineLength) {
  
  // get range for all alignment coordinates.
  int min_coord;
  int max_coord;
  vector<int> allCoords;
  for (int i=0; i < alignments.size(); i++) {
    struct coordset& coords = alignments[i].get_coords();
    allCoords.push_back(coords.lend);
    allCoords.push_back(coords.rend);
  }
  sort (allCoords.begin(), allCoords.end());
  min_coord = allCoords[0];
  max_coord = allCoords[allCoords.size()-1];
  
  int rel_max = max_coord - min_coord;
  ostringstream alignment_text;
  ostringstream assembly_summary;
  alignment_text << "Individual Alignments: (" << num_alignments << ")" << endl;
  for (int i=0; i < alignments.size(); i++) {
    alignment_text << alignments[i].toAlignIllustration(min_coord, rel_max, lineLength) << " index: [" << i << "]" << endl;
  }
  
  if (assemblies.size() != 0) {
    alignment_text << endl << "ASSEMBLIES: (" << assemblies.size() << ")" << endl;
    for (int i=0; i < assemblies.size(); i++) {
      alignment_text << assemblies[i].toAlignIllustration(min_coord, rel_max, lineLength) << " score: (" << assembly_containment_list[i].size() << ") contains [";
      assembly_summary << "assembly: (" << i << ") contains alignments: [";
      vector<int> assemblyIndexList = assembly_containment_list[i];
      for (int j=0; j < assemblyIndexList.size(); j++) {
        int alignmentIndex = assemblyIndexList[j];
        alignment_text << alignmentIndex;
        if (j != assemblyIndexList.size() -1) {
          alignment_text << ",";
        }
        assembly_summary << alignments[alignmentIndex].get_title(); 
        if (j != assemblyIndexList.size() -1) {
          assembly_summary << ",";
        }
      }
      alignment_text << "]" << endl;
      assembly_summary << "]" << " with structure [" << assemblies[i].toString() << "]" << " score: (" << assembly_containment_list[i].size() << ")" << endl;
    }
    alignment_text << assembly_summary.str();
  }
  return (alignment_text.str());
}


Lobject* CDNA_alignment_assembler::get_max_missing_Lobj (vector<Lobject*>& lobj_bin, map<int,bool>& accountedFor) {
  
  int max_missing = 0;
  Lobject* lobj = NULL;
  
  if (DEBUG) { 
    cout << "\n\nget_max_missing_Lobj()" << endl
         << "sifting thru bin of size: " << lobj_bin.size() << endl;
  }
  
  for (int i=0; i < lobj_bin.size(); i++) {
    Lobject* curr_lobj = lobj_bin[i];
    
    if (DEBUG) { cout << "\tanalyzing lobj[" << i << "] " << " at index: " << curr_lobj->index << endl; }
    
    vector<int> alignmentIndices = curr_lobj->getTraceIndices();
    // see if there are any unconsumed alignments in this assembly
    
    int num_missing = 0;
    for (int j=0; j < alignmentIndices.size(); j++) {
      int index = alignmentIndices[j];
      if (! accountedFor[index]) {
        num_missing++;
      }
    }
    if (num_missing > max_missing) {
      max_missing = num_missing;
      lobj = curr_lobj;
    }
  }

  if (DEBUG) { 
    if (max_missing > 0) {
      cout << "-reporting max missing as " << max_missing << " provided by Lobj index: " << lobj->index << endl;
    } else {
      cout << "-none missing in bin this round." << endl;
    }
  }
  
  
  return (lobj); // NULL returned if none found.
}

