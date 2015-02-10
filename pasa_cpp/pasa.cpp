// Program to Assemble Spliced Alignments
// software written by Brian Haas (bhaas@tigr.org)
// PASA algorithm devised by Arthur Delcher (adelcher@tigr.org)

#include "cdna_alignment.h"
#include "alignment_segment.h"
#include "cdna_alignment_assembler.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "stringFuncts.h"
#include "argProcessor.h"

using namespace std;


bool DEBUG = false;
const int MAX_INPUTLINE_LENGTH = 10000;


int main (int argc, char* argv[]) {
  if (argc < 2) {
    cerr << endl << endl << "\tUsage: " << argv[0] <<  " inputFile [opts]" << endl;
    cerr << endl << "\toptions: " << endl << endl
         << "\t-F fuzzlength (bp to discount at alignment termini during pairwise" << endl
         << "\t   compatibility checks.  (default: 20)" << endl
         << "\t-a illustrate incoming alignments only.  No assembly performed." << endl
         << "\t-v verbose" << endl << endl << endl ;
    exit(1);
  }
  
  processArgs(argc, argv);
  
  char* inputFile = argv[1];
  ifstream fileReader (inputFile);
  if (fileReader == 0) { // couldn't open file
    cerr << "Could not open " << inputFile << endl;
    exit(1);
  }
  
  if (argSet["-v"]) {
    DEBUG = true;
  }
  
  string fuzzlength = argVal["-F"];
  int i_fuzzlength = -1;
  if (fuzzlength != "") {
    int fuzzdist = stringToInt(fuzzlength);
    if (fuzzdist >= 0 && fuzzdist <= 100) {
      i_fuzzlength = fuzzdist;
    } else {
      cerr << "ERROR: specified fuzzdist of " << fuzzlength << " is not acceptable. " << endl;
      exit(2);
    }
  }
  
  vector<CDNA_alignment> cdnaList;
  cout << "//" << endl; // output record separator.
  char c_line[MAX_INPUTLINE_LENGTH];
  fileReader.getline(c_line, MAX_INPUTLINE_LENGTH);
  while (! fileReader.eof()) {
    if (fileReader.rdstate() & std::ios::failbit) {
      cerr << "Fatal Error reading input file." << endl;
      exit(5);
    }
    
    string line (c_line);
    if (line.find(",") != -1) {
      cout << "input: " << line << endl;
      vector<string> tokens = stringSplitter(line, ",");
      string acc = tokens[0];
      string orient = tokens[1];
      char orientation = orient[0];
      vector<Alignment_segment> seglist;
      for (int i=2; i< tokens.size(); i++) {
        string coordPair = tokens[i];
        vector<string> coords = stringSplitter(coordPair, "-");
        int lend = stringToInt(coords[0]);
        int rend = stringToInt(coords[1]);
        seglist.push_back(Alignment_segment(lend,rend));
      }
      CDNA_alignment align (seglist, orientation);
      align.set_title(acc);
      cdnaList.push_back(align);
      
      fileReader.getline(c_line, MAX_INPUTLINE_LENGTH);
    }
  }
  fileReader.close();
  
  CDNA_alignment_assembler assembler (cdnaList);
  if (! argSet["-a"]) {
    if (i_fuzzlength >= 0) {
      assembler.set_fuzzlength(i_fuzzlength);
    }
    assembler.assembleAlignments();
  }
  
  // print assemblies.
  cout << endl << assembler.toAlignIllustration(70) << endl << endl;
  
  
  return(0);
}



