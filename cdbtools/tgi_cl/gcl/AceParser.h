#ifndef ACEPARSER_H
#define ACEPARSER_H
#include "gcl/LayoutParser.h"

class AceParser : public LayoutParser {
 protected: 
  virtual LytSeqInfo* addSeq(char* s, LytCtgData* ctg);
  char* readSeq(); //assumes the next line is just sequence data
                   //reads everything until after the next empty line
  
 public:
  AceParser(const char* filename):LayoutParser(filename) {}
  virtual bool open();
  virtual bool parse(fnLytSeq* seqfn=NULL); //load all the file offsets
  virtual bool parseContigs(); //load contigs' file offsets
  virtual bool loadContig(int ctgidx, fnLytSeq* seqfn=NULL, 
                     bool re_pos=true); //for loading by browsing
  //sequence loading - only by request
  virtual char getFileType() { return 'A'; }
  virtual char* getSeq(LytSeqInfo* seqinfo);
  virtual char* getContigSeq(LytCtgData* ctgdata);  
};

#endif
