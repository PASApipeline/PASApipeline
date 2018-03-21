#ifndef LayoutParser_H
#define LayoutParser_H

#include "gcl/GBase.h"
#include "gcl/GList.hh"
#include "gcl/GHash.hh"
#include <stdio.h>
//hash data associated with a contig/sequence name
//a contig name key is always stored as its name plus .<length>

class LytCtgData;

struct LytSeqInterSeg {
 int  segEnd, segRClip;
 int  nextStart, nextLClip;
 char segRSplice, nextLSplice;
 LytSeqInterSeg(int end, int nextstart, int rclip=0, int nextlclip=0, 
                                 char segrsplice=0, char nextlsplice=0) {
  segEnd=end; segRClip=rclip;
  nextStart=nextstart; nextLClip=nextlclip;
  segRSplice=segrsplice;
  nextLSplice=nextlsplice;
  }

};
class LytSeqInfo { //info for a sequence within the file
   int xlen; //total sequence length (with all the added * within contig)
  public:
   char *name;
   LytCtgData* contig; //contig data containing this sequence, as above
   bool segmented;
   int numisegs; // number of intersegs[]
   LytSeqInterSeg* intersegs; 
   unsigned int fpos; //file position for the sequence data
   unsigned char reversed;
   int offs; //offset in contig (of the very left end)
   int left,right; //clear range (relative to sequence itself, max 1..xlen)
   LytSeqInfo(char* seqid, LytCtgData* ctg, int pos=0, unsigned char minus=0,
              int slen=0, int clpL=0, int clpR=0) {
     contig=ctg;
     offs=pos;
     reversed=minus;
     fpos=0;
     xlen=slen;
     left=clpL+1; //1 if no clpL given
     right=xlen-clpR; //0 if no len given
     segmented=false;
     numisegs=0;
     name=Gstrdup(seqid);
     intersegs=NULL;
     }
   ~LytSeqInfo() {
     GFREE(name);
     GFREE(intersegs);
     }
   void addInterSeg(int end, int nextstart, int rclip=0, int nextlclip=0,
                      char splice=0, char nextsplice=0) {    
    GREALLOC(intersegs,(numisegs+1)*sizeof(LytSeqInterSeg));
    intersegs[numisegs].segEnd=end; intersegs[numisegs].segRClip=rclip;
    intersegs[numisegs].nextStart=nextstart; intersegs[numisegs].nextLClip=nextlclip;
    intersegs[numisegs].segRSplice=splice;
    intersegs[numisegs].nextLSplice=nextsplice;
    numisegs++;
   }
   void setLength(int len) {
     //should only be called BEFORE setting the real clipping coordinates
     //(left,right)
     xlen=len;
     left=1;
     right=xlen;
     }
   int length() { return xlen; }
   bool operator==(const LytSeqInfo& s)  {
      return (offs+left-1==s.offs+s.left-1);
     }
   bool operator>(const LytSeqInfo& s)  {
      return (offs+left-1>s.offs+s.left-1);
     }
   bool operator<(const LytSeqInfo& s)  {
      return (offs+left-1<s.offs+s.left-1);
     }
  };

  
class LytCtgData {
   public:
    char* name; //contig name, as stored in file
    unsigned int len; //contig length (lsequence, from ACE file)
    int lpos, rpos;
    int numseqs;
    int offs; //some other type of user data that might be of use
    unsigned int fpos; //position in file for this contig's entry
    GList<LytSeqInfo> seqs;
    LytCtgData(int pos=0):seqs(false,false,false) {
      name=NULL;
      offs=0;
      len=0;
      numseqs=0;
      fpos=pos;
      }
    ~LytCtgData() {
      GFREE(name);
      seqs.Clear();
      }

   char* readName(char* s, GHash<int>& names);
   
   bool operator==(const LytCtgData& s)  {
      return (strcmp(name,s.name)==0);
     }
   bool operator>(const LytCtgData& s)  {
      return (strcmp(name,s.name)>0);
     }
   bool operator<(const LytCtgData& s)  {
      return (strcmp(name,s.name)<0);
     }
 };

typedef bool fnLytSeq(int ctgno, LytCtgData* d, LytSeqInfo* s, char* seq);

class LayoutParser {
 protected:
  FILE* f; //file stream
  long f_pos;
  char* fname;
  LytCtgData* currentContig; // currently loaded contig -- for browsing/loading 
  int numContigs; //total number of contigs found in this file
  //int numSeqs; //total number of (distinct) sequences found in this file
  GHash<LytSeqInfo> seqinfo; //sequence locations in the file
  GHash<int> ctgIDs; //contig IDs, to make them unique!

  GList<LytCtgData> contigs; //list of contig names with their size,
                       //number of sequences and filepos
 protected:
  GLineBuf* linebuf; //the line buffer
  int fskipTo(const char* linestart, char* butnot=NULL);
  bool startsWith(const char* s, const char* start, int tlen);
  virtual LytSeqInfo* addSeq(char* s, LytCtgData* ctg);
  int seek(int offset) {
      int r=fseek(f, offset, SEEK_SET);
      if (r==0) f_pos=offset;
      return r;
      }

 public:
  LayoutParser(const char* filename):contigs(false,true) {
   f=NULL;
   f_pos=0;
   numContigs=0;
   currentContig=NULL;
   if (filename==NULL) {
      f=stdin;
      fname=Gstrdup("stdin");
      }
    else
      fname=Gstrdup(filename);
   linebuf=new GLineBuf();
   }
  virtual ~LayoutParser() {
    ctgIDs.Clear();
    GFREE(fname);
    delete linebuf;
    close();
    numContigs=0;
    seqinfo.Clear();
    contigs.Clear();
    }
  virtual bool open();
  void close();
  virtual bool parse(fnLytSeq* seqfn=NULL); //load all the file offsets
  virtual bool parseContigs(); //load all the file offsets for contigs
  virtual bool loadContig(int ctgidx, fnLytSeq* seqfn=NULL, bool re_pos=true); //for loading by browsing
  //if parsefn is not NULL, it is executed, passing the sequence data(first time, with the contig sequence)
  //if parserfn returns true, the data is freed after it is processed
  virtual char getFileType() { return 'L'; }
  //sequence loading - only by request
  LytCtgData* getContig(int idx) { return contigs[idx]; }
  virtual char* getSeq(LytSeqInfo* sqinfo) { return NULL; }
  virtual char* getContigSeq(LytCtgData* ctgdata) { return NULL; }
  int getNumContigs() { return numContigs; }
  //int getNumSeqs() { return numSeqs; }
  int getFilePos() { return f_pos; }
  //sorting the list of contigs:
  void contigsByName();
  void contigsByLen();
  void contigsByNumSeqs();
};

#endif
