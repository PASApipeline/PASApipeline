#include "gcl/LayoutParser.h"

bool LayoutParser::startsWith(const char* s, const char* start, int tlen) {
 bool found=true;
 for (int i=0;i<tlen;i++)
  if (s[i]!=start[i])
     return false;
 return found;
 }



bool LayoutParser::open() { //also checks if it looks like a valid LYT file
   if (f==stdin) return true;
   if ((f=fopen(fname,"rb"))==NULL) return false;
     char first1[2];
     first1[1]='\0';
     if (fread(first1, 1,1, f)==1) {
       f_pos=1;
       return (strcmp(first1, "#")==0 || strcmp(first1, ">")==0);
       }
   return false;
   }


void LayoutParser::close() {
   if (f!=NULL) {
        fclose(f);
        f=NULL;
        }
   }

int LayoutParser::fskipTo(const char* linestart, char* butnot) {
  /* reads the file from the current position
  until the next occurence of a line
  starting with linestart
  returns the line in buf and the file offset
  of the beginning of the line;
  the file offset is -1 if linestart was not found
  or -2 if an unwanted line-start came out of order
  */
    int lastpos=getFilePos();
    int tlen=strlen(linestart);  
    int nlen=(butnot==NULL) ? 0 : strlen(butnot);
    while (linebuf->getLine(f, f_pos)!=NULL) {
     if (nlen>0 && startsWith(linebuf->line(), butnot, nlen)) {
        GMessage("fSkipTo: unwanted line '%s' encountered when searching for '%s'\n",
          linebuf->line(), linestart);
        return -2;
        }
     if (startsWith(linebuf->chars(), linestart, tlen)) return lastpos;
     lastpos=getFilePos();
     }
 return -1;   
 }

LytSeqInfo* LayoutParser::addSeq(char* s, LytCtgData* ctg) {
 LytSeqInfo* seq;
 //s must be the line with sequence data
 char* p=strchrs(s," \t");
 if (p==NULL) return NULL;
 p++;
 char c;
 int slen, soffs, clpL, clpR;
 clpL=0;clpR=0;
 if (sscanf(p,"%c %d %d %d %d", &c, &slen, &soffs, &clpL, &clpR)<3) return NULL;
 p--;
 *p='\0';
 if ((seq=seqinfo.Find(s))!=NULL) {
   GMessage("Sequence '%s' already found for contig '%s (%d nt)'\n"
      " so it cannot be added for contig '%s (%d nt)'\n",
     s, seq->contig->name, seq->contig->len,
     ctg->name, ctg->len);
   return NULL;
   }
 seq = new LytSeqInfo(s, ctg, soffs, (c=='-') ? 1 : 0, slen, clpL, clpR);
 seqinfo.shkAdd(seq->name, seq);
 ctg->seqs.Add(seq);
 //parse optional extensions, if any
 p+=strlen(s); //position p after the seqname
 char* m=NULL;
 int segEnd, segRclip,nextsegStart, nextsegLclip, prevSegStart;
 char segSplice, nextsegSplice;
 while ((m=strchr(p,':'))!=NULL) {
  switch (*(m-1)) {
    case 'G': //segmenting info
       prevSegStart=soffs+clpL-1;
       p=m+1;  //p to the beginning of G: data
       //accumulate the total length in lenSegs
       while (*p>='1' && *p<='9') {
         segEnd=0;
         segRclip=0;
         nextsegStart=0;
         nextsegLclip=0;
         segSplice=0;
         nextsegSplice=0;
         if (!parseInt(p,segEnd))
            GError("Error [segment] at LayoutParser for %s at: %s\n",
                 s, m-1);
         if (*p=='c') {
            p++;
            if (!parseInt(p,segRclip))
              GError("Error [segment] at LayoutParser for %s at: %s\n",
                      s, m-1);
            }
         if (*p=='S' || *p=='s') {
            segSplice=*p; p++;
            }
         if (*p!='-')
                GError("Error [segment] at LayoutParser for %s at: %s\n",
                      s, m-1);
            else p++;
         if (!parseInt(p,nextsegStart))
            GError("Error [segment] at LayoutParser for %s at: %s\n",
                 s, m-1);
         if (*p=='c') {
            p++;
            if (!parseInt(p,nextsegLclip))
              GError("Error [segment] at LayoutParser for %s at: %s\n",
                      s, m-1);
            }
         if (*p=='S' || *p=='s') {
            nextsegSplice=*p; p++;
            }
         seq->addInterSeg(segEnd,nextsegStart,segRclip,nextsegLclip, 
                                         segSplice, nextsegSplice);
         prevSegStart=nextsegStart;
         //
         if (*p==',') p++;
             else break;
         } //while inter-segment parsing
       break; // 'G:' case
    case 'L': //clone mates list
       p=m+1; //p to the beginning of L: data
       break;
    case 'D': //difference sequence
       p=m+1; //p to the beginning of D: data
       break;
    case 'S': //actual sequence
       p=m+1; //p to the beginning of S: data
       break;
    default:
       p=m+1;//next attribute
    }
  }

 return seq;
}

char* LytCtgData::readName(char* s, GHash<int>& names) {
   char* p=strchrs(s, " \t");
   if (p!=NULL) {
       char* tmp;
       char* tmp2;
       GMALLOC(tmp, (p-s+30)*sizeof(char));
       strncpy(tmp, s,p-s);
       tmp[p-s]='\0';
       GMALLOC(tmp2, (p-s+30)*sizeof(char));
       strcpy(tmp2, tmp);
       //make it unique (by simple versioning)
       int v=0;
       while (names.hasKey(tmp2)) {
         v++;
         sprintf(tmp2, "%s.%d", tmp, v);
         }
       name=Gstrdup(tmp2);
       GFREE(tmp);
       GFREE(tmp2);
       names.shkAdd(name, new int(1));  
       p++;
       }
      else {
       GMessage("LytCtgData::readName: Cannot find the token delimiter in:\n%s\n", s);
       }
     return p;
     }



bool LayoutParser::parse(fnLytSeq* seqfn) {
  //read all seqs and their positions from the file
  //also checks for duplicate seqnames (just in case)
  if (f!=stdin) seek(0);
  ctgIDs.Clear();
  //GHash<int> ctgIDs; //contig IDs, to make them unique!
  //
  int ctgpos; 
  numContigs=0;
  while ((ctgpos=fskipTo(">"))>=0) { //locate the contig line
    numContigs++;
    LytCtgData* ctgdata=new LytCtgData(ctgpos);
    char* p=ctgdata->readName(linebuf->chars()+1, ctgIDs);
    if (p==NULL) {
       GMessage("LayoutParser: error parsing contig name:\n%s\n", linebuf->chars());
       return false;
       }
    int ctg_lpos, ctg_rpos, numseqs;
    //p must be after contig name within linebuf!
    ctg_lpos=0;ctg_rpos=0;
    if (sscanf(p, "%d %d %d", &numseqs, &ctg_lpos, &ctg_rpos)<1) {
      GMessage("Error parsing contig len and seq count at:\n%s\n",
         p);
      return false;
      }
    //ctg_numSeqs=numseqs;
    ctgdata->numseqs=numseqs;
    ctgdata->rpos=ctg_rpos;
    ctgdata->lpos=ctg_lpos;
    ctgdata->len=ctg_rpos-ctg_lpos+1;
    ctgdata->offs=ctg_lpos;    
    int ctgidx=contigs.Add(ctgdata);
    //now look and load all the component sequences
    loadContig(ctgidx, seqfn, false);
    } //while lines
  contigs.setSorted(true);
  return true;
  }

//==============================================
/*
 Load contig data; can be called by parse - and then no fseek is needed and 
 the file position if right after parsing the contig summary data
*/
bool LayoutParser::loadContig(int ctgidx, fnLytSeq* seqfn, bool re_pos) {
    bool forgetCtg=false;
    char* r=NULL;
    if (ctgidx>=contigs.Count()) 
      GError("LayoutParser: invalid contig index '%d'\n", ctgidx);
      
    LytCtgData* ctgdata=contigs[ctgidx];
    if (re_pos && currentContig!=NULL) { //free previous contig data                                         
                                          //unless it was a parse() call
      currentContig->seqs.Clear();
      seqinfo.Clear();
      }
    currentContig=ctgdata;
    if (re_pos) {
       seek(ctgdata->fpos); //position right where the contig definition starts
       r=linebuf->getLine(f,f_pos);//skip the first line
       if (r==NULL) return false;
       }
    if (seqfn!=NULL)
       forgetCtg=(*seqfn)(numContigs, ctgdata, NULL, NULL);
    int ctg_numSeqs=ctgdata->numseqs;
    int numseqs=0;
    while ((r=linebuf->getLine(f,f_pos))!=NULL) {
       if (linebuf->length()<4) continue;
       if (linebuf->chars()[0]=='>') break; //reached next contig
       //sequence data parsing
       
       bool forgetSeq=false;
       LytSeqInfo* seq=NULL;
       if ((seq=addSeq(linebuf->chars(), ctgdata))==NULL) {
         GMessage("LayoutParser: error parsing sequence entry:\n%s\n",linebuf->chars());
         return false;
         }
        /*
        // Weird -- why would I MODIFY the given clipping of a sequence?
        //--
        bool ctg_clipping = (ctgdata->rpos>ctgdata->lpos);
        if (ctg_clipping) {
          if (ctgdata->lpos > seq->offs && ctgdata->lpos < seq->offs+seq->length())
             seq->left = ctgdata->lpos - seq->offs+1;
            if (ctgdata->rpos < seq->offs+seq->length() && ctgdata->rpos>seq->offs )
             seq->right = ctgdata->rpos-seq->offs+1;
          } */
        if (seqfn!=NULL)
          forgetSeq=(*seqfn)(numContigs, ctgdata, seq, NULL);
        if (forgetSeq) {
            ctg_numSeqs--;
            seqinfo.Remove(seq->name);
            ctgdata->seqs.RemovePtr(seq);
            }
          else {
            numseqs++;
            }
       } //while sequences
     if (forgetCtg) {
      ctgIDs.Remove(ctgdata->name);
      contigs.RemovePtr(ctgdata);
      }
    if (numseqs!=ctg_numSeqs) {
       GMessage("Mismatching number of sequences found (%d) for contig '%s' "
         "(length %d, numseqs %d)\n", numseqs,
                ctgdata->name, ctgdata->len, ctg_numSeqs);
       return false;
       }
return true;      
}


bool LayoutParser::parseContigs() { //load all the file offsets for contigs
  if (f!=stdin) seek(0);
  ctgIDs.Clear();
  //GHash<int> ctgIDs; //contig IDs, to make them unique!
  //
  int ctgpos; //locate the first contig line
  numContigs=0;
  while ((ctgpos=fskipTo(">"))>=0) {
    numContigs++;
    LytCtgData* ctgdata=new LytCtgData(ctgpos);
    char* p=ctgdata->readName(linebuf->chars()+1, ctgIDs);
    if (p==NULL) {
       GMessage("LayoutParser: error parsing contig name:\n%s\n", linebuf->chars());
       return false;
       }
    int ctg_lpos, ctg_rpos, numseqs;
    //p must be after contig name within linebuf!
    ctg_lpos=0;ctg_rpos=0;
    if (sscanf(p, "%d %d %d", &numseqs, &ctg_lpos, &ctg_rpos)<1) {
      GMessage("Error parsing contig len and seq count at:\n%s\n",
         p);
      return false;
      }
    ctgdata->len=ctg_rpos-ctg_lpos+1;
    ctgdata->numseqs=numseqs;
    ctgdata->lpos=ctg_lpos;
    ctgdata->rpos=ctg_rpos;
    ctgdata->offs=ctg_lpos;    
    contigs.Add(ctgdata);
    //ctgpos=fskipTo(">");
    } //while lines
  if (ctgpos==-2) return false;
  contigs.setSorted(true);
  return true;
}

//-- compare functions for contigs
int ctgByLen(void* p1, void* p2) {
 int c1=((LytCtgData*)p1)->len;
 int c2=((LytCtgData*)p2)->len;
 return (c1>c2)?-1:((c1<c2)?1:0);
}

int ctgByNumSeqs(void* p1, void* p2) {
 int c1=((LytCtgData*)p1)->numseqs;
 int c2=((LytCtgData*)p2)->numseqs;
 return (c1>c2)?-1:((c1<c2)?1:0);
}

void LayoutParser::contigsByName() {
 contigs.setSorted(true);
}

void LayoutParser::contigsByLen() {
 contigs.setSorted(&ctgByLen);
}

void LayoutParser::contigsByNumSeqs() {
 contigs.setSorted(&ctgByNumSeqs);
}
