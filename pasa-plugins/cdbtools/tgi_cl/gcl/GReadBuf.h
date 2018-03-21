//---------------------------------------------------------------------------
#ifndef GReadBufH
#define GReadBufH
//---------------------------------------------------------------------------
#include "gcl/GBase.h"
#include <fcntl.h>
#include <sys/types.h>
// If long is natively 64 bit, use the regular fseek and ftell
#ifdef _NATIVE_64
 #define ftello ftell
 #define fseeko fseek
#endif

class GReadBuf {
 protected:
  FILE* f;
  uchar* buf;
  int buflen;
  int bufused; //
  int bufpos;
  off_t fpos;
  bool eof;
  bool eob;

  int refill(bool repos=false) {
   //refill the buffer-----------
   if (repos && bufpos==0) return 0; //no need to repos
   if (eof) return 0;
   int fr=0;
   if (repos && bufpos<bufused) {
      int kept=bufused-bufpos;
      memmove((void*)buf, (void*)(buf+bufpos),kept);
      fr=(int)fread((void *)(buf+kept), 1, buflen-kept, f);
      if (fr<buflen-kept) eof=true;
      buf[kept+fr]='\0';
      bufused=kept+fr;
      }
     else {
      fr=(int)fread((void *)buf, 1, buflen, f);
      if (fr<buflen) eof=true;
      buf[fr]='\0'; //only for text record parsers
      bufused=fr;
      }
   if (feof(f)) eof=true;
   if (ferror(f)) {
     GMessage("GReadBuf::refill - error at fread!\n");
     eof=true;
     }
   bufpos=0;
   fpos+=fr; //bytes read from file so far
   return fr;
   }
 public:
  GReadBuf(FILE* fin, int bsize=4096) {
    f=fin;
    buflen=bsize;
    GMALLOC(buf,buflen+1);
    bufpos=0; //current pointer for get function
    bufused=0;
    fpos=0;
    eof=false;
    eob=false;
    refill();
    }
  ~GReadBuf() { GFREE(buf); }
  
  //reads len chars from stream into the outbuf
  //updates bufpos
  //->returns the number of bytes read
  int get(uchar *outbuf, int len) {
    if (eob) return 0;
    int rd=0; //bytes read
    while (!eob && rd<len) {
      int to_read=GMIN((bufused-bufpos),(len-rd));
      memcpy((void*)(outbuf+rd),(void*)(buf+bufpos), to_read);
      bufpos+=to_read;
      rd+=to_read;
      if (bufpos>=bufused) {
        if (eof) eob=true;
           else refill();
        }
      }//while
    return rd;
    }

  uchar* getStr(uchar *outbuf, int len) {
    int rd=get(outbuf,len);
    if (rd==0) return NULL;
      else {
       outbuf[rd]='\0';
       return outbuf;
       }
    }

  // getc equivalent
  int getch() {
    if (eob) return -1;
    int ch=(int)(uchar)buf[bufpos];
    bufpos++;
    if (bufpos>=bufused) {
        if (eof) eob=true;
           else refill();
        }
    return ch;
    }

  //---
  bool isEof() { return eob; }
  bool ended() { return eob; }
  off_t getPos() {
  //returns the virtual file position
  // = the actual file offset of the byte at bufpos
    return fpos-(bufused-bufpos);
    }
  //skip into the stream the specified number of bytes
  int skip(int skiplen) {
   if (eob) return 0;
   int r=0; //the actual number of bytes skipped
   while (skiplen && !eob) {
     int dif=GMIN(bufused-bufpos,skiplen);
     skiplen-=dif;
     bufpos+=dif;
     r+=dif;
     if (bufpos>=bufused) {
       if (eof) { eob=true; return r; }
       refill();
       }
     }
    return r;
   }
  //look ahead without updating the read pointer (bufpos)
  //Cannot peek more than buflen!
  int peek(uchar* outbuf, int len) {
    if (eob) return -1;
    //if (eob || len>buflen) return -1;
    if (len>bufused-bufpos) refill(true);
    int mlen=GMIN((bufused-bufpos),len);
    memcpy((void*)outbuf, (void*)(buf+bufpos), mlen);
    return mlen;
    }

  uchar* peekStr(uchar* outbuf, int len) {
   int rd=peek(outbuf,len);
   if (rd>0) { outbuf[rd]='\0'; return outbuf; }
        else return NULL;
   }
  //looks ahead to check if what follows matches 
  int peekCmp(char* cmpstr, int cmplen=0) {
    if (eob) //GError("GReadBuf::peekcmp error: eob!\n");
         return -2;
    if (!cmplen) cmplen=strlen(cmpstr);
    if (cmplen>bufused-bufpos) {
       refill(true);
       if (cmplen>bufused-bufpos) return -2;
       }
    //use memcmp
    return memcmp((void*)(buf+bufpos), cmpstr, cmplen);
    }
    
};

#endif
