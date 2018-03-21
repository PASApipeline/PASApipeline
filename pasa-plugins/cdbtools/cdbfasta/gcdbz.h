#ifndef _GCDBZ_H
#define _GCDBZ_H

#define GCDBZ_SBUF_LEN 8192
#define GCDBZ_LBUF_LEN 8192*2
#include <zlib.h>
#include <stdio.h>
#include "gcl/GReadBuf.h"

class GCdbz {
 private:
  char lbuf[GCDBZ_LBUF_LEN]; //larger buffer
  char sbuf[GCDBZ_SBUF_LEN]; //smaller buffer
  char* defline; //defline copy storage -- compression only
  int defline_cap; //currently allocated length of defline
  int defline_len; //currently      used length of defline
  z_stream zstream; // de/compression stream
  FILE* zf; //compressed file, could be input or output
  bool uncompress; // compression or decompression
  long zpos; //current position in zf
  int  zrecsize; // the size of the compressed record
  bool in_defline;
  bool zclosed; // if compress_end() was issued or not!
  void begin_defline() { defline_len=0;
                         in_defline=true; } // initialize the defline storage
  void extend_defline(int ch); //append character ch to defline
                          //reallocating as necessary
  void end_defline() { defline[defline_len]='\0';
                       in_defline=false; } // add \0
 public:
  GCdbz(FILE* af, bool dc = false, int zrsize=0);
  ~GCdbz();
  void compress_start();
  void compress_end();
  char* compress(GReadBuf *readbuf, char* delim);
  // returns a pointer to the defline copy or 
  // NULL if nothing was compressed;
  //   (getZRecSize should be called to find out the 
  //   actual number of compressed bytes written to zf)
  int getZRecSize() { return zrecsize; } 
      //to be called AFTER compress()
  long getZRecPos() {  return zpos; } 
      //to be called BEFORE compress()
  FILE* getZFile() { return zf; }   
  void decomp_start(int zrsize);
  void decomp_end();
  int decompress(FILE* outf, int csize=0, int zfofs=-1);
    // uncompress csize bytes from zf, from optional offset zfofs, 
    // and send the uncompressed stream to outf
};

 



#endif
