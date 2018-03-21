#ifndef _GCDBYANK_H
#define _GCDBYANK_H

#include "gcl/gcdb.h"
#include "gcl/GFastaFile.h"
// GFastaSeq class
// and *charFunc() callback type

#include <zlib.h>
#include <stdio.h>
#include "gcl/GReadBuf.h"
#include "gcl/GFastaFile.h"

#define GCDBZ_SBUF_LEN 8192
#define GCDBZ_LBUF_LEN 8192*2

#define DEF_CDBREC_DELIM ">"

class GCdbZFasta {
 private:
  char* recdelim;
  char lbuf[GCDBZ_LBUF_LEN]; //larger buffer
  char sbuf[GCDBZ_SBUF_LEN]; //smaller buffer
  char* defline; //defline copy storage -- compression only
  int defline_cap; //currently allocated length of defline
  int defline_len; //currently      used length of defline
  z_stream zstream; // de/compression stream
  FILE* zf; //compressed file
  long zpos; //current position in zf
  int  zrecsize; // the size of the compressed record
  GFastaCharHandler* chrhandler;
 public:
  GCdbZFasta(FILE* af, int zrsize=0, char* r_delim=DEF_CDBREC_DELIM);
  ~GCdbZFasta();
  FILE* getZFile() { return zf; }
  void decomp_start(int zrsize);
  void decomp_end();
  int decompress(FastaSeq& rec, int csize=0, int zfofs=-1, charFunc* seqCallBack=NULL);
    // uncompress csize bytes from file zf, from optional file offset zfofs,
    // and send the uncompressed stream to callbackFn
};

class GCdbYank {
  char* idxfile;
  char* dbfile;
  char* recdelim; //record delimiter -- typically ">"
  int warnings;
  bool is_compressed;
  char* dbname;
  char* info_dbname;
  off_t db_size;
  cdbInfo dbstat;
  GCdbRead* cdb;
  GCdbZFasta* cdbz;
  int fdb;
  int fd;
  FILE* fz; // if compressed
  GFastaCharHandler* fastahandler;
 protected:
   GCdbZFasta* openCdbz(char* p);

 public:
  GCdbYank(const char* fidx, char* recsep=DEF_CDBREC_DELIM);
  ~GCdbYank();
  int getRecord(char* key, FastaSeq& rec, charFunc* seqCallBack=NULL);
  off_t getRecordPos(char* key);
};

#endif
