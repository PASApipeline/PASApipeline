#ifndef __GCDB_H
#define __GCDB_H
#include <stdlib.h>
#include <stddef.h>
#include <fcntl.h>
#include <sys/stat.h>
#include "gcl/GBase.h"

// If long is natively 64 bit, use the regular fseek and ftell
#ifdef _NATIVE_64
 #define ftello ftell
 #define fseeko fseek
#endif



#ifdef __WIN32__
  #include <io.h>
  #define PROT_READ  1
  #define PROT_WRITE  2
  #define PROT_READWRITE  3
  #define MAP_SHARED  1
  #define MAP_PRIVATE  2
  #define F_OK 0
  #define R_OK 4
  #define W_OK 2
  #define RW_OK 6

  #define ftello ftell
  #define fseeko fseek
  
  #if !defined(MAP_FAILED)
  #define MAP_FAILED      ((void *) -1)
  #endif
  void *mmap(char *,size_t,int,int,int,off_t);
  int   munmap(void *,size_t);
#else
  #include <unistd.h>
  #include <sys/mman.h>
#endif

#define MAX_UINT 0xFFFFFFFFUL


//=====================================================
//-------------     buffer stuff    -------------------
//=====================================================
#define GCDBUFFER_INSIZE 8192
#define GCDBUFFER_OUTSIZE 8192


typedef int (*opfunc)(int, char*, size_t);

//typedef unsigned long gcdb_seek_pos;
typedef off_t gcdb_seek_pos;
typedef unsigned int (*uint_conv)(void*); //uint conversion function pointer
typedef off_t (*offt_conv)(void*); //uint conversion function pointer


//conversion function --> to platform independent uint
extern uint_conv gcvt_uint;
extern offt_conv gcvt_offt;

int endian_test(void);
unsigned int uint32_sun(void* x86int);
unsigned int uint32_x86(void* x86int);
//for file offsets: off_t runtime conversions:
off_t offt_sun(void* offt);
off_t offt_x86(void* offt);


class GCDBuffer {
 public:
  char *x;
  unsigned int p;
  unsigned int n;
  int fd;
  opfunc op;
//methods:
  GCDBuffer() {
    x=NULL;
    fd=0;
    op=NULL;
    n=0;
    //check endianness
    gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;
    gcvt_offt=(endian_test())? &offt_sun : &offt_x86;
    }
  GCDBuffer(opfunc aop,int afd,char *buf,unsigned int len) {
    //check endianness 
    gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;  
    gcvt_offt=(endian_test())? &offt_sun : &offt_x86;
    init(aop, afd, buf, len);
    }
  void init(opfunc aop,int afd,char *buf,unsigned int len) {
     x=buf;
     fd=afd;
     op=aop;
     p=0;
     n=len;
     }
  int  flush();
  int  write_all(char* buf, unsigned int pt);
  int  put(char* buf,unsigned int len);
  int  putalign(char* buf,unsigned int len);
  int  putflush(char* buf,unsigned int len);
  int  puts(char *buf);
  int  putsalign(char *buf);
  int  putsflush(char *buf);
  int  oneRead(char* buf, unsigned int len);
  int  getthis(char* buf,unsigned int len);
  int  get(char* buf,unsigned int len);
  int  bget(char* buf,unsigned int len);
  int  feed();
  char *peek();
  void seek(unsigned int len);
  int copy(GCDBuffer* bin);
};


//=====================================================
//-------------     cdb utils       -------------------
//=====================================================
#ifndef __WIN32__
 extern int errno;
#endif
extern int error_intr;
extern int error_nomem;
extern int error_proto;

//additional data to be appended to the cdb file:
#define CDBMSK_OPT_MULTI    0x00000001
#define CDBMSK_OPT_C        0x00000002
#define CDBMSK_OPT_CADD     0x00000004
#define CDBMSK_OPT_COMPRESS 0x00000008 
//creates a compressed version of the database
//uses plenty of unions for ensuring compatibility with
// the old 'CIDX' info structure

//damn, sun and 64bit machines
// align this to 64bit -- so sizeof() is misled!
#pragma pack(4)
// I wish, but stupid gcc 2.95.3 alpha-decosf version does not
// recognize this pragma directive !!?
// 
struct cdbInfo {
    uint32 num_keys;
    union {
     uint32 num_records;
     char oldtag[4]; // 'CIDX' for old tag style
     };
    // data file size -- used to be  uint32, now it could be 64bit
    union {
     off_t dbsize;
     uint32 oldnum[2]; //num_keys, num_records
     };
    union {  
     uint32 idxflags;
     uint32 old_dbsize;
     };
    union {  
     int dbnamelen;
     int old_idxflags;
     };
      // -- the actual db name precedes this fixed-size record
    union {
     char tag[4]; //'CDBX' for new files with LFS
     uint32 old_dbnamelen;
     };
   };
#pragma pack()   

extern int cdbInfoSIZE;

void uint32_pack(char *,uint32);
void uint32_pack_big(char *,uint32);
void uint32_unpack(char *,uint32 *);
void uint32_unpack_big(char *,uint32 *);

//=====================================================
//-------------     cdb index       -------------------
//=====================================================

#define CDB_HPLIST 1000

struct cdb_hp { uint32 h; uint32 p; } ;

struct cdb_hplist {
  struct cdb_hp hp[CDB_HPLIST];
  struct cdb_hplist *next;
  int num;
  };

//the index file should always be smaller than 4GB !

class GCdbWrite {
   GCDBuffer* cdbuf;
   char bspace[8192];
   char fname[1024];
   char final[2048];
   uint32 count[256];
   uint32 start[256];
   struct cdb_hplist *head;
   struct cdb_hp *split; /* includes space for hash */
   struct cdb_hp *hash;
   uint32 numentries;
   uint32 pos; //file position
   int posplus(uint32 len);
   int fd; //file descriptor
  public:
  //methods:
   GCdbWrite(int afd); //was: init
   GCdbWrite(char* fname);
   ~GCdbWrite();
   int addbegin(unsigned int keylen,unsigned int datalen);
   int addend(unsigned int keylen,unsigned int datalen,uint32 h);
   int addrec(char *key,unsigned int keylen,char *data,unsigned int datalen);
   int add(char *key, char *data, unsigned int datalen);
   int getNumEntries() { return numentries; }
   int finish();
   int close();
   int getfd() { return fd; }
   char* getfile() { return fname; }
};


//=====================================================
//-------------        cdb          -------------------
//=====================================================

#define CDB_HASHSTART 5381

uint32 cdb_hashadd(uint32,unsigned char);
uint32 cdb_hash(char *,unsigned int);

class GCdbRead {
  uint32 size; // initialized if map is nonzero
  uint32 loop; // number of hash slots searched under this key
  uint32 khash; // initialized if loop is nonzero
  uint32 kpos; // initialized if loop is nonzero
  uint32 hpos; // initialized if loop is nonzero
  uint32 hslots; // initialized if loop is nonzero
  uint32 dpos; // initialized if cdb_findnext() returns 1
  uint32 dlen; // initialized if cdb_findnext() returns 1
  char fname[1024];
  char *map; // 0 if no map is available
  int fd;
 public:
//methods:
  GCdbRead(int fd); //was cdb_init
  GCdbRead(char* afname); //was cdb_init
  ~GCdbRead(); //was cdb_free
  int read(char *,unsigned int,uint32);
  int match(char *key, unsigned int len, uint32 pos);
  void findstart() { loop =0; }
  int findnext(char *key,unsigned int len);
  int find(char *key);
  int datapos() { return dpos; }
  int datalen() { return dlen; }
  int getfd() { return fd; }
  char* getfile() { return fname; }
};

#endif
