#include "gcl/gcdb.h"
#include <errno.h>
#include <stdio.h>
#include <stdarg.h>
#ifdef __WIN32__
  #include <windows.h>
/*   m m a p           ===      got from imagick sources
%  Method mmap emulates the Unix method of the same name.
%  The format of the mmap method is:
%    void *mmap(char *address,size_t length,int protection,
%      int access,int file,off_t offset)
*/
void *mmap(char *address,size_t length,int protection,int access,
  int file, off_t offset) {
  void *map;
  HANDLE handle;
  map=(void *) NULL;
  handle=INVALID_HANDLE_VALUE;
  switch (protection)
  {
    case PROT_READ:
    default:
    {
      handle=CreateFileMapping((HANDLE) _get_osfhandle(file),0,PAGE_READONLY,0,
        length,0);
      if (!handle)
        break;
      map=(void *) MapViewOfFile(handle,FILE_MAP_READ,0,0,length);
      CloseHandle(handle);
      break;
    }
    case PROT_WRITE:
    {
      handle=CreateFileMapping((HANDLE) _get_osfhandle(file),0,PAGE_READWRITE,0,
        length,0);
      if (!handle)
        break;
      map=(void *) MapViewOfFile(handle,FILE_MAP_WRITE,0,0,length);
      CloseHandle(handle);
      break;
    }
    case PROT_READWRITE:
    {
      handle=CreateFileMapping((HANDLE) _get_osfhandle(file),0,PAGE_READWRITE,0,
        length,0);
      if (!handle)
        break;
      map=(void *) MapViewOfFile(handle,FILE_MAP_ALL_ACCESS,0,0,length);
      CloseHandle(handle);
      break;
    }
  }
  if (map == (void *) NULL)
    return((void *) MAP_FAILED);
  return((void *) ((char *) map+offset));
}

/*  =========== m u n m a p ===========================
%
%  Method munmap emulates the Unix method with the same name.
%  The format of the munmap method is:
%      int munmap(void *map,size_t length)
%  A description of each parameter follows:
%    > status:  Method munmap returns 0 on success; otherwise, it
%      returns -1 and sets errno to indicate the error.
%    > map: The address of the binary large object.
%    > length: The length of the binary large object.
%
*/
int munmap(void *map,size_t length) {
  if (!UnmapViewOfFile(map))
    return(-1);
  return(0);
}

#endif



int cdbInfoSIZE=offsetof(cdbInfo, tag)+4;


//=====================================================
//-------------     buffer stuff    -------------------
//=====================================================

//-------------------------------------
//--------- misc utility functions -----

static int gcdb_seek_set(int fd,gcdb_seek_pos pos) {
 if (lseek(fd, pos, 0) == -1)
         return -1;
 return 0;
 }

#define gcdb_seek_begin(fd) (gcdb_seek_set((fd),(gcdb_seek_pos) 0))

static unsigned int gcdb_strlen(char *s) {
  register char *t;
  t = s;
  for (;;) {
    if (!*t) return t - s; ++t;
    if (!*t) return t - s; ++t;
    if (!*t) return t - s; ++t;
    if (!*t) return t - s; ++t;
  }
}


static int byte_diff(char *s, unsigned int n,char *t) {
  for (;;) {
    if (!n) return 0; if (*s != *t) break; ++s; ++t; --n;
    if (!n) return 0; if (*s != *t) break; ++s; ++t; --n;
    if (!n) return 0; if (*s != *t) break; ++s; ++t; --n;
    if (!n) return 0; if (*s != *t) break; ++s; ++t; --n;
  }
  return ((int)(unsigned int)(unsigned char) *s)
       - ((int)(unsigned int)(unsigned char) *t);
}

static void gcdb_byte_copy(char *to, unsigned int n, char *from) {
  for (;;) {
    if (!n) return; *to++ = *from++; --n;
    if (!n) return; *to++ = *from++; --n;
    if (!n) return; *to++ = *from++; --n;
    if (!n) return; *to++ = *from++; --n;
  }
}

static void gcdb_byte_copyr(char *to, unsigned int n, char *from) {
  to += n;
  from += n;
  for (;;) {
    if (!n) return; *--to = *--from; --n;
    if (!n) return; *--to = *--from; --n;
    if (!n) return; *--to = *--from; --n;
    if (!n) return; *--to = *--from; --n;
  }
}

#define ALIGNMENT 16 /* XXX: assuming that this alignment is enough */
#define SPACE 4096 /* must be multiple of ALIGNMENT */

typedef union { char irrelevant[ALIGNMENT]; double d; } aligned;
static aligned realspace[SPACE / ALIGNMENT];
#define space ((char *) realspace)
      
static unsigned int avail = SPACE; /* multiple of ALIGNMENT; 0<=avail<=SPACE */

offt_conv gcvt_offt;
uint_conv gcvt_uint;


char *gcdb_alloc(unsigned int n) {
  char *x;
  n = ALIGNMENT + n - (n & (ALIGNMENT - 1)); /* XXX: could overflow */
  if (n <= avail) { avail -= n; return space + avail; }
  x = (char*) malloc(n);
  if (!x) return NULL;
  //if (!x) GError("Error: mgcdb_alloc(%d) failed !\n", n);
  return x;
}


int GCDBuffer::write_all(char* buf, unsigned int len) {
  int w;
  while (len) {
    w = op(fd,buf,len);
    if (w == -1) {
      if (errno == error_intr) continue;
      return -1; /* note that some data may have been written */
    }
    /*  if (w == 0) ;  luser's fault */
    buf += w;
    len -= w;
  }
  return 0;
}

int GCDBuffer::flush() {
  int pt=p;
  if (!pt) return 0;
  p = 0;
  //return allwrite(op,fd,x,pt);
  return write_all(x,pt);
}

int GCDBuffer::putalign(char *buf,unsigned int len) {
  unsigned int bn;

  while (len > (bn = n-p)) {
    gcdb_byte_copy(x + p,bn,buf); 
    p += bn; buf += bn; len -= bn;
    if (GCDBuffer::flush() == -1) return -1;
    }

  /* now len <= s->n - s->p */
  gcdb_byte_copy(x + p,len,buf);
  p += len;
  return 0;
}

int GCDBuffer::put(char *buf,unsigned int len) {
  unsigned int bn=n;
  if (len > bn - p) {
    if (GCDBuffer::flush() == -1) return -1;
    /* now s->p == 0 */
    if (bn < GCDBUFFER_OUTSIZE) bn = GCDBUFFER_OUTSIZE;
    while (len > n) {
      if (bn > len) bn = len;
      if (write_all(buf, bn) == -1) return -1;
      buf += bn;
      len -= bn;
    }
  }  
  /* now len <= s->n - s->p */
  gcdb_byte_copy(x + p,len,buf);
  p += len;
  return 0;
}

int GCDBuffer::putflush(char *buf,unsigned int len) {
  if (flush() == -1) return -1;
  return write_all(buf,len);
}

int GCDBuffer::putsalign(char *buf) {
  return GCDBuffer::putalign(buf, gcdb_strlen(buf));
}

int GCDBuffer::puts(char *buf) {
  return GCDBuffer::put(buf, gcdb_strlen(buf));
}

int GCDBuffer::putsflush(char *buf) {
  return GCDBuffer::putflush(buf, gcdb_strlen(buf));
}

static int oneread(opfunc op,int fd, char *buf,unsigned int len) {
  int r;
  for (;;) {
    r = op(fd,buf,len);
    if (r == -1 && errno == error_intr) continue;
    return r;
    }
}

int GCDBuffer::oneRead(char* buf, unsigned int len) {
  return op(fd,buf,len);
  /*int r;
  for (;;) {
    r = op(fd,buf,len);
    if (r == -1 && errno == error_intr) continue;
    return r;
    }*/
}

int GCDBuffer::getthis(char *buf,unsigned int len) {
  if (len > p) len = p;
  p -= len;
  gcdb_byte_copy(buf, len,x + n);
  n += len;
  return len;
}

int GCDBuffer::feed() {
  int r;
  if (p) return p;
  r = oneRead(x,n);
  if (r <= 0)
     return r;
  p = r;
  n -= r;
  if (n > 0) gcdb_byte_copyr(x + n,r,x);
  return r;
}

int GCDBuffer::bget(char *buf,unsigned int len) {
  int r;
  if (p > 0) return getthis(buf,len);
  if (n <= len) return oneRead(buf,n);
  r = GCDBuffer::feed(); if (r <= 0) return r;
  return getthis(buf,len);
}

int GCDBuffer::get(char *buf,unsigned int len) {
  int r; 
  if (p > 0) return getthis(buf,len);
  if (n <= len) return oneread(op,fd,buf,len);
  r = GCDBuffer::feed();
  if (r <= 0)
      return r;
  return getthis(buf,len);
}

char* GCDBuffer::peek() {
  return x + n;
}

void GCDBuffer::seek(unsigned int len) {
  n += len;
  p -= len;
}

int GCDBuffer::copy(GCDBuffer* bin) {
  int n_in;
  char *x_in;
  for (;;) {
    n_in = bin->feed();
    if (n_in < 0) return -2;
    if (!n_in) return 0;
    x_in = bin->peek();
    if (GCDBuffer::put(x_in,n_in) == -1) return -3;
    bin->seek(n_in);
   }
}

//=====================================================
//-------------     cdb utils       -------------------
//=====================================================

int error_intr =
#ifdef EINTR
EINTR;
#else
-1;
#endif

int error_nomem =
#ifdef ENOMEM
ENOMEM;
#else
-2;
#endif

int error_proto =
#ifdef EPROTO
EPROTO;
#else
-15;
#endif
//------------------------------------------------
//------------    allocation routines:

/*
 big/little endian check
*/
int endian_test(void) {
 unsigned short v=0x0001;
 unsigned char* b = (unsigned char*)&v;
 return b[1];
}
/* conversion of unsigned int offsets read from a file 
   can also be used to prepare unsigned integers to be written 
   into a file in an independent platform manner
*/

unsigned int uint32_sun(void* x86int) {
 unsigned char b[4];
 b[3]=((unsigned char*)x86int)[0];
 b[0]=((unsigned char*)x86int)[3];
 b[1]=((unsigned char*)x86int)[2];
 b[2]=((unsigned char*)x86int)[1];
 return *((unsigned int*)b);
}


unsigned int uint32_x86(void* offt) {
 return *((unsigned int*)offt);
}

//-------- 64bit types, if that's the case:

off_t offt_sun(void* offt) {
 unsigned char b[8];
 if (sizeof(off_t)==8) { //64 bit?
  // upper words:
  b[3]=((unsigned char*)offt)[4];
  b[0]=((unsigned char*)offt)[7];
  b[1]=((unsigned char*)offt)[6];
  b[2]=((unsigned char*)offt)[5];
  //--
  b[7]=((unsigned char*)offt)[0];
  b[4]=((unsigned char*)offt)[3];
  b[5]=((unsigned char*)offt)[2];
  b[6]=((unsigned char*)offt)[1];
  }
 else {
  b[3]=((unsigned char*)offt)[0];
  b[0]=((unsigned char*)offt)[3];
  b[1]=((unsigned char*)offt)[2];
  b[2]=((unsigned char*)offt)[1]; 
  }
 return *((off_t*)b);
}



off_t offt_x86(void* offt) {
 return *((off_t*)offt);
}



//------------------------ platform independent uint32 :

void uint32_pack(char s[4],uint32 u)
{
  s[0] = u & 255;
  u >>= 8;
  s[1] = u & 255;
  u >>= 8;
  s[2] = u & 255;
  s[3] = u >> 8;
}

void uint32_pack_big(char s[4],uint32 u)
{
  s[3] = u & 255;
  u >>= 8;
  s[2] = u & 255;
  u >>= 8;
  s[1] = u & 255;
  s[0] = u >> 8;
}

/* unpacking: */


void uint32_unpack(char s[4],uint32 *u)
{
  uint32 result;

  result = (unsigned char) s[3];
  result <<= 8;
  result += (unsigned char) s[2];
  result <<= 8;
  result += (unsigned char) s[1];
  result <<= 8;
  result += (unsigned char) s[0];

  *u = result;
}

void uint32_unpack_big(char s[4],uint32 *u)
{
  uint32 result;

  result = (unsigned char) s[0];
  result <<= 8;
  result += (unsigned char) s[1];
  result <<= 8;
  result += (unsigned char) s[2];
  result <<= 8;
  result += (unsigned char) s[3];

  *u = result;
}

//=====================================================
//-------------     cdb index       -------------------
//=====================================================

GCdbWrite::GCdbWrite(int afd) {
  //check endianness :)
  gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;
  gcvt_offt=(endian_test())? &offt_sun : &offt_x86;
  cdbuf=new GCDBuffer((opfunc)&write,(int) afd,(char*)bspace,sizeof bspace);
  head = NULL;
  split = 0;
  hash = 0;
  numentries = 0;
  fd = afd;
  pos = sizeof final;
  gcdb_seek_set(fd, pos);

  fname[0]='\0';
  //should return and test the result of gcdb_seek_set!!!
}

GCdbWrite::GCdbWrite(char* afname) {
#ifdef __WIN32__
   fd = open(afname,O_WRONLY | O_TRUNC | O_BINARY | O_CREAT, S_IREAD|S_IWRITE);
#else
   fd = open(afname,O_WRONLY | O_NDELAY | O_TRUNC | O_CREAT, 0664);
#endif
  if (fd == -1)
    GError("GCdbWrite: Error creating file '%s'\n", fname);
    
 //check endianness :)
  gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;
  gcvt_offt=(endian_test())? &offt_sun : &offt_x86;
    
  cdbuf=new GCDBuffer((opfunc)&write,(int) fd,(char*)bspace,sizeof bspace);
  head = NULL;
  split = 0;
  hash = 0;
  numentries = 0;
  pos = sizeof final;
  gcdb_seek_set(fd, pos);
  strcpy(fname, afname);

  //should return and test the result of gcdb_seek_set!!!
}

GCdbWrite::~GCdbWrite() {
  cdbuf->flush();
  #ifndef __WIN32__
   /* NFS silliness  */
   if (fsync(fd) == -1)
      GError("GCdbWrite: Error at fsync() for file '%s'\n",
                          fname);
  #endif
  if (::close(fd) == -1)
        GError("GCdbWrite: Error at closing file '%s'\n",
                          fname);
  delete cdbuf;
  if (head!=NULL) free(head);
  }

int GCdbWrite::posplus(uint32 len) {
  uint32 newpos = pos + len;
  if (newpos < len) { //errno = error_nomem; 
                     return -1; }
  pos = newpos;
  return 0;
}

int GCdbWrite::addend(unsigned int keylen,unsigned int datalen,uint32 h) {
  struct cdb_hplist *chead = head;
  if (!chead || (chead->num >= CDB_HPLIST)) {
    chead = (struct cdb_hplist *) gcdb_alloc(sizeof(struct cdb_hplist));
    if (!chead) return -1;
    chead->num = 0;
    chead->next = head;
    head = chead;
    }
  chead->hp[head->num].h = h;
  chead->hp[head->num].p = pos;
  ++chead->num;
  ++numentries;
  if (posplus(8) == -1) return -1;
  if (posplus(keylen) == -1) return -1;
  if (posplus(datalen) == -1) return -1;
  return 0;
}

int GCdbWrite::addbegin(unsigned int keylen,unsigned int datalen) {
  char buf[8];
  //if (keylen > MAX_UINT) { /* errno = error_nomem; */return -1; }
  // if (datalen > MAX_UINT) { /*errno = error_nomem;*/ return -1; } 
  uint32_pack(buf,keylen);
  uint32_pack(buf + 4,datalen);
  if (cdbuf->putalign(buf,8) == -1) return -1;
  return 0;
}

#define cdbuffer_PUTC(s,c) \
  ( ((s).n != (s).p) \
    ? ( (s).x[(s).p++] = (c), 0 ) \
    : (s).put(&(c),1) \
  )

int GCdbWrite::add(char* key, char* recdata, unsigned int datalen) {
 unsigned int i;
 unsigned int klen=strlen(key);
 if (klen<1) {
    GMessage("Warning: zero length key found\n");
    return 0;
    }
 //------------ adding record -----------------
 if (addbegin(klen,datalen)==-1)
     GError("GCdbWrite: Error at addbegin(%d, %d)\n",klen, datalen);
 uint32 h=CDB_HASHSTART;
 for (i = 0;i < klen; ++i) {
      //if (cdbuffer_PUTC(c.cdbuf,key[i]) == -1)
      if ( ((cdbuf->n!=cdbuf->p) ? (cdbuf->x[cdbuf->p++]=(key[i]),0 )
                      : cdbuf->put(&(key[i]),1) )==-1)
               GError("GCdbWrite: Error at cdbbuf.put, key '%s'\n", key);
      h = cdb_hashadd(h,key[i]);
      }
 if (cdbuf->put(recdata,datalen) == -1)
    GError("GCdbWrite: Error at final cdbuf.put() at key='%s', datalen=%d\n",
                    key, datalen);
 if (addend(klen,datalen,h) == -1)
    GError("GCdbWrite: Error at addend(%d, %d, h)\n", klen, datalen);
 return 1;
}

int GCdbWrite::addrec(char *key,unsigned int keylen,char *data,unsigned int datalen) {
  if (GCdbWrite::addbegin(keylen,datalen) == -1) return -1;
  if (cdbuf->putalign(key,keylen) == -1) return -1;
  if (cdbuf->putalign(data,datalen) == -1) return -1;
  return GCdbWrite::addend(keylen,datalen,cdb_hash(key,keylen));
}


int GCdbWrite::finish() {
  char buf[8];
  int i;
  uint32 len;
  uint32 u;
  uint32 memsize;
  uint32 icount;
  uint32 where;
  struct cdb_hplist *x;
  struct cdb_hp *hp;

  for (i = 0;i < 256;++i)
    count[i] = 0;

  for (x = head;x;x = x->next) {
    i = x->num;
    while (i--)
      ++count[255 & x->hp[i].h];
  }

  memsize = 1;
  for (i = 0;i < 256;++i) {
    u = count[i] * 2;
    if (u > memsize)
      memsize = u;
  }

  memsize += numentries; /* no overflow possible up to now */
  u = (uint32) 0 - (uint32) 1;
  u /= sizeof(struct cdb_hp);
  if (memsize > u) { /* errno = error_nomem;*/ return -1; }

  split = (struct cdb_hp *) gcdb_alloc(memsize * sizeof(struct cdb_hp));
  if (!split) return -1;

  hash = split + numentries;

  u = 0;
  for (i = 0;i < 256;++i) {
    u += count[i]; /* bounded by numentries, so no overflow */
    start[i] = u;
  }

  for (x = head;x;x = x->next) {
    i = x->num;
    while (i--)
      split[--start[255 & x->hp[i].h]] = x->hp[i];
  }

  for (i = 0;i < 256;++i) {
    icount = count[i];

    len = icount + icount; /* no overflow possible */
    uint32_pack(final + 8 * i,pos);
    uint32_pack(final + 8 * i + 4,len);

    for (u = 0;u < len;++u)
      hash[u].h = hash[u].p = 0;

    hp = split + start[i];
    for (u = 0;u < icount;++u) {
      where = (hp->h >> 8) % len;
      while (hash[where].p)
	if (++where == len)
	  where = 0;
      hash[where] = *hp++;
    }

    for (u = 0;u < len;++u) {
      uint32_pack(buf,hash[u].h);
      uint32_pack(buf + 4,hash[u].p);
      if (cdbuf->putalign(buf,8) == -1) return -1;
      if (posplus(8) == -1) return -1;
    }
  }

  if (cdbuf->flush() == -1) return -1;
  if (gcdb_seek_begin(fd) == -1) return -1;
  return cdbuf->putflush(final,sizeof final);
}

//=====================================================
//-------------        cdb          -------------------
//=====================================================
uint32 cdb_hashadd(uint32 h,unsigned char c) {
  h += (h << 5);
  return h ^ c;
}

uint32 cdb_hash(char *buf,unsigned int len) {
  uint32 h;
  h = CDB_HASHSTART;
  while (len) {
    h = cdb_hashadd(h,*buf++);
    --len;
    }
  return h;
}

//---------------------------------------------------------------
//-------------------------- cdb methods ------------------------

GCdbRead::GCdbRead(int afd) {
  struct stat st;
  char *x;
 //check endianness :)
  gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;
  gcvt_offt=(endian_test())? &offt_sun : &offt_x86;

  findstart();
  fd = afd;
  if (fstat(fd,&st) == 0)

    if (st.st_size <= MAX_UINT) {
      x = (char *) mmap(0,st.st_size,PROT_READ,MAP_SHARED,fd,0);
      if (x + 1) {
        size = st.st_size;
        map = x;
        }
       else {
         GError("Error mapping the file (size=%ld)!\n",st.st_size);
         }
       }
    else {
       GError("Error mapping the file (size %ld > MAX_UINT)\n",
           st.st_size);
       }
}

GCdbRead::GCdbRead(char* afname) {
  struct stat st;
  char *x;
  
  //check endianness :)
  gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;
  gcvt_offt=(endian_test())? &offt_sun : &offt_x86;

  
  findstart();
  #ifdef __WIN32__
    fd = open(afname, O_RDONLY|O_BINARY);
  #else
    fd = open(afname, O_RDONLY);
  #endif  
  if (fd == -1)
     GError("Error: cannot open file %s\n", afname);
  strcpy(fname, afname);
  if (fstat(fd,&st) == 0)
    if (st.st_size <= MAX_UINT) {
      x = (char *) mmap(0,st.st_size,PROT_READ,MAP_SHARED,fd,0);
      if (x + 1) {
        size = st.st_size;
        map = x;
        }
       else {
         GError("GCdbRead: Error mapping the file (size=%ld)!\n",st.st_size);
         }
      }
    else {
       GError("GCdbRead: Error mapping the file (size %ld > MAX_UINT)\n",
           st.st_size);
       }
}


GCdbRead::~GCdbRead() {
  if (map) {
    munmap(map,size);
    map = 0;
    }
}

int GCdbRead::read(char *buf,unsigned int len, uint32 pos) {
  if (map) {
    if ((pos > size) || (size - pos < len)) {
          /* errno = error_proto; */
          return -1;
          }
    gcdb_byte_copy(buf, len, map + pos);
    }
  else {
    if (gcdb_seek_set(fd,pos) == -1) return -1;
    while (len > 0) {
      int r;
      do {
        r = ::read(fd,buf,len);
        } while ((r == -1) && (errno == error_intr));
      if (r == -1) return -1;
      if (r == 0) {
          //errno = error_proto;
          return -1;
          }
      buf += r;
      len -= r;
    }
  }
  return 0;
}

int GCdbRead::match(char *key, unsigned int len, uint32 pos) {
  char buf[32];
  unsigned int n;
  while (len > 0) {
    n = sizeof buf;
    if (n > len) n = len;
    if (GCdbRead::read(buf,n,pos) == -1) return -1;
    if (byte_diff(buf,n,key)) return 0;
    pos += n;
    key += n;
    len -= n;
  }
  return 1;
}

int GCdbRead::findnext(char *key,unsigned int len) {
  char buf[8];
  uint32 pos;
  uint32 u;
  if (!loop) {
    u = cdb_hash(key,len);
    if (GCdbRead::read(buf,8,(u << 3) & 2047) == -1) return -1;
    uint32_unpack(buf + 4,&hslots);
    if (!hslots) return 0;
    uint32_unpack(buf,&hpos);
    khash = u;
    u >>= 8;
    u %= hslots;
    u <<= 3;
    kpos = hpos + u;
    }
  while (loop < hslots) {
    if (GCdbRead::read(buf,8,kpos) == -1) return - 1;
    uint32_unpack(buf + 4, &pos);
    if (!pos) return 0;
    loop += 1;
    kpos += 8;
    if (kpos == hpos + (hslots << 3)) kpos = hpos;
    uint32_unpack(buf,&u);
    if (u == khash) {
      if (GCdbRead::read(buf,8,pos) == -1) return -1;
      uint32_unpack(buf,&u);
      if (u == len)
        switch(GCdbRead::match(key,len,pos + 8)) {
          case -1:
            return -1;
          case 1:
            uint32_unpack(buf + 4,&dlen);
            dpos = pos + 8 + len;
            return 1;
        }
    }
  }
  return 0;
}

int GCdbRead::find(char *key) {
  GCdbRead::findstart();
  return GCdbRead::findnext(key,gcdb_strlen(key));
}
