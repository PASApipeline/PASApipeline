#include "gcl/GCdbYank.h"
#include "gcl/GBase.h"
#include <ctype.h>

#define ERR_READ "cdbyank: error reading from file.\n"
#define ERR_READFMT "cdbyank read error: incorrect file format.\n"
#define ERR_RANGEFMT "Sequence range parsing error for key '%s'\n"
#define ERR_RANGE_INVALID "Invalid range (%d-%d) specified for sequence '%s' of length %d\n"
// 1MB memory buffer:
#define MAX_MEM_RECSIZE 1048576 
#ifndef O_BINARY
 #define O_BINARY 0x0000
#endif

GCdbZFasta::GCdbZFasta(FILE* azf, int zrsize, char* r_delim) {
 zrecsize=-1;
 zpos=0;
 recdelim=Gstrdup(r_delim);
 zf=azf;
 decomp_start(zrsize);
 chrhandler=new GFastaCharHandler(recdelim);
}

GCdbZFasta::~GCdbZFasta() {
 //if (zf!=NULL && zf!=stdout && zf!=stdin) fclose(zf);
 // FULL_FLUSH method instead of finish
 delete chrhandler;
 GFREE(recdelim);
 decomp_end();
}

void GCdbZFasta::decomp_start(int zrsize) {
 zstream.zalloc = (alloc_func)0;
 zstream.zfree = (free_func)0;
 zstream.opaque = (voidpf)0;
 zstream.next_in  = (Bytef*)sbuf;
 zstream.avail_in = 0;
 zstream.next_out = (Bytef*)lbuf;
 int err = inflateInit(&zstream);
 if (err!=Z_OK)
     GMessage("Error at inflateInit()\n");
 //-- now read and discard the first record, so we can use random access later
 // (needed by zlib)
 int bytes_read=fread(sbuf, 1, zrsize, zf);
 if (bytes_read<zrsize)
     GError("Error reading 1st record from zrec file\n");
 zstream.next_in = (Bytef*)sbuf;
 zstream.avail_in = bytes_read;
//decompress first chunk
 zstream.next_out = (Bytef*)lbuf;
 zstream.avail_out = GCDBZ_LBUF_LEN;
 err = inflate(&zstream, Z_SYNC_FLUSH);
 if (err !=Z_OK && err!=Z_STREAM_END)
     GError("GCdbZFasta error: 1st record inflate failed! (err=%d)\n",err);     
}

void GCdbZFasta::decomp_end() {
  int err = inflateEnd(&zstream);
  if (err!=Z_OK)
     GError("Error at inflateEnd() (err=%d)\n", err);
}


//fasta record decompress
//returns: the number of bytes decompressed
int GCdbZFasta::decompress(FastaSeq& rec, int csize, int zfofs, charFunc* seqCallBack) {
 if (zfofs>=0) {
    if (fseek(zf, zfofs, 0))
      GError("GCdbZFasta::decompress: error fseek() to %d\n", zfofs);
    }
  else
     if (feof(zf)) return 0;
 bool in_rec=true;
 int err=0;
 int total_read=0;
 int total_written=0;
 chrhandler->init(&rec, seqCallBack);
 while (in_rec) { // main read loop
     int to_read=0;
     int bytes_read=0;
     if (csize<=0) { //read one byte at a time
        to_read=1;
        int c;
        if ((c =fgetc(zf))!=EOF) {
           bytes_read = 1;
           sbuf[0]=c;
           }
          else {
            //bytes_read=0;
            return 0; //eof
            }
        total_read+=bytes_read;
        }
      else {
        to_read = csize-total_read>GCDBZ_SBUF_LEN ?
                                 GCDBZ_SBUF_LEN : csize-total_read;
       // check for csize vs bytes_read match:
        if (to_read==0) return 0;
        bytes_read=fread(sbuf, 1, to_read, zf);
        if (bytes_read!=to_read)
            GError("Error reading from zrec file\n");
        total_read+=bytes_read;
        in_rec=(total_read<csize);
        }
     if (bytes_read==0) {
        //GMessage("bytes_read = 0\n");
        return 0;
        }
     if (in_rec && bytes_read<to_read) in_rec=false;
     zstream.next_in = (Bytef*)sbuf;
     zstream.avail_in = bytes_read;

     do { //decompression loop
        zstream.next_out = (Bytef*)lbuf;
        zstream.avail_out = GCDBZ_LBUF_LEN;
        uLong t_out=zstream.total_out;
        err = inflate(&zstream, Z_SYNC_FLUSH);
        uLong toWrite=zstream.total_out-t_out;
        if (toWrite>0) {
             /* if (fwrite(lbuf, 1, toWrite, outf)<toWrite) {
               GError("Error writing inflated chunk!\n");
               } */
             for (unsigned int i=0;i<toWrite;i++)
               chrhandler->processChar(lbuf[i]);

             total_written+=toWrite;
             }
        if (err==Z_STREAM_END) {
              in_rec=false;
              if (total_written==0) {
                GMessage("Z_STREAM_END found but total_written=0!\n");
                }
              break;
              }
         else if (err !=Z_OK)
                GError("GCdbZFasta error: inflate failed! (err=%d)\n",err);
        } while (zstream.avail_in!=0); //decompression loop
   } //read loop
  chrhandler->done(); 
 /*if (err!=Z_STREAM_END) {
   GError("decompress: Z_STREAM_END not found!\n");
   }*/
  return total_written;
}





void inplace_Lower(char* c) {
 char *p=c;
 while (*p!='\0') { *p=tolower(*p);p++; }
}

void buf_get(GCDBuffer* b, uint32& pos, char *buf, unsigned int len) {
  int r;
  while (len > 0) {
    r = b->get(buf,len);
    if (r == -1) GError(ERR_READ);
    if (r == 0)
       GError(ERR_READFMT);
    pos += r;
    buf += r;
    len -= r;
    }
}

void buf_getnum(GCDBuffer* b, uint32& pos, uint32 *num) {
  char buf[4];
  buf_get(b, pos, buf, 4);
  uint32_unpack(buf,num);
}


int read_dbinfo(int fd, char** fnameptr, cdbInfo& dbstat) {
//this is messy due to the need of compatibility with the
//old 32bit file-length
       char* dbname=*fnameptr;
       //read just the tag first: 4 bytes ID
       lseek(fd, -cdbInfoSIZE, SEEK_END);
       int r=read(fd, &dbstat, cdbInfoSIZE );
       if (r!=cdbInfoSIZE) return 2;
       //GMessage("Size of dbstat=%d\n", cdbInfoSIZE);
       if (strncmp(dbstat.oldtag, "CIDX", 4)==0) {
            //old dbstat structure -- convert it
            dbstat.num_keys=gcvt_uint(&dbstat.oldnum[0]);
            dbstat.num_records=gcvt_uint(&dbstat.oldnum[1]);
            dbstat.dbsize=gcvt_uint(&dbstat.old_dbsize);
            dbstat.idxflags = gcvt_uint(&dbstat.old_idxflags);
             //position on the dbnamelen entry
            dbstat.dbnamelen = gcvt_uint(&dbstat.old_dbnamelen);
            //GMessage("dbnamelen=%d\n", dbstat.dbnamelen);
            lseek(fd, -(off_t)(cdbInfoSIZE-4+dbstat.dbnamelen), SEEK_END);
            }
          else if (strncmp(dbstat.tag, "CDBX", 4)!=0) {
            GMessage("Error: this doesn't appear to be a cdbfasta created file!\n");
            return 1;
            }
           else { // new CDBX type:
            dbstat.dbsize = gcvt_offt(&dbstat.dbsize);
            dbstat.num_keys=gcvt_uint(&dbstat.num_keys);
            dbstat.num_records=gcvt_uint(&dbstat.num_records);
            dbstat.idxflags = gcvt_uint(&dbstat.idxflags);
            //position on the dbnamelen entry
            dbstat.dbnamelen = gcvt_uint(&dbstat.dbnamelen);
            //GMessage("dbnamelen=%d\n", dbstat.dbnamelen);
            lseek(fd, -(off_t)(cdbInfoSIZE+dbstat.dbnamelen), SEEK_END);
            }

       GMALLOC(dbname, dbstat.dbnamelen+1);
       dbname[dbstat.dbnamelen]='\0';
       r=read(fd, dbname, dbstat.dbnamelen);
       *fnameptr=dbname;
       if (r!=dbstat.dbnamelen)
         return 2;
       return 0;
       }

int parse_int(FILE* f, char* buf, char* key, int& e) {
   char* p, *q;
   while (e!=EOF && isspace(e)) { //skip any spaces
     if (e=='\n') GError(ERR_RANGEFMT, key);
     e=fgetc(stdin);
     }
   if (e==EOF) GError(ERR_RANGEFMT, key);
   //now e is the first non-space
   p=buf;
   q=p;
   while (e!=EOF && !isspace(e)) {
    *q=e;
    q++;
    e=fgetc(stdin);
    }
   *q='\0'; //now p is the starting coordinate string
   return atoi(p);
   //now the file pointer should be on the first space after the parsed value
}

int parse_int(char*& f, char* key, int& e) {
   char* p, *q;
   char buf[16];
   while (e!='\0' && isspace(e)) { //skip any spaces
     if (e=='\n') GError(ERR_RANGEFMT, key);
     f++;
     e=*f;
     } 
   if (e=='\0') GError(ERR_RANGEFMT, key);
   //now e is the first non-space char
   p=buf;
   q=p;
   while (e!='\0' && !isspace(e)) {
     *q=e;
     q++;
     f++;
     e=*f;     
     }
   *q='\0';
   return atoi(p);
   //now f and e should be on the first space after the parsed value (or '\0')
}

GCdbZFasta* GCdbYank::openCdbz(char* p) {
    //in case this was not done before:
    gcvt_uint=(endian_test())? &uint32_sun : &uint32_x86;
    FILE* zf=fopen(p, "rb");
    if (zf==NULL) {
          GMessage("Error: cannot open compressed file '%s'!\n",p);
          return NULL;
          }
    //check if the file is valid and read the length of the first record
    //
    char ztag[5];
    ztag[4]='\0';
    if (fread(ztag, 1, 4, zf)<4) {
       GMessage("Error reading header of compressed file '%s'\n",p);
       return NULL;
       }
    if (strcmp(ztag, "CDBZ")!=0) {
       GMessage("Error: file '%s' doesn't appear to be a zlib compressed cdb?\n",p);
       return NULL;
       }
    unsigned int zrecsize;
    if (fread((void*) &zrecsize,1,4,zf)<4) {
       GMessage("Error reading 1st compressed record size for file '%s'!\n",p);
       return NULL;
       }
   zrecsize=gcvt_uint(&zrecsize);
   return new GCdbZFasta(zf, zrecsize, recdelim);
}

GCdbYank::GCdbYank(const char* fidx, char* recsep) {
 is_compressed=false;
 fd=-1;
 cdb=NULL;
 cdbz=NULL;
 fdb=-1;
 fz=NULL;
 dbname=NULL;
 recdelim=Gstrdup(recsep);
 if (fidx==NULL) GError("GCdbYank Error: NULL index file name!");
 idxfile=Gstrdup(fidx);
 cdb=new GCdbRead(idxfile);
 fd=cdb->getfd();
 db_size=0;
 dbstat.dbsize=0;
 info_dbname=NULL;
 int r=read_dbinfo(fd, &info_dbname, dbstat);
 lseek(fd, 0, SEEK_SET);
 if (r==1) GError("This file does not seem to be a cdbfasta generated file.\n");
          else if (r==2)
                 GError("Error reading info chunk!\n");
  /*try to find the database file
     rules: if given, only the -d given filename is used
       otherwise:
        1) the same directory with the given index file(stripping the suffix)
        2) the dbstat filepath/name stored by cdbfasta
   */
 if (dbname==NULL) {
   char* p = rstrchr(idxfile, '.');
   if (p!=NULL) {
      /*GError("%s\ncdbyank error: cannot use %s as an index file. When no -d is\n\
      given, so the database file can be located in the same directory \n\
      by removing the index file suffix (.cidx)\n", USAGE, idxfile);*/
      int nlen=p-idxfile;
      char* namebuf=NULL;
      GMALLOC(namebuf, nlen+1);
      strncpy(namebuf, idxfile, nlen);
      namebuf[nlen]='\0';
      if (fileExists(namebuf))
         dbname=namebuf;

      }  // strip the index file extenstion
    // 2) try the stored dbstat name
   if (dbname==NULL) {
       if (fileExists(info_dbname)) dbname=info_dbname;
         else GError("Cannot locate the database file for this index\n");
      }
  }// database name not given
  is_compressed=(dbstat.idxflags & CDBMSK_OPT_COMPRESS);
       if (is_compressed)
            //try to open the dbname as a compressed file
             fz=fopen(dbname, "rb");
       else  fdb=open(dbname, O_RDONLY|O_BINARY);
     if (fdb==-1 && fz==NULL)
           GError("Error: cannot open database file %s\n",dbname);
     if (is_compressed) {
        fclose(fz);//just to start fresh here
        //determine size:
        int ftmp = open(dbname, O_RDONLY|O_BINARY);
        if (ftmp == -1) GError("Error reopening db '%s'?\n",dbname);
        struct stat fdbstat;
        fstat(ftmp, &fdbstat);
        db_size=fdbstat.st_size;
        close(ftmp);
        //-------- reopen here
        cdbz=openCdbz(dbname);
        if (cdbz==NULL)
           GError("Error opening the cdbz file '%s'\n");
        fz=cdbz->getZFile();
        }
       else {
        struct stat fdbstat;
        if (stat(dbname, &fdbstat)!=0) {
          perror("stat()");
          exit(1);
          }
        db_size=fdbstat.st_size;
        }
     //abort if the database size was read and it doesn't match the cdbfasta stored size
     if (dbstat.dbsize>0 && dbstat.dbsize!=db_size)
       GError("Error: invalid %d database size - (%lld vs %lld) please rerun cdbfasta for '%s'\n",
          fdb, dbstat.dbsize, db_size, dbname);
   fastahandler=new GFastaCharHandler(recdelim);
} //* GCdbYank constructor *//


GCdbYank::~GCdbYank() {
 if (is_compressed) {
     fclose(fz); delete cdbz;
     }
     else close(fdb);
 GFREE(info_dbname);
 delete fastahandler;
 GFREE(recdelim);
 GFREE(dbname);
 GFREE(idxfile);
 delete cdb;
 close(fd);
}


int GCdbYank::getRecord(char* key, FastaSeq& rec, charFunc* seqCallBack) {
//assumes fdb is open, cdb was created on the index file
 int r=cdb->find(key);
 if (r==0) return 0;
 if (r==-1)
   GError("cdbyank: error searching for key %s in %s\n", key, idxfile);
 /* while (r>0) { */
 off_t pos = cdb->datapos(); //position of this key's record in the index file
 unsigned int len=cdb->datalen(); // length of this key's record
 char bytes[32]; // data buffer -- should just accomodate fastarec_pos, fastarec_length
 if (cdb->read(bytes,len,pos) == -1)
       GError("cdbyank: error at GCbd::read (%s)!\n", idxfile);

 off_t fpos; //this will be the fastadb offset
 uint32 reclen;  //this will be the fasta record offset
 if (len>8) { //64 bit file offset was used
  fpos=gcvt_offt(bytes);
  reclen=gcvt_uint(&bytes[sizeof(uint32)<<1]);
  }
  else { //32bit offset used
  fpos=gcvt_uint(bytes);
  reclen=gcvt_uint(&bytes[sizeof(uint32)]);
  }
  //GMessage("reclen=%d\n", reclen);

 /* if (showQuery)
  fprintf(fout, "%c%s%c\t", delimQuery, key, delimQuery);*/
  /*========= FETCHING RECORD CONTENT ======= */
   if (is_compressed) {
     //for now: ignore special retrievals, just print the whole record
     return cdbz->decompress(rec, reclen, fpos, seqCallBack);
     }
    else { // not compressed -- position into the file and build an ad hoc GFastaFile
     lseek(fdb, fpos, SEEK_SET);
     // read it char by char and return it as output
     char c='\0';
     int charsread=0;
     fastahandler->init(&rec, seqCallBack);
     while (reclen-- && read(fdb, &c, 1)==1) {
       fastahandler->processChar(c);
       charsread++;
       }
     fastahandler->done();
     return charsread;
     } // not compressed
  /* if (many) r=cdb->findnext(key, strlen(key));
        else r=0;
   } */
}

off_t GCdbYank::getRecordPos(char* key) {
//assumes fdb is open, cdb was created on the index file
 int r=cdb->find(key);
 if (r==0 && warnings) {
   GMessage("cdbyank: key \"%s\" not found in %s\n", key, idxfile);
   return 0;
   }
 if (r==-1)
   GError("cdbyank: error searching for key %s in %s\n", key, idxfile);
 off_t pos = cdb->datapos(); //position of this key's record in the index file
 unsigned int len=cdb->datalen(); // length of this key's record
 char bytes[32]; // data buffer -- should just accomodate fastarec_pos, fastarec_length
 if (cdb->read(bytes,len,pos) == -1)
       GError("cdbyank: error at GCbd::read (%s)!\n", idxfile);

 off_t fpos; //this will be the fastadb offset
 //uint32 reclen;  this will be the fasta record length
 if (len>8) //64 bit file offset was used
      fpos=gcvt_offt(bytes);
     else  //32bit offset used
      fpos=gcvt_uint(bytes);
  return fpos;
}
