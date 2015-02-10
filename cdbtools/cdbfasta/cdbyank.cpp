#include "gcl/gcdb.h"
#include "gcl/GBase.h"
#include "gcl/GArgs.h"
#include "ctype.h"
#include <fcntl.h>
#include <string.h>
#include "gcdbz.h"
#ifdef __WIN32__
#define VERSION "cdbyank version 0.981w"
#else
#define VERSION "cdbyank version 0.981"
#endif
#define USAGE "Usage:\n\
  cdbyank <index_file> [-d <fasta_file>] [-a <key>|-n|-l|-s]\n\
      [-o <outfile>] [-q <char>|-Q][-F] [-R] [-P] [-x] [-w] \n\
      [-z <dbfasta.cdbz>\n\n\
    <index_file> is the index file created previously with cdbfasta\n\
       (usually having a \".cidx\" suffix)\n\
    -a <key> the sequence name (accession) for a fasta record to be\n\
       retrieved; if not given, a list of accessions is expected\n\
       at stdin\n\
    -d <fasta_file> is the fasta file to pull records from; \n\
       if not specified, cdbyank will look in the same directory\n\
       where <index_file> resides, for a file with the same name\n\
       but without the \".cidx\" suffix\n\
    -o the records found are written to file <outfile> instead of stdout\n\
    -x allows retrieval of multiple records per key, if the indexed \n\
       database had records with the same key (non-unique keys);\n\
       (without -x only one record for a given key is retrieved)\n\
    -i case insensitive query (expects the <index_file> to have been \n\
       created with cdbfasta -i option)\n\
    -Q output the query key surrounded by character '%' before the\n\
       corresponding record\n\
    -q same as -Q but use character <char> instead of '%'\n\
    -w enable warnings (sent to stderr) when a key is not found\n\
    -F pulls only the defline for each record (discard the sequence)\n\
    -P only displays the position(s) (file offset) within the \n\
       database file, for the requested record(s)\n\
    -R sequence range extraction: expects the input <key(s)> to have \n\
       the format: '<seq_name> <start> <end>'\n\
       and pulls only the specified sequence range\n\
    -z decompress the entire file <dbfasta.cdbz>\n\
       (assumes it was built using cdbfasta with '-z' option)\n\
    -v show version number and exit\n\
    \n\
    Index file statistics (no database file needed):\n\
    -n display the number of records indexed\n\
    -l list all keys stored in <index_file>\n\
    -s display indexing summary info\n\n"

/*
    -E same as -R but assumes FASTA records have a fixed line length\n\
       (faster extraction of distant ranges for long records)\n\
*/
#define ERR_READ "cdbyank: error reading from file.\n"
#define ERR_READFMT "cdbyank read error: incorrect file format.\n"
#define ERR_RANGEFMT "Sequence range parsing error for key '%s'\n"
#define ERR_RANGE_INVALID "Invalid range (%d-%d) specified for sequence '%s' of length %d\n"
// 1MB memory buffer:
#define MAX_MEM_RECSIZE 1048576 
#ifndef O_BINARY
 #define O_BINARY 0x0000
#endif


static char* idxfile;
static int warnings;
bool is_compressed=false;
bool defline_only=false;
bool rec_pos_only=false;
bool use_range=false;
bool fixed_linelen=false;
bool caseInsensitive=false;
bool showQuery=false;
char delimQuery='%';

off_t lastfpos=-1; //to avoid pulling the same record twice in a row..

FILE* fout=NULL;
GCdbRead* cdb=NULL;

GCdbz* cdbz=NULL;
int fdb=-1;
FILE* fz=NULL;

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


int fetch_record(char* key, char* dbname, int many, int r_start=0, int r_end=0) {
//assumes fdb is open, cdb was created on the index file
 if (caseInsensitive) inplace_Lower(key);
 int r=cdb->find(key);
 if (r==0 && warnings) {
   GMessage("cdbyank: key \"%s\" not found in %s\n", key, idxfile);
   return 0; 
   }
 if (r==-1)
   GError("cdbyank: error searching for key %s in %s\n", key, idxfile);
 while (r>0) {
   off_t pos = cdb->datapos(); //position of this key's record in the index file
   unsigned int len=cdb->datalen(); // length of this key's record 
   char bytes[32]; // data buffer -- should just accomodate fastarec_pos, fastarec_length
   if (cdb->read(bytes,len,pos) == -1)
       GError("cdbyank: error at GCbd::read (%s)!\n", idxfile);
       
   off_t fpos; //this will be the fastadb offset
   uint32 reclen;  //this will be the fasta record offset
   if (len>8) { //64 bit file offset was used
    fpos=gcvt_offt(bytes);
    if (rec_pos_only) {
      fprintf(fout, "%lld\n", fpos);
      return 1;
      }
    reclen=gcvt_uint(&bytes[sizeof(uint32)<<1]);
    }
    else { //32bit offset used
    fpos=gcvt_uint(bytes);
    if (rec_pos_only) {
      fprintf(fout, "%lld\n", fpos);
      return 1;
      }
    reclen=gcvt_uint(&bytes[sizeof(uint32)]);
    }
    //GMessage("reclen=%d\n", reclen);


   if (fpos == lastfpos) {
      if (many) r=cdb->findnext(key, strlen(key));
           else r=0;
      continue;
      }
   lastfpos=fpos;
   if (showQuery)
    fprintf(fout, "%c%s%c\t", delimQuery, key, delimQuery);
   if (is_compressed) {
     //for now: ignore special retrievals, just print the whole record
     cdbz->decompress(fout, reclen, fpos);
     if (many) r=cdb->findnext(key, strlen(key));
         else r=0;
     continue;
     }
   lseek(fdb, fpos, SEEK_SET);
   if (reclen<=MAX_MEM_RECSIZE) {
       char* p;
       GMALLOC(p,reclen+1);
       //errno=0;
       r=read(fdb, p, reclen);
       if (r<=0)
          GError("cdbyank: Error reading from database file [%s] for %s (returned %d, offset %d) !\n",
                  dbname, idxfile, r, fpos);
       p[reclen]='\0';
       //--- now we have the whole record, check if some special options were given:
       if (defline_only) {
         char* q=strchr(p,'\n');
         if (q!=NULL) *q='\0';
         //skip '>' char
         fprintf(fout, "%s\n",p+1);
         }
        else
         if (use_range && r_start + r_end >2) { //range case
           //extract only a substring of the sequence
           char* r=strchr(p,'\n');
           if (r!=NULL) *r='\0'; //now p only has the defline
           fprintf(fout, "%s\n", p); //output the defline
           r++;
           unsigned int recpos=r-p; //p[recpos] MUST be a nucleotide or aminoacid now!
           int seqpos=0;
           char linebuf[61];
           int linelen=0;
           while (recpos<reclen) {
              if (isspace(p[recpos])) recpos++; //skip newlines, etc. in the fasta sequence
                 else {
                    seqpos++;
                    if (seqpos>=r_start && seqpos<=r_end) {
                      linebuf[linelen]=p[recpos];
                      linelen++;
                      if (linelen==60 || seqpos==r_end) {
                         linebuf[linelen]='\0';
                         linelen=0;
                         fprintf(fout, "%s\n", linebuf);
                         }
                      }
                    recpos++;
                    }
              }
           }
        else { //not range display
         fprintf(fout, "%s\n",p);
         }
       GFREE(p);
      } //small record
    else { //large record, read it char by char and return it as output
     char c='\0';
     if (defline_only) {
     	  reclen--;
     	  read(fdb, &c, 1);
     	  }
     while (reclen-- && read(fdb, &c, 1)==1) {
       fprintf(fout, "%c", c);
       if (c=='\n') break;
       }
     //defline written
     if (!defline_only) {
      int seqpos=1;
      if (use_range) {
         while (reclen-- && read(fdb, &c, 1)==1 && seqpos<=r_end) {
           if (isspace(c)) continue;
           if (seqpos>=r_start) {
              int written=seqpos-r_start;
              if (written && written%60 == 0)
                   fprintf(fout,"\n");
              fprintf(fout, "%c", c);
              }
           seqpos++;
           }//while
         } //range case
       else { //no range, just copy all chars to output
         while (reclen-- && read(fdb, &c, 1)==1) {
            fprintf(fout, "%c", c);
            }
         }
       fprintf(fout, "\n");
       }
     }
   if (many) r=cdb->findnext(key, strlen(key));
        else r=0;
   }
 return 1;
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

GCdbz* openCdbz(char* p) {
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
   return new GCdbz(zf, true, zrecsize);
}


int main(int argc, char **argv) {
  char namebuf[1024];
  int r_start, r_end;
  char* p;
  char* dbname=NULL;
  int result=0;
  int r=0;
  cdbInfo dbstat;
  dbstat.dbsize=0;
  GArgs args(argc, argv, "a:d:o:z:q:nlsxwvFREiPQ");
  int e=args.isError();
  if (e>0)
     GError("%s Invalid argument: %s\n", USAGE, argv[e]);
  if (args.getOpt('v')!=NULL) {
    printf("%s\n",VERSION); 
    return 0;
    }
  char* outfile=(char*)args.getOpt('o');
  if (outfile!=NULL) {
      if ((fout=fopen(outfile, "wb"))==NULL)
         GError("Cannot create file '%s'!", outfile);
      }
    else fout=stdout;  

  if ((p=(char*)args.getOpt('z'))!=NULL) { //simply stream-decompress cdbz
    
    GCdbz* cdbz=openCdbz(p);
    if (cdbz==NULL)
       GError("Error opening the cdbz file '%s'\n");
    FILE* zf=cdbz->getZFile();
    int numrecs=0;
    int xcode;
    while ((xcode=cdbz->decompress(fout))>0) numrecs++;
    delete cdbz;
    fclose(zf);
    return 0;
    }
  int numfiles = args.startNonOpt();
  if (numfiles==0)
    GError("%s Error: an index file must be provided !\n", USAGE);
  idxfile=(char*)args.nextNonOpt(); //first fasta file given
  char* key=(char*)args.getOpt('a');
  
  defline_only=(args.getOpt('F')!=NULL);
  rec_pos_only=(args.getOpt('P')!=NULL);
  showQuery=(args.getOpt('Q')!=NULL);
  const char* q;
  if ((q=args.getOpt('q'))!=NULL) {
   delimQuery=*q;
   showQuery=true;
   }
  use_range=((args.getOpt('R')!=NULL) || (args.getOpt('E')!=NULL));
  fixed_linelen=(args.getOpt('E')!=NULL);
  caseInsensitive=(args.getOpt('i')!=NULL);
  /*is_compressed=((args.getOpt('Z')!=NULL) || 
                (strstr(idxfile,".cidxz")!=NULL));*/
  int listQuery=(args.getOpt('l')!=NULL);
  warnings=(args.getOpt('w')!=NULL);
  int dataQuery=(!listQuery && args.getOpt('n')==NULL
           && args.getOpt('l')==NULL &&args.getOpt('s')==NULL);
        //exclude the possibility of index-only stats query
 dbname=(char*)args.getOpt('d');
 int fd;
 cdb=new GCdbRead(idxfile);
 fd=cdb->getfd();
 char* info_dbname=NULL;
 off_t db_size=0;
 dbstat.dbsize=0;

 r=read_dbinfo(fd, &info_dbname, dbstat);
 lseek(fd, 0, SEEK_SET);
 if (r==1) GError("This file does not seem to be a cdbfasta generated file.\n");
          else if (r==2)
                 GError("Error reading info chunk!\n");
 if (dataQuery) {
   //--------------- DB QUERY MODE: (always read the cdb stored info!)
   /*try to find the database file
     rules: if given, only the -d given filename is used
       otherwise:
        1) the same directory with the given index file(stripping the suffix)
        2) the dbstat filepath/name stored by cdbfasta
   */

   if (!rec_pos_only && dbname==NULL) { // no -d database given, find it
    // 1) try to rip the suffix:
    p = rstrchr(idxfile, '.');
    if (p!=NULL) {
     /*GError("%s\ncdbyank error: cannot use %s as an index file. When no -d is\n\
     given, so the database file can be located in the same directory \n\
     by removing the index file suffix (.cidx)\n", USAGE, idxfile);*/
     int nlen=p-idxfile;
     strncpy(namebuf, idxfile, nlen);
     namebuf[nlen]='\0';
     if (fileExists(namebuf))
        dbname=namebuf;
     }
    // 2) try the stored dbstat name
    if (dbname==NULL) {
      if (fileExists(info_dbname)) dbname=info_dbname;
       else GError("Cannot locate the database file for this index\n");
      }
    }
   if (!rec_pos_only) {
     if (!is_compressed) {
       if (r==0 && (dbstat.idxflags & CDBMSK_OPT_COMPRESS))
         is_compressed=true;
       }
     if (is_compressed)
            //try to open the dbname as a compressed file
             fz=fopen(dbname, "rb");
       else  fdb=open(dbname, O_RDONLY|O_BINARY);
     if (fdb==-1 && fz==NULL)
           GError("Error: cannot open database file %s\n",dbname);
     if (is_compressed) {
        fclose(fz);//just to start fresh here
        if (use_range)
           GError("Error: cannot use range extraction with compressed records, sorry.\n");
        if (defline_only)
          GError("Error: cannot use defline-only retrieval with compressed records (sorry).\n");
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
     }
   int many=(args.getOpt('x')!=NULL);
   int keypos=0;
   if (key==NULL) { //key not given
       GMALLOC(key, 2048);
       //get the keys at stdin
       if (use_range) {
           //expects the key and its sequence range on a single line!
         while ((e=fgetc(stdin)) != EOF) {
          if (isspace(e)) { //word end, close it
            key[keypos]='\0';
            if (keypos==0) continue;
            r_start=parse_int(stdin, &key[keypos+1], key, e);
            if (r_start<=0) GError(ERR_RANGEFMT, key);
            if (e==EOF || e=='\n') GError(ERR_RANGEFMT, key);
            r_end=parse_int(stdin, &key[keypos+1], key, e);
            if (r_end<=0 || r_end<=r_start) GError(ERR_RANGEFMT, key);
            fetch_record(key, dbname, many, r_start, r_end);
            //if (rec_pos_only) break;
            if (e==EOF) break;
            keypos=0;
            }
           else { //extend the key string
            key[keypos]=e;
            keypos++; 
            }
          } //while
         } //range case
       else { //no range, accept any space delimiter
         while ((e=fgetc(stdin)) != EOF) {
          if (isspace(e)) { //word end, close it
            key[keypos]='\0';
            fetch_record(key, dbname, many);
            //if (rec_pos_only) break;
            keypos=0;
            }
           else { //extend the key string
            key[keypos]=e;
            keypos++; 
            }
          } //while
       }
       GFREE(key);
      } //stdin case
    else { //key given already on command line
       //get only the first word of it:
       size_t keylen=strlen(key);
       p=key; while (!isspace(*p) && *p!='\0') p++;
       if (*p!='\0') *p='\0';
       if (use_range) {
           //parse the range from the query string
           if (keylen==strlen(p)) GError(ERR_RANGEFMT, key);
           p++;e=*p;
           r_start=parse_int(p, key, e);
           if (r_start<=0) GError(ERR_RANGEFMT, key);
           if (e=='\0' || e=='\n') GError(ERR_RANGEFMT, key);
           r_end=parse_int(p, key, e);
           if (r_end<=0 || r_end<=r_start) GError(ERR_RANGEFMT, key);
           }
         else {
          r_start=0;
          r_end=0;
          }
        if (fetch_record(key, dbname, many, r_start, r_end)==0)
               result=1; //the only key given not found
       }
    //end data query:
    if (!rec_pos_only) {
        if (is_compressed) {
         fclose(fz); delete cdbz;
         }
        else close(fdb);
       }
    if (fout!=NULL) fclose(fout);
    }
  //--------------- INDEX ONLY QUERY MODE:
  else { //index query mode: just retrieve some statistics or key names
    if (listQuery) { //request for list keys
       uint32 eod;
       uint32 pos=0;
       uint32 klen;
       uint32 dlen;
       char* bufspace;
       GMALLOC(bufspace, GCDBUFFER_INSIZE);
       GCDBuffer* readbuf=new GCDBuffer((opfunc)&read,
           fd, bufspace, GCDBUFFER_INSIZE);

       buf_getnum(readbuf, pos, &eod);
       GMALLOC(key, 1024); //!!! hopefully we don't have keys larger than that
       while (pos < 2048)
         buf_getnum(readbuf, pos, &dlen);
       while (pos < eod) {
          buf_getnum(readbuf, pos,&klen);
          buf_getnum(readbuf, pos,&dlen);
          //read key:
          buf_get(readbuf, pos, key, klen);
          key[klen]='\0';
          printf("%s\n", key);
          //read data (and ignore it)
          //assume that data is always shorter than 1K (should be just 4 bytes)
          buf_get(readbuf, pos, key, dlen);
          }
       GFREE(key);
       GFREE(bufspace);
       delete readbuf;
       }
     else { //dig up the info written at the end of the database file
       if (args.getOpt('n')!=NULL) {
          printf("%d\n",dbstat.num_records);
          }
        else {//must be -s
            printf("-= Indexing information: =-\n");
            printf("Number of records:%12d\n", dbstat.num_records);
            printf("Number of keys   :%12d\n", dbstat.num_keys);
            if (dbstat.idxflags & CDBMSK_OPT_COMPRESS)
                printf("Database records are compressed.\n");
            if (dbstat.idxflags & CDBMSK_OPT_MULTI)
                printf("Index was built with \"multi-key\" option enabled.\n");
            if (dbstat.idxflags & CDBMSK_OPT_C)
                printf("Index was built with \"shortcut keys\" only.\n");
               else if (dbstat.idxflags & CDBMSK_OPT_CADD)
                printf("The index was built with full keys and \"shortcut keys\".\n");
            printf("Database file: %s\n", info_dbname);
            printf("Database size: %lld bytes\n", dbstat.dbsize);
            }
       }
    }
 GFREE(info_dbname);
 delete cdb;
 close(fd);
 //getc(stdin);
 return result;
}
