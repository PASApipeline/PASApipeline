#include <ctype.h>
#include <fcntl.h>
#include <string.h>
#include <sys/stat.h>
#include "gcl/gcdb.h"
#include "gcl/GBase.h"
#include "gcl/GArgs.h"
#include "gcl/GReadBuf.h"
#include "gcl/GHash.hh"

#include "gcdbz.h"
#ifdef __WIN32__
#define VERSION "cdbfasta version 0.982w"
#else
#define VERSION "cdbfasta version 0.982"
#endif

#define USAGE "Usage:\n\
  cdbfasta <fastafile> [-o <index_file>] [-r <record_delimiter>]\n\
   [-z <compressed_db>] [-i] [-m|-n <numkeys>|-f<LIST>]|-c|-C]\n\
    [-w <stopwords_list>] [-s <stripendchars>] [-v]\n\
   \n\
   Creates an index file for records from a multi-fasta file.\n\
   By default (without -m/-n/-c/-C option), only the first \n\
   space-delimited token from the defline is used as a key.\n\
  \n\
   <fastafile> is the multi-fasta file to index; \n\
   -o the index file will be named <index_file>; if not given,\n\
      the index filename is database name plus the suffix '.cidx'\n\
   -r <record_delimiter> a string of characters at the beginning of line\n\
      marking the start of a record (default: '>')\n\
   -z database is compressed into the file <compressed_db>\n\
      before indexing (<fastafile> can be \"-\" or \"stdin\" \n\
      in order to get the input records from stdin)\n\
   -s strip extraneous characters from *around* the space delimited\n\
      tokens, for the multikey options below (-m,-n,-f);\n\
      Default <stripendchars> set is: '\",`.(){}/[]!:;~|><+-\n\
   -m (\"multi-key\" option) create hash entries pointing to \n\
      the same record for all tokens found in\n\
      the defline\n\
   -n <numkeys> same as -m, but only takes the first <numkeys>\n\
      tokens from the defline\n\
   -f indexes *space* delimited tokens (fields) in the defline as given\n\
      by LIST of fields or fields ranges (the same synatx as UNIX 'cut')\n\
   -w <stopwordslist> exclude from indexing all the words found\n\
      in the file <stopwordslist> (for options -m, -n and -k)\n\
   -i do case insensitive indexing (i.e. create additional keys for \n\
      all-lowercase tokens used for indexing from the defline \n\
   -c for deflines in the format: db1|accession1|db2|accession2|...,\n\
      only the first db-accession pair ('db1|accession1') is taken as key\n\
   -C like -c, but also subsequent db|accession constructs are indexed,\n\
      along with the full (default) token; additionally,\n\
      all nrdb concatenated accessions found in the defline \n\
      are parsed and stored (assuming 0x01 or '^|^' as separators)\n\
   -v show program version and exit\n"


#define ERR_TOOMANYFIELDS "Error: too many fields for -f option\n"
//16K initial defline buffer
#define KBUFSIZE 0x4000
#ifndef O_BINARY
 #define O_BINARY 0x0000
#endif
#define MAX_KEYLEN 1024
//64K input buffer
#define GREADBUF_SIZE 0x010000


typedef void (*addFuncType)(char*, off_t, uint32);

// for passing around index data:
struct IdxData32 {
   uint32 fpos;
   uint32 reclen;
   };

struct IdxData {
   off_t fpos; //32bit on Unix systems
   uint32 reclen;
   };

static int IdxDataSIZE=offsetof(IdxData, reclen)+sizeof(uint32);
static int IdxDataSIZE32=offsetof(IdxData32, reclen)+sizeof(uint32);

static char ftmp[365];
static char fztmp[365];
static char record_marker[127]; //record delimiter
static int  record_marker_len=1;
static int num_recs;
static int num_keys;

static int compact_plus; //shortcut key and

static bool do_compress=false; // compression used
FILE* zf=NULL; //compressed file handle

//store just offset and record length
char* defWordJunk="'\",`.(){}/[]!:;~|><+-";
static char* wordJunk;
static bool caseInsensitive=false; //case insensitive storage
static bool useStopWords=false;
static unsigned int numFields=0;
// we have fields[numFields-1]=MAX_UINT as defined in gcdb.h -- 
// as an indicator of taking every single token in the defline, 
// or for open ended ranges (e.g. -f5- )
static unsigned int fields[255]; //array of numFields field indices (1-based)
GHash<int> stopList;
//static int datalen=sizeof(uint32)+sizeof(off_t);

static char lastKey[MAX_KEYLEN]; //keep a copy of the last valid written key

GCdbWrite* cdbidx;
addFuncType addKeyFunc;

#define ERR_W_DBSTAT "Error writing the database statististics!\n"


void die_write(char* fname) {
  GError("Error: cdbhash was unable to write into file %s\n",
      fname);

 }

void die_read(char* infile) {
  GError("Error: cdbhash was unable to read the input file %s\n",infile);
}


void die_readformat(char* infile) {
  GError("Error: bad format for file %s; is it a fastA file?\n",
       infile);
}


bool add_cdbkey(char* key, off_t fpos, uint32 reclen) {

 unsigned int klen=strlen(key);
 if (klen<1) {
    GMessage("Warning: zero length key found following key '%s'\n",
              lastKey);
    return false;
    }
  //------------ adding record -----------------
 num_keys++;
 strncpy(lastKey, key, MAX_KEYLEN-1);
 lastKey[MAX_KEYLEN-1]='\0';
 if ((uint64)fpos>(uint64)MAX_UINT) { //64 bit file offset
  IdxData recdata;
  uint64 v= (uint64) fpos; //needed for Solaris' off_t issues with gcc/32
  recdata.fpos=gcvt_offt(&v);
  recdata.reclen=gcvt_uint(&reclen);  
  if (cdbidx->addrec(key,klen,(char*)&recdata,IdxDataSIZE)==-1)
    GError("Error adding cdb record with key '%s'\n",key);
  }
 else {//32 bit file offset
  IdxData32 recdata;
  uint32 v=(uint32) fpos;
  recdata.fpos=gcvt_uint(&v);
  recdata.reclen=gcvt_uint(&reclen);
  //GMessage("Adding 32bit: '%s' reclen=%d\n", key, recdata.reclen);
  if (cdbidx->addrec(key,klen,(char*)&recdata,IdxDataSIZE32)==-1)
    GError("Error adding cdb record with key '%s'\n",key);
  }

 return true;
}

//default indexing: key directly passed -- 
// as the first space delimited token
void addKey(char* key, off_t fpos,
             uint32 reclen) {
 num_recs++;
 add_cdbkey(key, fpos, reclen);
 if (caseInsensitive) {
   char* lckey=loCase(key);
   if (strcmp(lckey, key)!=0) 
      add_cdbkey(lckey, fpos, reclen);
   GFREE(lckey);
   }
}

//the whole defline is passed 
void addKeyMulti(char* defline,
                 off_t fpos, uint32 reclen) {
 char* p=defline;
 unsigned int fieldno=0;
 char* pn;
 num_recs++;
 bool stillParsing=true;
 unsigned int fidx=0; //index in fields[] array
 while (stillParsing) {
       while ((*p)==' ' || (*p)=='\t') p++;
       if (*p == '\0') break;
       //skip any extraneous characters at the beginning of the token
       while (chrInStr(*p, wordJunk)) p++;
       //skip any padding spaces or other extraneous characters 
       //at the beginning of the word
       if (*p == '\0') break;
       pn=p;
       while (*pn!='\0' && !isspace(*pn)) pn++; 
       //found next space or end of string
       fieldno++;
       while (fields[fidx]<fieldno && fidx<numFields-1) fidx++;
       //GMessage("> p=>'%s' [%d, %d, %d] (next='%s')\n",p, numFields, 
        //      fieldno, fields[numFields-1], pn+1);
       stillParsing = (((*pn)!='\0') && (fieldno+1<=fields[numFields-1]));
       char* pend = pn-1; //pend is on the last non-space in the current token
       //--strip the ending junk, if any
       while (chrInStr(*pend, wordJunk) && pend>p) pend--;
       if (pend<pn-1) *(pend+1)='\0';
                 else *pn='\0';
       
       if (strlen(p)>0) {
         if (fields[fidx]==MAX_UINT || fields[fidx]==fieldno) {
           if (useStopWords && stopList.hasKey(p)) {
             p=pn+1;
             continue;
             }
           //--- store this key with the same current record data:
           add_cdbkey(p, fpos, reclen);
           //---storage code ends here
           if (caseInsensitive) {
              char* lcp=loCase(p);
              if (strcmp(lcp,p)!=0)
                  add_cdbkey(lcp, fpos, reclen);
              GFREE(lcp);
              }
           }
         //if (isEnd) break; //if all the token were stored
         }
       p=pn+1;//p is now on the next token's start
       } //token parsing loop
}


static int qcmpInt(const void *p1, const void *p2) {
 //int n1=*((int*)p1);
 //int n2=*((int*)p2);
 const unsigned int *a = (unsigned int *) p1;
 const unsigned int *b = (unsigned int *) p2;
 if (*a < *b) return -1;
      else if (*a > *b) return 1;
                    else return 0;
 }


char* parse_dbacc(char* pstart, char*& end_acc) {
 if (pstart==NULL || *pstart==0) return NULL;
 bool hasDigits=false;
 char* pend=pstart;
 end_acc=NULL; //end of first accession
 for(char* p=pstart;;p++) {
   if (*p>='1' && *p<='9') hasDigits=true;
   if (*p==0) { //end of seq_name
      pend=p; //doesn't matter if it's accession or "db"-like
      break;
      }
   if (*p=='|') {      
      if (hasDigits) {//previous token was an accession or the last part
        pend=p; //advance pend
        if (end_acc==NULL && p-pstart>4) end_acc=p;
        }
       else //no digits -- we just got over a 'db|'
         if (pend!=pstart) break; //new database
      hasDigits=false;//reset flag
      } // pipe char
   }
 if (pend!=pstart) return pend;
              else return NULL;
}

char* parseSpToken(char* str) {
 if (str==NULL) return NULL;
 char* p=str;
 while (*p!=' ' && *p!='\t' && *p!='\v' && *p!=0) p++;
 *p=0;
 return p;
}

#define NRDB_CHARSEP "\1\2\3\4"
#define NRDB_STRSEP "^|^"

//nrdbp is positioned at the end of current nrdb concatenated
//defline, or NULL if there is no more
inline void NRDB_Rec(char* &nrdbp, char* defline) {
 nrdbp=strpbrk(defline, NRDB_CHARSEP);
 if (nrdbp==NULL) {
   nrdbp=strstr(defline, NRDB_STRSEP);
   if (nrdbp!=NULL) {
     *nrdbp='\0';
     nrdbp+=2;
     }
   }
  else *nrdbp='\0';
}

//-c/-C indexing: key up to the first space or second '|'
//receives the full defline as an unprocessed "key"
void addKeyCompact(char* defline,
              off_t fpos, uint32 reclen) {
 //we got the first token found on the defline
 num_recs++;
 char* nrdb_end;
 //breaks defline at the next nrdb concatenation point
 NRDB_Rec(nrdb_end, defline);
 //isolate the first token in this nrdb concatenation
 char* token_end=parseSpToken(defline);
 if (!compact_plus) { //shortcut key
   //only the first db|accession construct will be indexed, if found
   char* end_acc1=NULL; //end of first accession
   char* dbacc_end=parse_dbacc(defline, end_acc1);
   if (end_acc1!=NULL) { //has acceptable shortcut
     *end_acc1=0;
     add_cdbkey(defline, fpos, reclen);
     return;
     }
   if (dbacc_end!=NULL) {
      *dbacc_end=0;
      add_cdbkey(defline, fpos, reclen);
      return;
      }
   //store this whole non-space token as key:
   add_cdbkey(defline, fpos, reclen);
   return;
   }
 //from now on, only -C treatment
 for(;;) {
    //defline is on the first token    
    if (strlen(defline)>0) //add whole seq_name as the "full key"
       add_cdbkey(defline, fpos, reclen);
    //add the db|accession constructs as keys
    char* dbacc_start=defline;
    char* firstacc_end=NULL;
    char* dbacc_end=parse_dbacc(dbacc_start, firstacc_end);
    while (dbacc_end!=NULL) {
      if (firstacc_end!=NULL && firstacc_end<dbacc_end) {
        char c=*firstacc_end;
        *firstacc_end=0;
        add_cdbkey(dbacc_start, fpos, reclen);
        *firstacc_end=c;
        }
      if (dbacc_start==defline && dbacc_end==token_end) 
           break; //the whole seq_name was only one db entry
      *dbacc_end=0; //end key here
      add_cdbkey(dbacc_start, fpos, reclen);      
      if (dbacc_end==token_end) 
           break; //reached the end of this whole seq_name (non-space token)
      dbacc_start=dbacc_end+1;
      firstacc_end=NULL;
      dbacc_end=parse_dbacc(dbacc_start, firstacc_end);
      }
    // -- get to next concatenated defline, if any:
    if (compact_plus && nrdb_end!=NULL) { 
      defline=nrdb_end+1; //look for the next nrdb concatenation
      NRDB_Rec(nrdb_end, defline);
      //isolate the first token in this nrdb record
      token_end=parseSpToken(defline);
      }
     else break;
   } //for
 }

int readWords(FILE* f, GHash<int>& xhash) {
  int c;
  int count=0;
  char name[256];
  int len=0;
  while ((c=getc(f))!=EOF) {
    if (isspace(c) || c==',' || c==';') {
      if (len>0) {
        name[len]='\0';
        xhash.Add(name, new int(1));
        count++;
        len=0;
        }
      continue;
      }
    //a non-space
    name[len]=(char) c;
    len++;
    if (len>255) {
      name[len-1]='\0';
      GError("Error reading words file: token too long ('%s') !\n",name);
      }
    }
 if (len>0) {
   name[len]='\0';
   xhash.Add(name, new int(1));
   count++;
   }
 return count;
}



//========================== MAIN ===============================
int main(int argc, char **argv) {
  FILE* f_read=NULL;
  off_t fdbsize;
  int ch;
  char* zfilename;
  char* fname;
  char* marker; //record marker
  int maxkeys=0;

  int multikey;
  strcpy(record_marker,">");
  GArgs args(argc, argv, "mn:o:r:z:w:f:s:icvC");
  int e=args.isError();
  if  (e>0)
     GError("%s Invalid argument: %s\n", USAGE, argv[e] );
  if (args.getOpt('v')!=NULL) {
    printf("%s\n",VERSION);
    return 0;
    }
  multikey = (args.getOpt('m')!=NULL);
  if (multikey) {
    fields[numFields]=1;
    numFields++;
    fields[numFields]=MAX_UINT;
    numFields++;  
    }
  caseInsensitive = (args.getOpt('i')!=NULL);
  compact_plus=(args.getOpt('C')!=NULL);
  wordJunk = (char *)args.getOpt('s');
  if (wordJunk==NULL) wordJunk=defWordJunk;
  int compact=(args.getOpt('c')!=NULL || compact_plus);
  if (compact && multikey) {
    GError("%s Error: invalid flags combination.\n", USAGE);
    }
  char* s = (char*)args.getOpt('n');
  if (s!=NULL) {
    maxkeys = atoi(s);
    if (maxkeys<=1 || compact || multikey)
        GError("%s Error: invalid options (-m, -c/C, -n and -f are exclusive)\n", USAGE);
    multikey=1;
    numFields=maxkeys;
    if (numFields>254) GError(ERR_TOOMANYFIELDS);
    for (unsigned int i=1;i<=numFields;i++) fields[i-1]=i;
    }
  char* argfields = (char*)args.getOpt('f');
  if (argfields!=NULL) { //parse all the field #s
    char* pbrk;
    int prevnum=0;
    char prevsep='\0';
    numFields=0;
    char sep;
    char *p=argfields;
    do {
     pbrk=strpbrk(p,",-");
     if (pbrk==NULL) {
         sep='\0';
         pbrk=p+strlen(p);
         if (prevsep == '-' && *p=='\0' && prevnum>0) {
           //open ended range -- ending with '-'
           fields[numFields]=prevnum;
           numFields++;
           if (numFields>253) GError(ERR_TOOMANYFIELDS);
           fields[numFields]=MAX_UINT;
           numFields++;
           //GMessage("--- stored %d, %d\n",prevnum, MAX_UINT);
           break;
           }// ending with '-'
         } // '\0'
      else { sep=*pbrk; *pbrk = '\0'; }
     int num = atoi(p);
     if (num<=0 || num>254 )
              GError("%s Error: invalid syntax for -f option.\n", USAGE);
     if (prevsep == '-') { //store a range
       for (int i=prevnum;i<=num;i++) {
          fields[numFields]=i;
          numFields++;
          if (numFields>254) GError(ERR_TOOMANYFIELDS);
          }
       }
      else if (sep!='-') {
         fields[numFields]=num;
         numFields++;
         if (numFields>254) GError(ERR_TOOMANYFIELDS);
         }
      
     prevsep=sep;
     prevnum=num;
     p=pbrk+1;
     } while (sep != '\0'); //range parsing loop
    if (numFields<=0 || numFields>254 )
              GError("%s Error at parsing -f option.\n", USAGE);
    //GMessage("[%d] Fields parsed (%d values):\n", sizeof(fields[0]), numFields);
    qsort(fields, numFields, sizeof(fields[0]), &qcmpInt);    
    multikey=1;
    /*-- --------debug:
    for (unsigned int i=0;i<numFields-1;i++) {
      GMessage("%d,", fields[i]);
      }
    GMessage("%d\n",fields[numFields-1]);
    exit(0); */ 
    } //fields

  if (args.getOpt('r')!=NULL) {//non-FASTA record delimiter?
   marker=(char*)args.getOpt('r'); //
   int v=0;
   if (strlen(marker)>126) 
      GError("Error: the specified record delimiter is too long. "
        "Maximum accepted is 126\n");
   //special case: hex (0xXX) and octal codes (\XXX) are accepted, only if by themselves
   if (strlen(marker)==4 && marker[0]=='\\' || (marker[0]=='0' && (toupper(marker[1])=='X'))) {
       if (marker[0]=='\\') {
           marker++;
           v=strtol(marker, NULL, 8);
           }
          else v=strtol(marker, NULL, 16);
       if (v==0 || v>255)
         GError("Invalid record delimiter: should be only one character,\n"
                "'\\NNN' (octal value), '0xNN' (hexadecimal value)");
       record_marker[0]=v;
       record_marker_len=1;
       }
     else {
      strcpy(record_marker, marker);
      record_marker_len=strlen(record_marker);
      }
   }
 char* stopwords=(char*)args.getOpt('w'); //stop words filename given?
 if (stopwords!=NULL) {
  FILE* fstopwords=NULL;
  if ((fstopwords=fopen(stopwords, "r"))==NULL)
       GError("Cannot open stop words file '%s'!\n", stopwords);
  int c=readWords(fstopwords, stopList);
  GMessage("Loaded %d stop words.\n", c);
  fclose(fstopwords);
  useStopWords=(c>0);
  }
  if ((zfilename=(char*)args.getOpt('z')) !=NULL) {
    do_compress=true;
    strcpy(fztmp,zfilename);
    strcat(fztmp,"_ztmp");
    zf=fopen(fztmp,"wb");
    if (zf==NULL)
      GError("Error creating file '%s'\n'", fztmp);
    }
  char* outfile=(char*) args.getOpt('o');
  int numfiles = args.startNonOpt();
  if (numfiles==0)
    GError("%sError: no fasta file given.\n", USAGE);
  fname=(char*) args.nextNonOpt(); //first fasta file given
  if (do_compress)  { //-------- compression case -------------------
     if (strcmp(fname, "-")==0 || strcmp(fname, "stdin")==0)
           f_read=stdin;
       else f_read= fopen(fname, "rb");
     if (f_read == NULL) die_read(fname);
     fname=zfilename; //forget the input file name, keep the output
     }
  else {//
    int fdread= open(fname, O_RDONLY|O_BINARY);
    if (fdread == -1) die_read(fname);
    struct stat dbstat;
    fstat(fdread, &dbstat);
    fdbsize=dbstat.st_size;
    close(fdread);
    f_read= fopen(fname, "rb");
    if (f_read == NULL) die_read(fname);
    }

  char idxfile[365];
  if (outfile==NULL) {
    if (do_compress) {
      strcpy(ftmp, zfilename);
      strcat(ftmp, ".cidx");
      strcpy(idxfile, ftmp);
      strcat(ftmp, "_tmp");
      }
    else {
      strcpy(ftmp, fname);
      strcat(ftmp, ".cidx");
      strcpy(idxfile, ftmp);
      strcat(ftmp, "_tmp");
      }
    //should add the process ID, time and user to make this unique?
    }
  else {
    strcpy(ftmp, outfile);
    strcpy(idxfile, outfile);
    strcat(ftmp, "_tmp");
    }

  cdbidx=new GCdbWrite(ftmp); //test if this was successful?

  if (compact)
         addKeyFunc=&addKeyCompact;
    else if (multikey)
               addKeyFunc = &addKeyMulti;
          else addKeyFunc = &addKey;

  off_t recpos=0;
  off_t r=0;
  unsigned int recsize=0;
  char* key=NULL;
  bool fullDefline=(multikey || compact_plus);
  GReadBuf *readbuf = new GReadBuf(f_read, GREADBUF_SIZE);
  if (do_compress) { //---------------- compression case -------------
     fdbsize=0;
     GCdbz cdbz(zf); // zlib interface
     recpos=cdbz.getZRecPos();
     while ((key=cdbz.compress(readbuf,record_marker))!=NULL) {
       recsize=cdbz.getZRecSize();
       if (!fullDefline) {
         //find first space after the record_marker and place a '\0' there
         for (int i=record_marker_len; key[i]!='\0';i++) {
           if (isspace(key[i])) { key[i]='\0';break; }
           }
         }
       addKeyFunc(key, recpos, recsize);
       recpos = cdbz.getZRecPos();
       }
     remove(zfilename);
     cdbz.compress_end();
     fclose(zf);
     //determine the size of this file:
     int ftmp= open(fztmp, O_RDONLY|O_BINARY);
     if (ftmp == -1) die_read(fztmp);
     struct stat dbstat;
     fstat(ftmp, &dbstat);
     fdbsize=dbstat.st_size;
     //rename it to the intended file name
     if (rename(fztmp,zfilename) != 0) {
       GMessage("Error: unable to rename '%s' to '%s'\n",fztmp,zfilename);
       perror("rename");
       }
    }
  else { // not compressed -- buffered file access
     bool defline=false;
     int kbufsize=KBUFSIZE;
     if (fullDefline) { GMALLOC(key, KBUFSIZE); }//large defline storage buffer, just in case
                 else { GMALLOC(key, 1024); }
     int kidx=-1;
     num_recs=0;
     num_keys=0;
     char lastchar=0;
     //first iteration -- for the beginning of file case
     if (readbuf->peekCmp(record_marker, record_marker_len)==0) {
           //new record start found (defline)
           recpos=readbuf->getPos(); //new record pos
           defline=true; //we're in defline
           readbuf->skip(record_marker_len);
           kidx=0;
           }//new record start
     while ((ch=readbuf->getch())>0) {
      if (defline && kidx>=0) { //on the defline here, still parsing possible keys
          key[kidx]=(char)ch;
          kidx++;
          if (kidx>=kbufsize) {
             kbufsize+=KBUFSIZE;
             GREALLOC(key, kbufsize);
             }
          if (((isspace(ch) || ch<31) && fullDefline==false) || ch=='\n' || ch=='\r') {
               //end key here, don't care about the rest
               key[kidx-1]='\0';
               kidx=-1;
               }
        }
      if (ch=='\n') { // newline!
         //check ahead if a record delimiter follows
        if (readbuf->peekCmp(record_marker, record_marker_len)==0) {
           //new record start (defline)
           recsize = readbuf->getPos()-recpos-1; //previous recsize
           if (recsize>off_t(record_marker_len+1) && key[0]!='\0') {
             //add previous record, if there
             addKeyFunc(key, recpos, recsize);
             //GMessage("adding key=%s\n",key);
             }
           recpos=readbuf->getPos(); //new record pos
           defline=true; //we're in defline
           readbuf->skip(record_marker_len);
           //if (r<0) die_readformat(fname);
           kidx=0;
           } //new record start
        else { //after newline but not a new record start
            if (defline) { //we just finished a defline
                if (fullDefline) { //close the defline string
                  if (kidx>0) key[kidx-1]='\0';
                  kidx=-1;
                  }
               defline=false;
               }
           }
        } // was newline 
      lastchar=ch;
      }//while char
     recsize=readbuf->getPos()-recpos;
     if (recsize>0) {//add last record, if there
         if (lastchar=='\n') recsize--;
         if (fullDefline && kidx>0) {//close the defline string
               if (lastchar!='\n') kidx++;
               key[kidx-1]='\0';
               }
         addKeyFunc(key, recpos, recsize);
         //GMessage("adding key=%s\n",key);
         }
   delete readbuf;
   }
  if (f_read!=stdin) fclose(f_read);
  if (cdbidx->finish() == -1) die_write("");

  // === add some statistics at the end of the cdb index file!
  r=lseek(cdbidx->getfd(), 0, SEEK_END);
  cdbInfo info;
  memcpy((void*)info.tag, (void*)"CDBX", 4);
  info.idxflags=0;
  if (multikey) info.idxflags |= CDBMSK_OPT_MULTI;
  if (do_compress) {
      info.idxflags |= CDBMSK_OPT_COMPRESS;
      GMessage("Input data were compressed into file '%s'\n",fname);
      }
  if (compact) {
      if (compact_plus)
            info.idxflags |= CDBMSK_OPT_CADD;
         else
            info.idxflags |= CDBMSK_OPT_C;
     }
  info.num_records=gcvt_uint(&num_recs);
  info.num_keys=gcvt_uint(&num_keys);
  info.dbsize=gcvt_offt(&fdbsize);
  info.idxflags=gcvt_uint(&info.idxflags);
  int nlen=strlen(fname);
  info.dbnamelen=gcvt_uint(&nlen);
  r=write(cdbidx->getfd(), fname, nlen);
  if (r!=nlen)
        GError(ERR_W_DBSTAT);
  r=write(cdbidx->getfd(), &info, cdbInfoSIZE);
  if (r!=cdbInfoSIZE)
        GError(ERR_W_DBSTAT);
  delete cdbidx;
  GFREE(key);
  remove(idxfile);
  if (rename(ftmp,idxfile) == -1)
    GError("Error: unable to rename %s to %s",ftmp,idxfile);
  GMessage("%d entries from file %s were indexed in file %s\n",
      num_recs, fname, idxfile);
  return 0;
}
