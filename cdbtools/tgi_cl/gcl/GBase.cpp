#include <stdarg.h>
#include <ctype.h>
#include "gcl/GBase.h"
#ifdef __WIN32__
  #include <windows.h>
#endif

#define IUPAC_DEFS "AaCcTtGgUuMmRrWwSsYyKkVvHhDdBbNnXx-*"
#define IUPAC_COMP "TtGgAaCcAaKkYyWwSsRrMmBbDdHhVvNnXx-*"

unsigned char ntCompTable[256];

void ntCompTableInit() {
       if (ntCompTableReady) return;
       char n[]=IUPAC_DEFS;
       char c[]=IUPAC_COMP;
       int l=strlen(IUPAC_DEFS);
       ntCompTable[0]=0;
       for (int ch=1;ch<256;ch++) {
          ntCompTable[ch]=0;
          for (int i=0;i<l;i++)
                if (ch==n[i]) {
                  ntCompTable[ch]=c[i];
                  break;
                  }
          if (ntCompTable[ch]==0)
              ntCompTable[ch]='N';
          }
      ntCompTableReady=true;
     }

char ntComplement(char c) {
 return ntCompTable[(int)c];
 }

bool ntCompTableReady=false;

static char msg[4069];
//************************* Debug helpers **************************
// Assert failed routine
void GAssert(const char* expression, const char* filename, unsigned int lineno){
  sprintf(msg,"%s(%d): GASSERT(%s) failed.\n",filename,lineno,expression);
  fprintf(stderr,"%s",msg);
  }
// Error routine (prints error message and exits!)
void GError(const char* format,...){
  #ifdef __WIN32__	
    va_list arguments;
    va_start(arguments,format);
    vsprintf(msg,format,arguments);
    va_end(arguments);
    OutputDebugString(msg);
    fprintf(stderr,"%s",msg); // if a console is available
    MessageBox(NULL,msg,NULL,MB_OK|MB_ICONEXCLAMATION|MB_APPLMODAL);
  #else
    va_list arguments;
    va_start(arguments,format);
    vfprintf(stderr,format,arguments);
    va_end(arguments);
    //abort();
  #endif
    exit(1);
  }
// Warning routine (just print message without exiting)
void GMessage(const char* format,...){
  va_list arguments;
  va_start(arguments,format);
  vsprintf(msg,format,arguments);
  va_end(arguments);
  #ifdef __WIN32__
    OutputDebugString(msg);
  #endif
  fprintf(stderr,"%s",msg);fflush(stderr);
  }

/*************** Memory management routines *****************/
// Allocate memory
bool GMalloc(pointer* ptr,unsigned long size){
  //GASSERT(ptr);
  if (size!=0) *ptr=malloc(size);
  return *ptr!=NULL;
  }
  
// Allocate cleaned memory (0 filled)
bool GCalloc(pointer* ptr,unsigned long size){
  GASSERT(ptr);
  *ptr=calloc(size,1);
  return *ptr!=NULL;
  }
// Resize memory
bool GRealloc(pointer* ptr,unsigned long size){
  //GASSERT(ptr);
  if (size==0) {
    GFree(ptr);
    return true;
    }
  if (*ptr==NULL) {//simple malloc
   void *p=malloc(size);
   if (p != NULL) {
     *ptr=p;
     return true;
     }
    else return false;
   }//malloc
  else {//realloc  
   void *p=realloc(*ptr,size);
   if (p) {
       *ptr=p;
       return true;
       }
   return false;
   }
 }
// Free memory, resets ptr to NULL afterward
void GFree(pointer* ptr){
  GASSERT(ptr);
  if (*ptr) free(*ptr);
  *ptr=NULL;
  }

char* Gstrdup(const char* str) {
  char *copy;
  GMALLOC(copy, strlen(str)+1);
  strcpy(copy,str);
  return copy;
  }

char* Gsubstr(const char* str, char* from, char* to) {
 //extract (and allocate) a substring, including boundaries (from/to)
 int len=strlen(str);
 if (to==NULL) to=(char *)(str+len-1); //extract tail from 'from' char, including it
 if (!(from>=str && to>from && to<str+len)) return NULL;
 int newlen=to-from+1;
 char* subs;
 GMALLOC(subs, newlen);
 memcpy(subs, str, newlen-1);
 subs[newlen]='\0';
 return subs;
 }

char* replaceStr(char* &str, char* newvalue) {
 if (str!=NULL) GFREE(str);
 if (newvalue==NULL) { return NULL; }
 GMALLOC(str, strlen(newvalue)+1);
 strcpy(str,newvalue);
 return str;
 }

void* Gmemscan(void *mem, unsigned int len,
                   void *part, unsigned int partlen) {
char* p;
unsigned int restlen=len-partlen+1;
void* oldp=mem;
while ( (p=(char*)memchr(oldp, ((char*)part)[0], restlen))!=NULL) {
  //located first char, try to match the rest:
  p++;
  if (memcmp(p, &((char*)part)[1], partlen-1)==0) return p-1;
  //no string match, prepare next iteration
  restlen-=(p-(char*)oldp);
  oldp=p;
  }//while
return NULL;
}

//rindex function is missing on some platforms ?
char* rstrchr(char* str, char ch) {  /* returns a pointer to the rightmost
  occurence of ch in str  */
 char *p;
 if (str==NULL) return NULL;
 p=str+strlen(str)-1;
 while (p>=str) {
    if (*p==ch) return p;
    p--;
    }
 return NULL;
 }


/* DOS/UNIX safer fgets : reads a text line from a (binary) file and
  update the file position accordingly and the buffer capacity accordingly.
  The given buf is resized to read the entire line in memory 
    -- even when it's abnormally long
  */
char* fgetline(char *buf, int& buf_cap, FILE *stream, long& f_pos) {
  //reads a char at a time until \n and/or \r are encountered
  int i=0;
  int c=0;
  while ((c=getc(stream))!=EOF) {
    if (i>=buf_cap-1) {
       buf_cap+=1024;
       GREALLOC(buf,buf_cap);
    }
    if (c=='\n' || c=='\r') {
      buf[i]='\0';
      if (c=='\r') {
        if ((c=getc(stream))!='\n') ungetc(c,stream);
                               else f_pos++;
        }
      f_pos++;
      return buf;
      }
    f_pos++;
    buf[i]=(char)c;
    i++;
    } //while i<buf_cap-1
  if (c==EOF && i==0) return NULL;
  buf[i]='\0';
  return buf;
  }

char* GLineBuf::getLine(FILE* stream, off_t& f_pos) {
   //reads a char at a time until \n and/or \r are encountered
   len=0;
   int c=0;
   while ((c=getc(stream))!=EOF) {
     if (len>=allocated-1) {
        allocated+=1024;
        GREALLOC(buf, allocated);
     }
     if (c=='\n' || c=='\r') {
       buf[len]='\0';
       if (c=='\r') { //DOS file -- special case
         if ((c=getc(stream))!='\n') ungetc(c,stream);
                                else f_pos++;
         }
       f_pos++;
       return buf;
       }
     f_pos++;
     buf[len]=(char)c;
     len++;
     } //while i<buf_cap-1
   if (c==EOF) {
     isEOF=true;
     if (len==0) return NULL;
     }
   buf[len]='\0';
   return buf;
}


//strchr but with a set of chars instead of only one
char* strchrs(char* s, const char* chrs) {
  if (s==NULL || chrs==NULL || *chrs=='\0' || *s=='\0')
         return NULL;
  unsigned int l=strlen(s);
  unsigned int r=strcspn(s, chrs);
  if (r==l) return NULL;
  return (s+r);
}

char* upCase(const char* str) {
 if (str==NULL) return NULL;
 int len=strlen(str);
 char* upstr;
 GMALLOC(upstr, len+1);
 upstr[len]='\0';
 for (int i=0;i<len;i++) upstr[i]=toupper(str[i]);
 return upstr;
 }

char* loCase(const char* str) {
 if (str==NULL) return NULL;
 int len=strlen(str);
 char* lostr;
 GMALLOC(lostr, len+1);
 lostr[len]='\0';
 for (int i=0;i<len;i++) lostr[i]=tolower(str[i]);
 return lostr;
 }

//test if a char is in a given string (set)
bool chrInStr(char c, char* str) {
 if (str==NULL || *str=='\0') return false;
 for (char* p=str; (*p)!='\0'; p++) {
   if ((*p)==c) return true;
   }
 return false;
 }



char* rstrfind(char* str, char* substr) {
/* like rindex() for a string */
 int l,i;
 if (str==NULL || *str=='\0') return NULL;
 if (substr==NULL || *substr=='\0') return NULL;
 char* p=str+strlen(str)-strlen(substr); 
   //rightmost position that could match
 l=strlen(substr);  
 while (p>=str) {
    for (i=0; i<l && *(p+i) == *(substr+i); i++);
    if (i==l) return p; //found!
    p--;
    }
 return NULL;
}

// tests if string s has the given prefix
bool startsWith(char* s, const char* prefix) {
 if (prefix==NULL || s==NULL) return false;
 int i=0;
 while (prefix[i]!='\0' && prefix[i]==s[i]) i++;
 return (prefix[i]=='\0');
 }


char* reverseChars(char* str, int slen) {
  if (slen==0) slen=strlen(str);
  int l=0;
  int r=slen-1;
  register char c;
  while (l<r) {
     c=str[l];str[l]=str[r];
     str[r]=c;
     //swap(str[l],str[r]);
     l++;r--;
     }
  return str;
}


char* rstrstr(char* rstart, char *lend, char* substr) {  /*like strstr, but starts searching
 from right end, going up to lend and returns a pointer to the last (right) 
 matching character in str */
 char *p;
 int l,i;
 l=strlen(substr);
 p=rstart-l+1;
 while (p>=lend) {
    for (i=0;i<l;i++) if (*(p+i) != *(substr+i)) break;
    if (i==l) return p+l-1;
    p--;
    }
 return NULL;
 }


//hash function used for strings in GHash
int strhash(const char* str){
  register int h=0;
  register int g;
  while (*str) {
    h=(h<<4)+*str++;
    g=h&0xF0000000;
    if(g) h^=g>>24;
    h&=0x0fffffff;
    }
  GASSERT(h<=0x0fffffff);
  return h;
  }

// removes the directory part from a full-path file name
// this is a destructive operation for the given string!!!
// the trailing '/' is guaranteed to be there
void delFileName(char* filepath) {
 char *p, *sep;
 if (filepath==NULL) return;
 for (p=filepath, sep=filepath;*p!='\0';p++)
     if (*p==CHPATHSEP) sep=p+1;
 *sep='\0'; // truncate filepath
}

// returns a pointer to the file name part in a full-path filename
char* getFileName(char* filepath) {
 char *p, *sep;
 if (filepath==NULL) return NULL;
 for (p=filepath, sep=filepath;*p!='\0';p++)
     if (*p==CHPATHSEP) sep=p+1;
 return sep;
}


bool fileExists(char* filepath) {
  if (filepath==NULL) return false;
  FILE* ft=fopen(filepath, "rb");
  if (ft==NULL) return false;
  fclose(ft);
  return true;
}


bool parseNumber(char* &p, double& v) {
 //skip any spaces..
 while (*p==' ' || *p=='\t') p++;
 char* start=p;
 if (*p=='-') p++;
       else if (*p=='+') { p++;start++; }
 while ((*p>='1' && *p<='9') || *p=='0' || *p=='.') p++;
 //now p is on a non-digit;
 if (*start=='-' && p==start+1) return false;
 char saved=*p;
 *p='\0';
 char* endptr=p;
 v=strtod(start,&endptr);
 *p=saved;
 if (endptr!=p) return false;
 return true;
}

bool parseInt(char* &p, int& i) {
 while (*p==' ' || *p=='\t') p++;
 char* start=p;
 if (*p=='-') p++;
       else if (*p=='+') { p++;start++; }
 while ((*p>='1' && *p<='9') || *p=='0') p++;
 //now p is on a non-digit;
 if (*start=='-' && p==start+1) return false;
 char saved=*p;
 *p='\0';
 char* endptr=p;
 long l=strtol(start,&endptr,10);
 i=(int)l;
 *p=saved;
 if (endptr!=p || i!=l) return false;
 return true;
}

