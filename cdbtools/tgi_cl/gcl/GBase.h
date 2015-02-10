#ifndef G_BASE_DEFINED
#define G_BASE_DEFINED

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>

typedef unsigned int uint32;
typedef int int32;
typedef unsigned char uchar;


#if defined(_NATIVE_64) || defined(_LP64) || defined(__LP64__)
 typedef long int64;
 typedef unsigned long uint64;
#else
 //assume 32bit environment with long long for int64 stuff
 typedef long long int64;
 typedef unsigned long long uint64;
#endif

/****************************************************************************/

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

/****************************************************************************/
#define ERR_ALLOC "Error allocating memory.\n"
#ifdef __WIN32__
  #define CHPATHSEP '\\'
 #else
  #define CHPATHSEP '/'
#endif

//-------------------




// Debug helpers
#ifndef NDEBUG
 #define GASSERT(exp) ((exp)?((void)0):(void)GAssert(#exp,__FILE__,__LINE__))
 #ifdef TRACE
  #define GTRACE(exp)  (GMessage exp)
 #else
  #define GTRACE(exp)  ((void)0)
 #endif 
#else
 #define GASSERT(exp) ((void)0)
 #define GTRACE(exp)  ((void)0)
#endif

#define GERROR(exp) (GError exp)
/**********************************  Macros  ***********************************/
// Abolute value
#define GABS(val) (((val)>=0)?(val):-(val))

// Min and Max
#define GMAX(a,b) (((a)>(b))?(a):(b))
#define GMIN(a,b) (((a)>(b))?(b):(a))

// Min of three
#define GMIN3(x,y,z) ((x)<(y)?GMIN(x,z):GMIN(y,z))

// Max of three
#define GMAX3(x,y,z) ((x)>(y)?GMAX(x,z):GMAX(y,z))

// Return minimum and maximum of a, b
#define GMINMAX(lo,hi,a,b) ((a)<(b)?((lo)=(a),(hi)=(b)):((lo)=(b),(hi)=(a)))

// Clamp value x to range [lo..hi]
#define GCLAMP(lo,x,hi) ((x)<(lo)?(lo):((x)>(hi)?(hi):(x)))

typedef void* pointer;
typedef unsigned int uint;

typedef int GCompareProc(const pointer item1, const pointer item2);
typedef void GFreeProc(pointer item); //usually just delete,
      //but may also support structures with embedded dynamic members

#define GMALLOC(ptr,size)  if (!GMalloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GCALLOC(ptr,size)  if (!GCalloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GREALLOC(ptr,size) if (!GRealloc((pointer*)(&ptr),size)) \
                                     GError(ERR_ALLOC)
#define GFREE(ptr)       GFree((pointer*)(&ptr))

inline char* min(char *arg1, char *arg2) {
    return (strcmp(arg1, arg2) < 0)? arg1 : arg2;
}

inline int iround(double x) {
   return (int)floor(x + 0.5);
}


/****************************************************************************/

inline char* max(char *arg1, char *arg2) {
    return (strcmp(arg2, arg1) < 0)? arg1 : arg2;
}

inline int Gintcmp(int a, int b) {
 return (a>b)? 1 : ((a==b)?0:-1);
}

inline void swap(int &arg1, int &arg2){
 arg1 ^= arg2 ^= arg1 ^= arg2;
 }

inline void swap(char* &arg1, char* &arg2){
 register char* swp=arg1;
 arg1=arg2; arg2=swp; 
 }
 
inline void swap(unsigned int &arg1, unsigned int &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(short &arg1, short &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(unsigned short &arg1, unsigned short &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(long &arg1, long &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(unsigned long &arg1, unsigned long &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(char &arg1, char &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(unsigned char &arg1, unsigned char &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }

inline void swap(bool &arg1, bool &arg2)
{ arg1 ^= arg2 ^= arg1 ^= arg2; }


/**************** Memory management ***************************/

bool GMalloc(pointer* ptr, unsigned long size); // Allocate memory
bool GCalloc(pointer* ptr, unsigned long size); // Allocate and initialize memory
bool GRealloc(pointer* ptr,unsigned long size); // Resize memory
void GFree(pointer* ptr); // Free memory, resets ptr to NULL

/********************* debug functions *********************/

void GError(const char* format,...); // Error routine (aborts program)
void GMessage(const char* format,...);// Log message to stderr
// Assert failed routine:- usually not called directly but through GASSERT
void GAssert(const char* expression, const char* filename, unsigned int lineno);


// ****************** string manipulation *************************
char *Gstrdup(const char* str);
//duplicate a string by allocating a copy for it and returning it
char* Gsubstr(const char* str, char* from, char* to=NULL);
//extracts a substring, allocating it, including boundaries (from/to)

char* replaceStr(char* &str, char* newvalue);

//conversion: to Lower/Upper case
char* upCase(const char* str);
char* loCase(const char* str);

//strstr but for memory zones: scans a memory region
//for a substring:
void* Gmemscan(void *mem, unsigned int len,
                   void *part, unsigned int partlen);

// test if a char is in a string:
bool chrInStr(char c, char* str);
 
char* rstrchr(char* str, char ch);
/* returns a pointer to the rightmost
  occurence of ch in str - like rindex for platforms missing it*/

char* strchrs(char* s, const char* chrs);
//strchr but with a set of chars instead of only one

char* rstrfind(char* str, char *substr); /* like rindex() but for strings
or like the right side version of strstr()
*/
//reverse character string or 
char* reverseChars(char* str, int slen=0);

char* rstrstr(char* rstart, char *lend, char* substr);
/*the reversed, rightside equivalent of strstr: starts searching
 from right end (rstart), going back to left end (lend) and returns 
 a pointer to the last (right) matching character in str */


//Determines if a string begins with a given prefix
//(returns false when any of the params is NULL, 
// but true when prefix is '' (empty string)!) 
bool startsWith(char* s, const char* prefix);

// ELF hash function for strings
int strhash(const char* str);

//--------------------------------------------------------
// ************** simple line buffer class for text files

//GLineBuf -- text line reading/buffering class
class GLineBuf {
   int len;
   int allocated;
   char* buf;
   bool isEOF;
   FILE* file;
   off_t filepos; //current position
 public:
   char* chars() { return buf; }
   char* line() { return buf; }
   int length() { return len; }
   int size() { return len; } //same as size();
   bool isEof() {return isEOF; }
   bool eof() { return isEOF; }
   off_t getfpos() { return filepos; }
   off_t getFpos() { return filepos; }
   char* getLine() { return getLine(file);  }
   char* getLine(FILE* stream) { return getLine(stream, filepos); }
   char* getLine(FILE* stream, off_t& f_pos); //read a line from a stream and update
                                              // the given file position
   GLineBuf(FILE* stream=NULL, off_t fpos=0) {
     len=0;
     isEOF=false;
     allocated=1024;
     GMALLOC(buf,allocated);
     file=stream;
     filepos=fpos;
     }
   ~GLineBuf() {
     GFREE(buf);
     }
};


/* DOS/UNIX safer (?) fgets -  to read one line from a file and
  update the file position correctly ! */
char* fgetline(char* buf, int& n, FILE* stream, long& f_pos);

/*********************** File management functions *********************/

// removes the directory part from a full-path file name
// this is a destructive operation for the given string!
void delFileName(char* filepath);

// returns a pointer to the file name part in a full-path filename
char* getFileName(char* filepath);

bool fileExists(char* filepath);

//parses the next number found in a string at the current position
//until a non-digit (and not a '.') is encountered;
//updates the char* pointer to be after the last digit parsed
bool parseNumber(char* &p, double& v);
bool parseInt(char* &p, int& i);

extern bool ntCompTableReady;
char ntComplement(char c);
void ntCompTableInit();


#endif /* G_BASE_DEFINED */ 
