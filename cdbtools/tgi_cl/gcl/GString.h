//---------------------------------------------------------------------------
#ifndef GStringH
#define GStringH
//---------------------------------------------------------------------------
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "gcl/GBase.h"

// This class uses reference counting and copy-on-write semantics 

// All indexes are zero-based.  For all functions that accept an index, a
// negative index specifies an index from the right of the string.  Also,
// for all functions that accept a length, a length of -1 specifies the rest
// of the string.
enum enTokenizeMode {
 tkFullString,
 tkCharSet
 };

class GString {
        friend GString operator+(const char* s1, const GString& s2);
        friend bool operator==(const char* s1, const GString& s2);
        friend bool operator<(const char* s1, const GString& s2);
        friend bool operator<=(const char* s1, const GString& s2);
        friend bool operator>(const char* s1, const GString& s2);
        friend bool operator>=(const char* s1, const GString& s2);
        friend bool operator!=(const char* s1, const GString& s2);
        friend void swap(GString& s1, GString& s2);
    public:
        GString();
        GString(const GString& s);
        GString(const char* s);
        GString(const int i);
        GString(const double f);        
        GString(char c, int n = 1);
        ~GString();
        operator const char* () const { return my_data->chars;} //inline here
        char& operator[](int index);
        char operator[](int index) const;
        GString& operator=(const GString& s);
        GString& operator=(const char* s);
        GString& operator=(const int i);
        GString& operator=(const double f);
        GString operator+(const GString& s) const;
        GString operator+(const char* s) const;
        GString operator+(const char c) const;
        GString operator+(const int i) const;
        GString operator+(const double f) const;
        bool operator==(const GString& s) const;
        bool operator==(const char* s) const;
        bool operator<(const GString& s) const;
        bool operator<(const char* s) const;
        bool operator<=(const GString& s) const;
        bool operator<=(const char* s) const;
        bool operator>(const GString& s) const;
        bool operator>(const char* s) const;
        bool operator>=(const GString& s) const;
        bool operator>=(const char* s) const;
        bool operator!=(const GString& s) const;
        bool operator!=(const char* s) const;
        GString& operator+=(const GString& s);
        GString& operator+=(const char* s);
        GString& operator+=(const char c);
        GString& operator+=(const int i);
        GString& operator+=(const double f);        
      //interface:
      public:
        int length() const;
        bool is_empty() const;
        bool is_space() const;
        GString substr(int index = 0, int len = -1) const;
        GString to(char c); //return the first part up to first occurence of c
                           //or whole string if c not found
        GString from(char c); //same as to, but starting from the right side
        GString copy() const;
        GString& format(const char *fmt,...);
        GString& appendfmt(const char *fmt,...);
        GString& cut(int index = 0, int len = -1); //delete a specified length
        //paste a string at the specified position
        GString& paste(const GString& s, int index = 0, int len=-1);
        GString& paste(const char* s, int index = 0, int len = -1);
        GString& replace(const char* from, const char* to=NULL);
        GString& insert(const GString& s, int index = 0);
        GString& insert(const char* s, int index = 0);
        GString& append(const char* s);
        GString& append(const GString& s);
        GString& upper();
        GString& lower();
        GString& clear();//make empty
        //character translation or removal:
        GString& tr(char* from, char* to=NULL);
        //number of occurences of a char in the string:
        int count(char c);
        void startTokenize(const char* delimiter, enTokenizeMode tokenizemode=tkCharSet);
        bool nextToken(GString& token);
        int asInt(int base=10);
        double asReal();
        int index(const GString& s, int start_index = 0) const;
        int index(const char* s, int start_index = 0) const;
        int index(char c, int start_index = 0) const;
        int rindex(char c) const;
        int rindex(char* str) const;
        bool contains(const GString& s) const;
        bool contains(const char* s) const;
        bool contains(char c) const;
        GString split(char* delim);
        GString split(char c);
           /* splits "this" in two parts, at the first (leftmost) 
                 encounter of delim:
                 1st would stay in "this" 
                 (which this way is truncated)
                 2nd will go to the returned string
           */
        GString splitr(char* delim);
        GString splitr(char c);
           /* splits "this" in two parts, at the last (rightmost) 
                 encounter of delim:
                 1st would stay in "this"
                 2nd will be returned
           */

        int peelInt() const; //extract an integer, (left to right), from a
                //mixed alphanumeric string, e.g. 'T24HC1234b'=> 2
        int peelIntR() const; //same as above, but starts from the right side
        //e.g. 'T2HC1234b'=> 1234
        GString& trim(char c);
        GString& trim(char* c=" \t\n\r");
        GString& trimR(char* c=" \t\n\r"); //trim only right end
        GString& trimR(char c=' ');
        GString& trimL(char* c=" \t\n\r"); //trim only left end
        GString& trimL(char c=' ');
        GString& padR(int len, char c=' '); //align it in len spaces to the right
        GString& padL(int len, char c=' '); //align it in len spaces to the left
        GString& padC(int len, char c=' '); //center it
        size_t read(FILE* stream, char* delimiter="\n", size_t bufsize=4096);
          //read next token from stream, using the given string as
          //a marker where the block should stop

        static const int max_token_size = 200;
        static const int max_line_size = 600;
        const char* chars() const;
        const char* text() const;
    protected:
        char* fTokenDelimiter;
        int fLastTokenStart;
        enTokenizeMode fTokenizeMode;
        void* readbuf; //file read buffer for the read() function
        size_t readbufsize; //last setting for the readbuf 
        static void invalid_args_error(const char* fname);
        static void invalid_index_error(const char* fname);
        struct Data {//structure holding actual
                     //string data and reference count information
               Data() { ref_count=0; length=0; chars[0] = '\0'; }
               unsigned int ref_count;
               int length;
               char chars[1];
              };
        static Data* new_data(int length); //alloc a specified length string's Data
        static Data* new_data(const char* str); //alloc a copy of a specified string
        void replace_data(int length);
        void replace_data(Data* data);
        void make_unique();
        char* chrs(); // this is dangerous, length should not be affected
        static Data null_data; //a null (empty) string Data is available here
        Data* my_data; //pointer to a Data object holding actual string data
};

/***************************************************************************/

inline int GString::length() const {
 return my_data->length;
 }


inline const char *GString::chars() const {
 return my_data->chars;
 }

inline char *GString::chrs() { //protected version, allows modification of the chars
 return my_data->chars;
 }

inline const char *GString::text() const {
 return my_data->chars;
 }


inline bool operator>=(const char *s1, const GString& s2) {
 return (strcmp(s1, s2.chars()) >= 0);
 }

inline bool operator!=(const char *s1, const GString& s2) {
 return (strcmp(s1, s2.chars()) != 0);
 }

inline void swap(GString& s1, GString& s2) {
 GString::Data *tmp = s1.my_data; s1.my_data = s2.my_data;
 s2.my_data = tmp;
 }


#endif
