//*******************************************************************
//* File: vfgk.h
//* Description: This is a Vitter implementation of the FGK (dinamic 
//* huffman coding) algorithm.  
//*
//* Notes: Algorithm V has been proven to get better results than original
//* fgk and is comparable with what bzip2 can achieve for our data. However,
//* future development of this class
//* will include usage of the skeleton trees which will increase the 
//* computational speed more than 50%. Hopefully the class will became 
//* part of the tgi class library and will support serialization, RTT and
//* c++ exception style handling. But I guess you already know that ;) 
//*
//* Modifications:
//* 
//*  v1.0 Jul-20-2002  - first version has been born
//*  v2.0 Oct-11-2002  - initial dictionary added
//*  v3.0 Oct-28-2002  - capital letters made an subset of the others
//*                      halfing in this way the alphabet.
//*                    - init bug fixed
//*******************************************************************

#ifndef VFGK_HEADER_INCLUDED
#define VFGK_HEADER_INCLUDED

#include <stdio.h>

const int   ctAlphabetSize =  256;
//Simple debugprint macro
#ifdef _DEBUG
#define DBGPRINT(exp) do{ \
    printf("At line %d in file %s DEBUG says:\n" __LINE__,__FILE__); \
    printf(exp); \
    }while(0)
#else
#define DBGPRINT(exp) ((void) 0)
#endif
// Simple assert macro
#ifdef _DEBUG
#define ASSERT(exp) ((exp)||(printf("Assert Failure at Line %d in file %s",__LINE__,__FILE__))) 
#else
#define ASSERT(exp) ((void) 0) 
#endif

/* / just a base class
class CGenericClass {
public:
   CGenericClass(){};
   ~CGenericClass(){};
// we don't use exception in this implementation
#ifdef _DEBUG
   virtual void		Error(char* err){printf("got the error:%s",err);};
#else
   virtual void		Error(char* err) {exit(100);};
#endif  

};
*/
// Compressor
 class Cvfgk {
  public:
    //construction/destruction
    Cvfgk();
    ~Cvfgk();
    // Methods
    void	Reset();
    int			Compress(FILE* fin, FILE* fout);
    int			Compress(FILE* fin, FILE* fout, char);
    int			CompressFasta(FILE* fin, FILE* fout, char);
    int			Decompress(FILE* fin, FILE* fout);
    // For byte-to-byte compression
    int			BeginByteCompression(FILE* fout);// return 0 if success
    int			CompressNextByte(char);// return 0 if success
    int			EndByteCompression();// return No. of bytes written or 0 if fail
  private:
    // main data structures to hold our tree
    int		root;
    int		iNextNode;
    // child[n] is the left child of the node n
    // child[n]+1 is the right child of node n
    int		child[2*ctAlphabetSize]; 
    int		weight[2*ctAlphabetSize];
    int		parent[2*ctAlphabetSize]; 
    // leaf[ch] is the leaf node for char ch
    int		leaf[ctAlphabetSize];
    // NodeChar[n] is the char code of leaf n
    char	NodeChar[2*ctAlphabetSize];
    // a primitive stack
    char	stack[ctAlphabetSize+1];
    int		iStackPtr;
    int		nBits; // bit counter for code
    char	code;
    // used as a state only in fasta file compressions
    bool	fEndOfLineReach;
    // for byte-to-byte compression
    FILE*	fByteOut;
    int		iBytesWritten;
    // state of letters in output
    bool	fCapitalLetters;
  private:
    // Helpers
    int Encode(FILE*,unsigned char);
    void Update(char);
    int Flush(FILE*);
    void AddNode(char);
    void Normalize();
    int SendStack(FILE*);
    unsigned char ValidateInput(unsigned char);
 };

#endif 

//*******************************************************************
//***********************  END OF FILE  *****************************
//*******************************************************************





