//*******************************************************************
//* File: vfgk.cpp
//* Description: This is a Vitter implementation of the FGK (dinamic 
//* huffman coding) algorithm. We want to keep this as simple as 
//* possible for better understanding. 
//*
//* Note: check header file for version details.
//*******************************************************************

#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include "gcl/gcompress.h"

const unsigned char  ctEndChar  = 168;//0xA8;0x07;250-168
const unsigned char  ctZeroNode = 157;//0x9D;0x00;235-157
const unsigned char  ctChangeLettersMode = 166;//0xA6;0x06;246-166


Cvfgk::Cvfgk() {
  int i;
  Reset();
}

Cvfgk::~Cvfgk() {
//  if(fByteOut!=NULL)
//    fclose(fByteOut);
}

void Cvfgk::Reset() {
  char  Dict[]={'A','C','G','T','N','\n','M','S','H','D','I',
                'W','P','V','Y','F','L','K','R','Q','E','B',
                'Z','X',' ','>','0','1','2','3','4','5','6','7','8','9',
                '-','_','+','\t','@','#',',','.','"','(',')','[',']',
                ctChangeLettersMode};
  //char  Dict[]={'A','C','G','T','N','\n',ctChangeLettersMode}; 
  for(int i=0;i<2*ctAlphabetSize;i++){
    child[i]=weight[i]=parent[i]=0;
  }
  root = 1; iNextNode = 2;
  for(int i=0;i<ctAlphabetSize;i++)
     leaf[i] = 0;
  child[root]  = 0;
  weight[root] = 1;
  parent[root] = 0;
  leaf[ctZeroNode] = root;
  iStackPtr = 0;
  nBits = 0;
//  if(fByteOut!=NULL) fclose(fByteOut);
  fByteOut = NULL;
  iBytesWritten = 0;
  // start adding the characters with high probability
  for(int i=0;i<sizeof(Dict);i++){
     AddNode(Dict[i]);
     Update(Dict[i]);
  }
  // increase p for nucleotides
  for(int i=0;i<5;i++) Update(Dict[i]);
  // By default its starts with capital letters
  fCapitalLetters = true;
}

int Cvfgk::Compress(FILE* fin, FILE* fout) {
  unsigned char ch;
  int ret;
  if((fin==NULL)||(fout==NULL)) return 1;
  if(!feof(fin)) ch = ValidateInput((unsigned char)fgetc(fin));
  while(!feof(fin)) {
    if(Encode(fout,ch)) return 1;
    Update(ch);
    ch = ValidateInput((unsigned char)fgetc(fin));
  }
  if(Encode(fout,ctEndChar)) return 1;
  if(Flush(fout)) return 1;
  return 0;
}

int Cvfgk::Decompress(FILE* fin, FILE* fout) {
  unsigned char ch,chout;
  int i,n,bit,nBitLeft,ret;
  int bEndStreamNotFound = 1;
  if((fin==NULL)||(fout==NULL)) return 1;
  ret=fgetc(fin);
  if(ret==EOF)
     return 1;
  else
     ch=(unsigned char)ret;  
  nBitLeft = 8;
  while(bEndStreamNotFound) {
    n = root;
    while(child[n]!=0){
      bit = (ch>127)?1:0; ch = ch<<1;
      if(--nBitLeft==0){
        ret=fgetc(fin);
        if(ret==EOF)
           return 1;
        else
           ch=(unsigned char)ret;  
        nBitLeft = 8;
      }
      n = child[n]+bit;
    }
    chout = NodeChar[n];
    if(chout==ctZeroNode){//a new char is following
      for(i=0;i<8;i++){
         chout=(chout<<1)+((ch>127)?1:0);
         ch = ch<<1;
         if(--nBitLeft==0){
           ret=fgetc(fin);
           if(ret==EOF)
              return 1;
           else
              ch=(unsigned char)ret;  
           nBitLeft = 8;
         }
      }
      AddNode(chout);
    }
    if(chout==ctChangeLettersMode){//change the state capitals/lowers
       fCapitalLetters=!fCapitalLetters;
       Update(chout);
       continue;
    }
    if(chout==ctEndChar)// end of compressed stream
      bEndStreamNotFound = 0;
    else { 
      Update(chout);
      if(isalpha(chout))
         chout = (fCapitalLetters)?toupper(chout):tolower(chout);
      if(fputc(chout,fout)==EOF) return 1;
    }  
  }
  return 0;  
}

int Cvfgk::Compress(FILE* fin, FILE* fout, char byTermChar) {
  unsigned char ch;
  int ret;
  if((fin==NULL)||(fout==NULL)) return 1;
  if(!feof(fin)) ch = ValidateInput((unsigned char)fgetc(fin));
  while(!feof(fin)) {
    if(Encode(fout,ch)) return 1;
    Update(ch);
    ch = ValidateInput((unsigned char)fgetc(fin));
    if(ch==byTermChar){
      fseek(fin,-1,SEEK_CUR);
      break;
    }
  }
  if(Encode(fout,ctEndChar)) return 1;
  if(Flush(fout)) return 1;
  return 0;
}

int Cvfgk::CompressFasta(FILE* fin, FILE* fout, char byRecStartChar) {
  unsigned char ch;
  int ret;
  fEndOfLineReach=false;
  if((fin==NULL)||(fout==NULL)) return 1;
  if(!feof(fin)) ch = ValidateInput((unsigned char)fgetc(fin));
  while(!feof(fin)) {
    if(Encode(fout,ch)) return 1;
    Update(ch);
    ch = ValidateInput((unsigned char)fgetc(fin));
    if((ch==byRecStartChar)&&fEndOfLineReach){
      fseek(fin,-1,SEEK_CUR);
      break;
    }
    fEndOfLineReach=(ch=='\n')?true:false;
  }
  if(Encode(fout,ctEndChar)) return 1;
  if(Flush(fout)) return 1;
  return 0;
}

void Cvfgk::AddNode(char ch) {
  // transform leaf[ctZeroNode] into 
  // an internal node with weight 1
  //  a left child of weight 0 for leaf[ch]
  //  a right child of weight 1 for leaf[ctZeroNode]
  int  n = leaf[ctZeroNode];
  weight[n] = 1;
  leaf[ch] = child[n] = iNextNode;
  parent[iNextNode] = n;
  NodeChar[iNextNode] = ch;
  weight[iNextNode++] = 0;
  
  leaf[ctZeroNode] = iNextNode;
  parent[iNextNode] = n;
  NodeChar[iNextNode] = ctZeroNode;
  weight[iNextNode++] = 1;
}

void Cvfgk::Update(char ch) {
  int n;
  int m,tmp;
  unsigned char c,cn,cm;
  int bnisleaf,bmisleaf;
  int bNeedNormalize = 0;
  if(isalpha(ch)) ch=toupper(ch);
  n = leaf[ch];
  while(n!=root){
    if(weight[n]>INT_MAX-10) bNeedNormalize = 1;
    weight[n]++;
    m = n;
    while(weight[m-1]<weight[n])
       m = m-1;
    // swap nodes m,n
    bnisleaf = bmisleaf = 0;
    tmp = weight[m];
    weight[m] = weight[n];
    weight[n] = tmp;
    tmp = child[m];
    child[m] = child[n];
    child[n] = tmp;
    if(tmp) {
      parent[tmp] = parent[tmp+1] = n;
    }  
    else {
      bmisleaf = 1;cm = NodeChar[m];
    } 
    tmp = child[m];
    if(tmp) { 
      parent[tmp] = parent[tmp+1] = m;
    }  
    else {
      bnisleaf = 1;cn = NodeChar[n];
    }
    if(bmisleaf){
      leaf[cm] = n;
      NodeChar[n] = cm;
    }
    if(bnisleaf){
      leaf[cn] = m;
      NodeChar[m] = cn;
    }
    n = parent[m];   
  }
  if(weight[root]>INT_MAX-10) bNeedNormalize = 1;
  weight[root]++;
  if(bNeedNormalize) Normalize();
}

// Try to encode the next char into the internal buffer m_buff. If fail
// return an error, otherwise return 0 (NO_ERROR).
int Cvfgk::Encode(FILE* fout,unsigned char ch) {
  int  n;
  unsigned char c;
  //Check only alpha chars
  if(isalpha(ch)){
     if(islower(ch)&&fCapitalLetters){// signal the state "lowers"
        fCapitalLetters=false;
        if(Encode(fout,ctChangeLettersMode)) return 1;
        Update(ctChangeLettersMode);
     }
     if(isupper(ch)&&(!fCapitalLetters)){//signal the state "capitals"
        fCapitalLetters=true;
        if(Encode(fout,ctChangeLettersMode)) return 1;
        Update(ctChangeLettersMode);
     }
     ch=toupper(ch);
  }
  n=leaf[ch];
  iStackPtr = 0;
  if(n == 0) 
    n = leaf[ctZeroNode];
  while(n!=root){
    if(n%2) // if its odd
      stack[iStackPtr++] = 1; // is a right child
    else // if its even
      stack[iStackPtr++] = 0; // is a left child
    n = parent[n];    
  }
  // Send now stack
  if(SendStack(fout)) return 1;
  if(leaf[ch] == 0){
    // Send original char code
    c = ch;
    for(int i=0;i<8;i++){
      code = (code<<1)|((c>127)?1:0);
      c = c<<1;
      if(++nBits==8){
        if(fputc(code,fout)==EOF) return 1;
        nBits = 0;iBytesWritten++;
      }
    }
    AddNode(ch); // Split the zero-node too add the new character
  }  
  return 0;
}

int Cvfgk::Flush(FILE* fh) {
  // flush remained bits from code
  code = code << (8-nBits);
  nBits = 0;
  if(fputc(code,fh)==EOF) return 1;
  iBytesWritten++;
  return 0;
}

void Cvfgk::Normalize() {
  for(int i=0;i<2*ctAlphabetSize;i++)
     weight[i] = weight[i] >> 1;
}

int Cvfgk::SendStack(FILE* fh) {
  while(iStackPtr) {
    code = (code<<1)|stack[--iStackPtr];
    if(++nBits==8){
       if(fputc(code,fh)==EOF) return 1;
       nBits = 0;iBytesWritten++;
    }
  }
  return 0;
}

// return 1 if error
int Cvfgk::BeginByteCompression(FILE* fout){
  if(fout==NULL)
     return 1;
  Reset();    
  fByteOut = fout; 
  return 0;  
}

// return 1 if error
int Cvfgk::CompressNextByte(char ch){
  int ret;
  ch=ValidateInput((unsigned char)ch);
  if(Encode(fByteOut,ch)) return 1;
  Update(ch);
  return 0;
}

// return 0 if error or no. of bytes writte if success
int Cvfgk::EndByteCompression(){
  if(Encode(fByteOut,ctEndChar)) return 0;
  if(Flush(fByteOut)) return 0;
  fByteOut = NULL;
  return iBytesWritten;
}

unsigned char Cvfgk::ValidateInput(unsigned char c){
  if((c==ctEndChar)||(c==ctZeroNode)||(c==ctChangeLettersMode)){
    fprintf(stderr,"Input contain binary data. Reserved char encountered.\n");
    exit(1);
  }
  return c;
}


//*******************************************************************
//***********************  END OF FILE  *****************************
//*******************************************************************
