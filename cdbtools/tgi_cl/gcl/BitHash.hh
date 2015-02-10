/********************************************************************************
*                  Bit operations and a hash table classes+template 
*********************************************************************************/

#ifndef BitHash_HH
#define BitHash_HH
#include "gcl/GBase.h"


/* GNumHash numeric hash template, based on a very basic hashing function,
   borrowed from bit-arrays. There are exactly n/32 buckets,
   and each bucket might have max 32 entries.
   This particular hash has a FIXED size (never grows!) allocated from the beginning
   */
class WHashNode {//structure holding actual data for a key
     public:
       int key;
       OBJ* data;
       HashNode* next;
       HashEntry(int ky=-1, OBJ* pdata=NULL, HashEntry* from=NULL) {
          key=ky;
          data=pdata;
          next=NULL;
          if (from!=NULL)
             from->next=this;
          }
       };

template <class OBJ> class  GNumHash {
 protected:
  HashNode**     hash;       //  allocated array of pointers
                              // to HashNode linked lists 
  int         BCount;     // Bucket count fCapacity/32
  int         fCapacity;  // table size - fixed!
  int         fCount;     // number of valid entries
  int  fCurrentEntry;
 public:
private:
  GFreeProc* fFreeProc; //procedure to free item data
public:
  static void DefaultFreeProc(pointer item) {
      delete (OBJ*)item;
      }
public:
  // constructs of an empty hash
  GNumHash(int cap, GFreeProc* freeProc);
  GNumHash(int cap); // allocates an empty hash with a fixed size
  void setFreeItem(GFreeProc *freeProc) { fFreeProc=freeProc; }
  void setFreeItem(bool doFree) { fFreeProc=(doFree)? &DefaultFreeProc : NULL; }
  int Capacity() const { return fCapacity; } // table's size, including the empty slots.
  int Count() const { return fCount; } // the total number of entries in the table.
  // Insert a new entry into the table given key
  // If there is already an entry with that key, leave it unchanged,
  const OBJ* Add(int ky, const OBJ* ptr);
  // Remove a given key and its data
  OBJ* Remove(int ky);
  // Find data OBJ* given key.
  OBJ* Find(int ky);
  bool hasKey(int ky) { return Find(ky); } 
  OBJ* operator[](int ky); { return Find(ky); }
  void startIterate(); //iterator-like initialization
  int NextKey(); //returns next valid key in the table (-1 if no more)
  OBJ* NextData(); //returns next valid hash[].obj or NULL if no more 
  OBJ* NextData(int& key); //returns next obj and its key
  /// Clear all entries
  void Clear();
  /// Destructor
  virtual ~GNumHash();
  };

#define S_ERR_HASHVALUE "Invalid hash value %d - capacity overflow!\n"
#define FREEDATA (fFreeProc!=NULL)
#define TEST_INDEX(x) if (x<0 || x>=fCapacity) \
     GError(S_ERR_HASHVALUE, x)

/*******************************************************************************/
// Construct empty hash

template <class OBJ> GNumHash<OBJ>::GNumHash(int cap) {
  if (cap%32!=0) cap+=32;
  fCapacity=(cap/32)*32;
  BCount=fCapacity/32;
  GCALLOC(hash, BCount*sizeof((HashEntry*)) );
  //hopefully NULL is indeed 0 on all platforms
  fCount=0;
  fFreeProc=&DefaultFreeProc;
  /*
  for (uint i=0; i<DEF_IHASH_SIZE; i++)
         hash[i].hash=-1; //this will be an indicator for 'empty' entries
  */
  }

template <class OBJ> GNumHash<OBJ>::GNumHash(int cap, GFreeProc* freeProc) {
  if (cap%32!=0) cap+=32;
  BCount=fCapacity/32;
  GCALLOC(hash, sizeof(HashEntry*)*BCount);
    //hopefully NULL is indeed 0 on all platforms
  fCount=0;
  fFreeProc=freeProc;
  }

// add a new entry, leave it alone if already existing
template <class OBJ> const OBJ* GNumHash<OBJ>::Add(int ky,
                      const OBJ* pdata){
  register int h;
  HashNode* node;
  GASSERT(fCount<fCapacity);
  h=ky/32;
  TEST_INDEX(h);
  if (hash[h]==NULL) {
     // new bucket
     hash[h] = new HashNode(ky, pdata);
     return hash[h]->pdata;
     }
  for (node=hash[h];node!=NULL;node=node->next) {
   if (node->key==ky) { //key exists, replace?
      node=new HashNode(ky, pdata, node)
      return
      }

  while(n && hash[p].hash!=-1) {
    if ((i==-1)&&(hash[p].hash==-2)) i=p;
    if (hash[p].hash==h && hash[p].key==ky) {
      //replace hash data for this key!
      hash[p].data = (void*) pdata;
      return (OBJ*)hash[p].data;
      }
    p=(p+x)%fCapacity;
    n--;
    }
  if(i==-1) i=p;
  GTRACE(("GNumHash::insert: key=\"%d\"\n",ky));
  GASSERT(0<=i && i<fCapacity);
  GASSERT(hash[i].hash<0);
  hash[i].hash=h;
  hash[i].mark=mrk;
  hash[i].key=ky;
  hash[i].data= (void*) pdata;
  fCount++;
  if((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
  GASSERT(fCount<fCapacity);
  return pdata;
  }


// Add or replace entry
template <class OBJ>  OBJ* GNumHash<OBJ>::Replace(int ky,const OBJ* pdata, bool mrk){
  register int p,i,x,h,n;
  GASSERT(fCount<fCapacity);
  //h=strhash(ky);
  h=ky;
  GASSERT(0<=h);
  p=HASH1(h,fCapacity);
  GASSERT(0<=p && p<fCapacity);
  x=HASH2(h,fCapacity);
  GASSERT(1<=x && x<fCapacity);
  i=-1;
  n=fCapacity;
  while(n && hash[p].hash!=-1){
    if((i==-1)&&(hash[p].hash==-2)) i=p;
    if(hash[p].hash==h && hash[p].key==ky) {
      if(hash[p].mark<=mrk){
        GTRACE(("GNumHash::replace: %08x: replacing: \"%d\"\n",this,ky));
        if (FREEDATA) (*fFreeProc)(hash[p].data);
        hash[p].mark=mrk;
        hash[p].data=pdata;
        }
      return hash[p].data;
      }
    p=(p+x)%fCapacity;
    n--;
    }
  if(i==-1) i=p;
  GTRACE(("GNumHash::replace: %08x: inserting: \"%s\"\n",this,ky));
  GASSERT(0<=i && i<fCapacity);
  GASSERT(hash[i].hash<0);
  hash[i].hash=h;
  hash[i].mark=mrk;
  hash[i].key=ky;
  hash[i].data=pdata;
  fCount++;
  if ((100*fCount)>=(MAX_LOAD*fCapacity)) Resize(fCount);
  GASSERT(fCount<fCapacity);
  return pdata;
  }


// Remove entry
template <class OBJ> OBJ* GNumHash<OBJ>::Remove(int ky){
  register int p,x,h,n;
  //if(!ky){ GError("GNumHash::remove: NULL key argument.\n"); }
  if(0<fCount){
    //h=strhash(ky);
    h=ky;
    GASSERT(0<=h);
    p=HASH1(h,fCapacity);
    GASSERT(0<=p && p<fCapacity);
    x=HASH2(h,fCapacity);
    GASSERT(1<=x && x<fCapacity);
    GASSERT(fCount<fCapacity);
    n=fCapacity;
    while(n && hash[p].hash!=-1){
      if(hash[p].hash==h && hash[p].key==ky){
        GTRACE(("GNumHash::remove: %08x removing: \"%d\"\n",this,ky));
        hash[p].hash=-2;
        hash[p].mark=false;
        if (FREEDATA) (*fFreeProc)(hash[p].data);
        hash[p].key=-1;
        hash[p].data=NULL;
        fCount--;
        if((100*fCount)<=(MIN_LOAD*fCapacity)) Resize(fCount);
        GASSERT(fCount<fCapacity);
        return NULL;
        }
      p=(p+x)%fCapacity;
      n--;
      }
    }
  return NULL;
  }


// Find entry
template <class OBJ> bool GNumHash<OBJ>::hasKey(int ky) {
  register int p,x,h,n;
  //if(!ky){ GError("GNumHash::find: NULL key argument.\n"); }
  if(0<fCount){
    //h=strhash(ky);
    h=ky;
    GASSERT(0<=h);
    p=HASH1(h,fCapacity);
    GASSERT(0<=p && p<fCapacity);
    x=HASH2(h,fCapacity);
    GASSERT(1<=x && x<fCapacity);
    GASSERT(fCount<fCapacity);
    n=fCapacity;
    while(n && hash[p].hash!=-1){
      if(hash[p].hash==h && hash[p].key==ky){
        return true;
        }
      p=(p+x)%fCapacity;
      n--;
      }
    }
  return false;
}

template <class OBJ> OBJ* GNumHash<OBJ>::Find(int ky){
  register int p,x,h,n;
  //if(!ky){ GError("GNumHash::find: NULL key argument.\n"); }
  if(0<fCount){
    //h=strhash(ky);
    h=ky;
    GASSERT(0<=h);
    p=HASH1(h,fCapacity);
    GASSERT(0<=p && p<fCapacity);
    x=HASH2(h,fCapacity);
    GASSERT(1<=x && x<fCapacity);
    GASSERT(fCount<fCapacity);
    n=fCapacity;
    while(n && hash[p].hash!=-1) {
      if(hash[p].hash==h && hash[p].key==ky){
        return (OBJ*)hash[p].data;
        }
      p=(p+x)%fCapacity;
      n--;
      }
    }
  return NULL;
  }


template <class OBJ> void GNumHash<OBJ>::startIterate() {// initialize a key iterator; call
 fCurrentEntry=0;
}

template <class OBJ> int GNumHash<OBJ>::NextKey() {
 register int pos=fCurrentEntry;
 while (pos<fCapacity && hash[pos].hash<0) pos++;
 if (pos==fCapacity) {
                 fCurrentEntry=fCapacity;
                 return NULL;
                 }
              else {
                 fCurrentEntry=pos+1;
                 return hash[pos].key;
                 }
}

template <class OBJ> const OBJ* GNumHash<OBJ>::NextData() {
 register int pos=fCurrentEntry;
 while (pos<fCapacity && hash[pos].hash<0) pos++;
 if (pos==fCapacity) {
                 fCurrentEntry=fCapacity;
                 return NULL;
                 }
              else {
                 fCurrentEntry=pos+1;
                 return (OBJ*)hash[pos].data;
                 }

}

template <class OBJ> const GIHashEntry* GNumHash<OBJ>::NextEntry() {
 register int pos=fCurrentEntry;
 while (pos<fCapacity && hash[pos].hash<0) pos++;
 if (pos==fCapacity) {
                 fCurrentEntry=fCapacity;
                 return NULL;
                 }
              else {
                 fCurrentEntry=pos+1;
                 return &hash[pos];
                 }
}


// Get first non-empty entry
template <class OBJ> int GNumHash<OBJ>::First() const {
  register int pos=0;
  while(pos<fCapacity){ if(0<=hash[pos].hash) break; pos++; }
  GASSERT(fCapacity<=pos || 0<=hash[pos].hash);
  return pos;
  }

// Get last non-empty entry
template <class OBJ> int GNumHash<OBJ>::Last() const {
  register int pos=fCapacity-1;
  while(0<=pos){ if(0<=hash[pos].hash) break; pos--; }
  GASSERT(pos<0 || 0<=hash[pos].hash);
  return pos;
  }


// Find next valid entry
template <class OBJ> int GNumHash<OBJ>::Next(int pos) const {
  GASSERT(0<=pos && pos<fCapacity);
  while(++pos <= fCapacity-1){ if(0<=hash[pos].hash) break; }
  GASSERT(fCapacity<=pos || 0<=hash[pos].hash);
  return pos;
  }


// Find previous valid entry
template <class OBJ> int GNumHash<OBJ>::Prev(int pos) const {
  GASSERT(0<=pos && pos<fCapacity);
  while(--pos >= 0){ if(0<=hash[pos].hash) break; }
  GASSERT(pos<0 || 0<=hash[pos].hash);
  return pos;
  }


// Remove all
template <class OBJ> void GNumHash<OBJ>::Clear(){
  register int i;
  for(i=0; i<fCapacity; i++){
    if(hash[i].hash>=0){
      if (hash[i].keyalloc) GFREE((hash[i].key));
      if (FREEDATA)
         (*fFreeProc)(hash[i].data);
      }
    }
  GFREE(hash);
  GMALLOC(hash, sizeof(GIHashEntry)*DEF_IHASH_SIZE);
  //reinitialize it
  for (i=0; i<DEF_IHASH_SIZE; i++)
         hash[i].hash=-1; //this will be an indicator for 'empty' entries
  fCapacity=DEF_IHASH_SIZE;
  fCount=0;
  }


// Destroy table
template <class OBJ> GNumHash<OBJ>::~GNumHash(){
  register int i;
  for(i=0; i<fCapacity; i++){
    if(hash[i].hash>=0){
      if (FREEDATA) (*fFreeProc)(hash[i].data);
      }
    }
  GFREE(hash);
  }

#endif
