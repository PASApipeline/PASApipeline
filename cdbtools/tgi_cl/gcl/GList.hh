//---------------------------------------------------------------------------
/*
Sortable collection of pointers to objects
*/

#ifndef GListHH
#define GListHH
#include "gcl/GBase.h"


#define SLISTINDEX_ERR "GList error:Invalid list index: %d"
#define SLISTCAPACITY_ERR "GList error: invalid capacity: %d"
#define SLISTCOUNT_ERR "GList error: invalid count: %d"
#define SLISTSORTED_ERR "Cannot do that for a sorted list!"
#define SLISTUNSORTED_ERR "Cannot do that for an unsorted list!"
// ------ macros:
#define BE_UNSORTED if (fCompareProc!=NULL) { GError(SLISTSORTED_ERR); return; }
#define BE_SORTED if (fCompareProc==NULL) { GError(SLISTUNSORTED_ERR); return; }

#define MAXLISTSIZE INT_MAX-1

#define SORTED (fCompareProc!=NULL)
#define UNSORTED (fCompareProc==NULL)
#define FREEDATA (fFreeProc!=NULL)
#define TEST_INDEX(x) if (x<0 || x>=fCount) \
     GError(SLISTINDEX_ERR, x)

template <class OBJ> class GList {
  protected:
    OBJ** fList; //pointer to an array of pointers to objects
    int fCount; //total number of entries in list
    int fCapacity; //currGent allocated size
    bool fUnique;
    GCompareProc* fCompareProc; //a pointer to a Compare function
    GFreeProc* fFreeProc; //useful for deleting objects
    static int DefaultCompareProc(const pointer item1, const pointer item2) {
      //the comparison operators MUST be defined for OBJ class!
      if (*((OBJ*)item1) > *((OBJ*)item2)) return 1;
        else if (*((OBJ*)item2) > *((OBJ*)item1)) return -1;
                                             else return  0;
      }
  protected:
    void sortInsert(int idx, OBJ* item);
  public:
    static void DefaultFreeProc(pointer item) {
      delete (OBJ*)item;
      }
    GList(GCompareProc* compareProc=NULL); //free by default
    GList(GCompareProc* compareProc, //unsorted by default
        GFreeProc *freeProc,
        bool beUnique=false);
    GList(bool sorted, bool free_elements=true, bool beUnique=false);
    GList(int init_capacity, bool sorted, bool free_elements=true, bool beUnique=false);
    GList(GList<OBJ>& list); //copy constructor..
    virtual ~GList();
    void freeItem(int idx);
    void setSorted(GCompareProc* compareProc);
       //sorted if compareProc not NULL; sort the list if compareProc changes !
    void setFreeItem(GFreeProc *freeProc) { fFreeProc=freeProc; }
    void setFreeItem(bool doFree) {
       if (doFree) fFreeProc=DefaultFreeProc;
             else  fFreeProc=NULL;
       }
    bool Sorted() { return fCompareProc!=NULL; }
    void setSorted(bool sorted) {
     if (sorted) {
         if (fCompareProc!=&DefaultCompareProc) {
             fCompareProc=&DefaultCompareProc;
             Sort();
             }
          }
      else fCompareProc=NULL;
      }
    int Add(OBJ* item); //-- specific implementation if sorted
    void Clear();
    void Delete(int index);
    void Forget(int idx);
    void Exchange(int idx1, int idx2);
    void Expand();
    OBJ* First() { return Get(0); }
    OBJ* Last()  { return (fCount>0)?fList[fCount-1]:NULL; }
    int Capacity() { return fCapacity; }
    int Unique() { return fUnique; }
    //this will reject identical items in sorted lists only!
    void setUnique(bool beUnique) { fUnique = beUnique; };

    void setCapacity(int NewCapacity);
    int Count() { return fCount; }
    void setCount(int NewCount);
    GCompareProc* GetCompareProc() {return fCompareProc;}
    OBJ* Get(int idx);
    OBJ* operator[](int i) {
          TEST_INDEX(i);
          return fList[i];
          }
    void Grow();
    void Grow(int idx, OBJ* item);
    int IndexOf(OBJ* item); //this has a specific implementation for sorted lists
               //if list is sorted, item data is located by binary search
               //based on the Compare function
               //if not, a linear search is performed, but
               //this needs the == operator to be defined for OBJ class! 
    bool Found(OBJ* item, int & idx); // sorted only;
               //search by content; if found, returns true and idx will be the index
               //of the first item found matching for which GTCompareProc returns 0
    void Insert(int idx, OBJ* item); //unsorted only, place item at position idx
    void Move(int curidx, int newidx);
    void Put(int idx, OBJ* item);
    int Remove(OBJ* item);
    int RemovePtr(OBJ* item); //use always linear search to find the pointer!
    void Pack();
    void Sort(); //explicit sort may be requested using this function
    void QuickSort(int L, int R);
    const GList<OBJ>& operator=(GList& list);
};


//-------------------- TEMPLATE IMPLEMENTATION-------------------------------

template <class OBJ> GList<OBJ>::GList(GList& list) { //copy constructor
 fCount=0;
 fCapacity=0;
 fList=NULL;
 fCompareProc=list.fCompareProc;
 fFreeProc=list.fFreeProc;
 for (int i=0;i<list.Count();i++) Add(list[i]);
}

template <class OBJ> GList<OBJ>::GList(GCompareProc* compareProc, 
       GFreeProc* freeProc, bool beUnique) {
  fCount=0;
  fCapacity=0;
  fList=NULL;
  fCompareProc = compareProc;
  fFreeProc    = freeProc;
  fUnique = beUnique; //only affects sorted lists
}

template <class OBJ> GList<OBJ>::GList(GCompareProc* compareProc) {
  fCount=0;
  fCapacity=0;
  fList=NULL;
  fCompareProc = compareProc;
  fFreeProc    = &DefaultFreeProc;
  fUnique = false; //only affects sorted lists
}

template <class OBJ> GList<OBJ>::GList(bool sorted, 
    bool free_elements, bool beUnique) {
  fCount=0;
  fCapacity=0;
  fList=NULL;
  if (sorted) {
     if (free_elements) {
        fCompareProc=&DefaultCompareProc;
        fFreeProc=&DefaultFreeProc;
        fUnique=beUnique;
        }
       else {
        fCompareProc=&DefaultCompareProc;
        fFreeProc=NULL;
        fUnique=beUnique;
        } 
     }
   else {
     if (free_elements) {
        fCompareProc=NULL;
        fFreeProc=&DefaultFreeProc;
        fUnique=beUnique;       
        }
      else {
        fCompareProc=NULL;
        fFreeProc=NULL;
        fUnique=beUnique;
        }
     }
}

template <class OBJ> GList<OBJ>::GList(int init_capacity, bool sorted,
    bool free_elements, bool beUnique) {
  fCount=0;
  fCapacity=0;
  fList=NULL;
  if (sorted) {
     if (free_elements) {
        fCompareProc=&DefaultCompareProc;
        fFreeProc=&DefaultFreeProc;
        fUnique=beUnique;
        }
       else {
        fCompareProc=&DefaultCompareProc;
        fFreeProc=NULL;
        fUnique=beUnique;
        }
     }
   else {
     if (free_elements) {
        fCompareProc=NULL;
        fFreeProc=&DefaultFreeProc;
        fUnique=beUnique;
        }
      else {
        fCompareProc=NULL;
        fFreeProc=NULL;
        fUnique=beUnique;
        }
     }
 setCapacity(init_capacity);
}



template <class OBJ> GList<OBJ>::~GList() {
 Clear();//this will free the items if fFreeProc is defined
}

template <class OBJ> void GList<OBJ>::setCapacity(int NewCapacity) {
  if (NewCapacity < fCount || NewCapacity > MAXLISTSIZE)
    GError(SLISTCAPACITY_ERR, NewCapacity);
    //error: capacity not within range
  if (NewCapacity!=fCapacity) {
   if (NewCapacity==0)
      GFREE(fList);
    else  
      GREALLOC(fList, NewCapacity*sizeof(OBJ*));
   fCapacity=NewCapacity;
   }
}

template <class OBJ> void GList<OBJ>::freeItem(int idx) {
  TEST_INDEX(idx);
  (*fFreeProc)(fList[idx]);
  fList[idx]=NULL;
}

template <class OBJ> void GList<OBJ>::Clear() {
 if (FREEDATA) {
   for (int i=0; i<fCount; i++) {
     (*fFreeProc)(fList[i]);
     }
   }
 GCompareProc* fcmp=fCompareProc;
 fCompareProc=NULL;
 setCount(0);
 setCapacity(0); //so the array itself is deallocated too!
 fCompareProc=fcmp;
}


template <class OBJ> void GList<OBJ>::Exchange(int idx1, int idx2) {
 BE_UNSORTED; //cannot do that in a sorted list!
 TEST_INDEX(idx1);
 TEST_INDEX(idx2);
 OBJ* item=fList[idx1];
 fList[idx1]=fList[idx2];
 fList[idx2]=item;
}


template <class OBJ> void GList<OBJ>::Expand()
{
 if (fCount==fCapacity) Grow();
 //return this;
}


template <class OBJ> OBJ* GList<OBJ>::Get(int idx)
{
 TEST_INDEX(idx);
 return fList[idx];
}

template <class OBJ> const GList<OBJ>& GList<OBJ>::operator=(GList& list) {
 if (&list!=this) {
     Clear();
     fCompareProc=list.fCompareProc;
     fFreeProc=list.fFreeProc;
     //attention: the pointers are copied directly,
     //but their content is not!!! which is rather stupid
     for (int i=0;i<list.Count();i++) Add(list[i]);
     }
 return *this;
}

template <class OBJ> void GList<OBJ>::setSorted(GCompareProc* compareProc) {
 GCompareProc* old_proc=fCompareProc;
 fCompareProc=compareProc;
 if (fCompareProc!=old_proc && fCompareProc!=NULL) 
       Sort(); //new compare method
}

template <class OBJ> void GList<OBJ>::Grow() {
 int delta;
 if (fCapacity > 64) delta = fCapacity/4;
   else if (fCapacity > 8) delta = 16;
                      else delta = 4;
  setCapacity(fCapacity + delta);
}

template <class OBJ> void GList<OBJ>::Grow(int idx, OBJ* newitem) {
 int delta;
 if (fCapacity > 64) delta = fCapacity/4;
   else if (fCapacity > 8) delta = 16;
                      else delta = 4;
 // setCapacity(fCapacity + delta);
 int NewCapacity=fCapacity+delta;
  if (NewCapacity <= fCount || NewCapacity > MAXLISTSIZE)
    GError(SLISTCAPACITY_ERR, NewCapacity);
    //error: capacity not within range
  if (NewCapacity!=fCapacity) {
    if (NewCapacity==0)
      GFREE(fList);
    else  {//add the new item
      if (idx==fCount) {
         GREALLOC(fList, NewCapacity*sizeof(OBJ*));
         fList[idx]=newitem;
         }
       else {
        OBJ** newList;
        GMALLOC(newList, NewCapacity*sizeof(OBJ*));
        //copy data before idx
        memmove(&newList[0],&fList[0], idx*sizeof(OBJ*));
        newList[idx]=newitem;
        //copy data after idx
        memmove(&newList[idx+1],&fList[idx], (fCount-idx)*sizeof(OBJ*));
        memset(&newList[fCount+1], (long)NULL, NewCapacity-fCount-1);
        //data copied:
        GFREE(fList);
        fList=newList;
        }
      fCount++;  
      }
   fCapacity=NewCapacity;
   } 
}


template <class OBJ> int GList<OBJ>::IndexOf(OBJ* item) {
 int result=0;
 if (Found(item, result)) return result;
                     else return -1;
 }

template <class OBJ> int GList<OBJ>::Add(OBJ* item) {
 int result;
 if (SORTED) {
   if (Found(item, result))
      if (fUnique) return -1; //cannot add a duplicate!
   //Found sets result to the position where the item should be!
   sortInsert(result, item);
   }
  else {
   if (fUnique && Found(item,result)) return -1; //set behaviour
   result = fCount;
   if (result==fCapacity) Grow();
   fList[result]=item;
   fCount++;
   }
 return result;
}


template <class OBJ> bool GList<OBJ>::Found(OBJ* item, int& idx) {
 //search the list by using CompareProc (if defined)
 //or == operator for a non-sortable list
 //for sorted lists, even when the result is false, the idx is
 //set to the closest matching object!
 int i;
 idx=-1;
 if (fCount==0) { idx=0;return false;}
 if (SORTED) { //binary search based on CompareProc
   //do the simple test first:
   
   if ((*fCompareProc)(fList[0],item)>0) {
                       idx=0;
                       return false;
                       }
   if ((*fCompareProc)(item, fList[fCount-1])>0) {
                       idx=fCount;
                       return false;
                       }

   int l, h, c;
   l = 0;
   h = fCount - 1;
   while (l <= h) {
       i = (l + h) >> 1;
       c = (*fCompareProc)(fList[i], item);
       if (c < 0)  l = i + 1;
         else {
            h = i - 1;
            if (c == 0) {
                 idx=i;
                 return true;
                }
            }
       } //while
   idx = l;
   return false;
   }
 else {//not sorted: use linear search
   // needs == operator to compare user defined objects !
   i=0;
   while (i<fCount) {
      if (*fList[i]==*item) {
         idx=i;
         return true;
         }
      i++;
      }
   return false;
   }
}

template <class OBJ> void GList<OBJ>::sortInsert(int idx, OBJ* item) {
 //idx must be the new position this new item must have
 //so the allowed range is [0..fCount]
 //the old idx item all the above will be shifted to idx+1
 if (idx<0 || idx>fCount) GError(SLISTINDEX_ERR, idx);
 if (fCount==fCapacity) {
    Grow(idx, item); 
    //expand and also copy/move data and insert the new item
    return;
    }
 //room still left, just move data around and insert the new one
 if (idx<fCount) //copy/move pointers only!
      memmove(&fList[idx+1], &fList[idx], (fCount-idx)*sizeof(OBJ*));
 fList[idx]=item;
 fCount++;
}

template <class OBJ> void GList<OBJ>::Insert(int idx, OBJ* item) {
 //idx can be [0..fCount] so an item can be actually added
 BE_UNSORTED; //cannot do that with a sorted list!
 if (idx<0 || idx>fCount) GError(SLISTINDEX_ERR, idx);
 if (fCount==fCapacity) {
   Grow(idx, item);
   return;
   }
 if (idx<fCount)
      memmove(&fList[idx+1], &fList[idx], (fCount-idx)*sizeof(OBJ*));
 fList[idx]=item;
 fCount++;
}

template <class OBJ> void GList<OBJ>::Move(int curidx, int newidx) {
 BE_UNSORTED; //cannot do that in a sorted list!
 if (curidx!=newidx || newidx>=fCount)
     GError(SLISTINDEX_ERR, newidx);
 OBJ* p;
 p=Get(curidx);
 //this is a delete:
 fCount--;
 if (curidx<fCount)
    memmove(&fList[curidx], &fList[curidx+1], (fCount-curidx)*sizeof(OBJ*));
 //-this was instead of delete
 Insert(newidx, p);
}


template <class OBJ> void GList<OBJ>::Put(int idx, OBJ* item) {
 //WARNING: this will never free the replaced item!!!
 TEST_INDEX(idx);
 fList[idx]=item; 
 if (SORTED) Sort(); //re-sort 
  //! if NULL was given here, CompareProc will crash!
}

template <class OBJ> void GList<OBJ>::Forget(int idx) {
 TEST_INDEX(idx);
 fList[idx]=NULL;
}

template <class OBJ> void GList<OBJ>::Delete(int index) {
 TEST_INDEX(index);
 if (fFreeProc!=NULL) {
   (*fFreeProc)(fList[index]); //freeItem
   }
 fList[index]=NULL;
 fCount--;
 if (index<fCount)
   memmove(&fList[index], &fList[index+1], (fCount-index)*sizeof(OBJ*));
}

template <class OBJ> int GList<OBJ>::Remove(OBJ* item) {
//removes an item if it's in our list
 int result=IndexOf(item);
 if (result>=0) Delete(result);
 return result;
}

//linear search for the pointer and mark it as invalid
template <class OBJ> int GList<OBJ>::RemovePtr(OBJ* item) {
int i;
if (item==NULL) return -1;
for (i=0;i<fCount;i++)
   if (fList[i]==item) break;
if (i==fCount) return -1; //not found
Delete(i);
return i;
}


template <class OBJ> void GList<OBJ>::Pack()  {//also frees items!
 for (int i=fCount-1; i>=0; i--) if (fList[i]==NULL) Delete(i);
}

template <class OBJ> void GList<OBJ>::setCount(int NewCount) {
  if (NewCount<0 || NewCount > MAXLISTSIZE)
     GError(SLISTCOUNT_ERR, NewCount);
  if (NewCount > fCapacity) setCapacity(NewCount);
  if (NewCount > fCount)
    memset(fList[fCount], 0, (NewCount - fCount) * sizeof(OBJ*));
  fCount = NewCount;
}

template <class OBJ> void GList<OBJ>::QuickSort(int L, int R) {
 int I, J;
 OBJ* P;
 OBJ* T;
 do {
    I = L;
    J = R;
    P = fList[(L + R) >> 1];
    do {
      while (fCompareProc(fList[I], P) < 0) I++;
      while (fCompareProc(fList[J], P) > 0) J--;
      if (I <= J) {
        T = fList[I];
        fList[I] = fList[J];
        fList[J] = T;
        I++;
        J--;
        }
      }
    while (I <= J);
    if (L < J) QuickSort(L, J);
    L = I;
    }
 while (I < R);
 
}

template <class OBJ> void GList<OBJ>::Sort() {
 if (fList!=NULL && fCount>0 && fCompareProc!=NULL)
     QuickSort(0, fCount-1);
}

//---------------------------------------------------------------------------
#endif
