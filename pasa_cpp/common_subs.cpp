#include "common_subs.h"
#include <iostream>
using namespace std;

void swap_ints(int& i, int& j) {
  int temp = i;
  i = j;
  j = temp;
}


int min (int a, int b) {
  if (a < b) {
    return (a);
  } else {
    return (b);
  }
}


int min (vector<int>& myList) {
  if (myList.size() == 0) {
    cerr << "Can't determine minimum of an empty vector." << endl;
    exit(1);
  }
  int Min = myList[0];
  for (int i = 1; i < myList.size(); i++) {
    if (myList[i] < Min) {
      Min = myList[i];
    }
  }
  return (Min);
}

int max (int a, int b) {
  if (a > b) {
    return (a);
  } else {
    return (b);
  }
}




int max (vector<int>& myList) {
  if (myList.size() == 0) {
    cerr << "Can't determine maximum of an empty vector." << endl;
    exit(1);
  }
  int Max = myList[0];
  for (int i = 1; i < myList.size(); i++) {
    if (myList[i] > Max) {
      Max = myList[i];
    }
  }
  return (Max);
}


bool overlap (struct coordset& a, struct coordset& b) {
  if (a.lend <= b.rend && a.rend >= b.lend) { //overlap
    return (true);
  } else {
    return (false);
  }
}


int** twoDarray (int x, int y) {
  
  // make array x rows and y columns
  int** int2D = new int*[x];
  for (int i = 0; i < y; i++) {
    int2D[i] = new int[y];
    // init entries to zero
    for (int j = 0; j < y; j++) {
      int2D[i][j] = 0;
    }
  }

  return (int2D);
}
 
bool** twoDarray (int x, int y, bool init) {
  
  bool** twoD = new bool*[x];
  for (int i = 0; i < y; i++) {
    twoD[i] = new bool[y];
    // init values to init
    for (int j=0; j < y; j++) {
      twoD[i][j] = init;
    }
  }
  
  return (twoD);

}
   

void free2Darray (int** twoD) {
  
  /*  temporarily disabling while I recode the 2D-array
      
      // first, free the full space alloc:
      int* firstPtr = twoD[0];
      delete [] firstPtr;
      
      // now, delete the ** array
      delete [] twoD;

  */

}

void free2Darray (bool** twoD) {
  
  /*  temporarily disablling while I recode the 2D-array
      // first, free the full space alloc:
      bool* firstPtr = twoD[0];
      delete [] firstPtr;
      
      // now, delete the ** array
      delete [] twoD;
      
  */
  
}



