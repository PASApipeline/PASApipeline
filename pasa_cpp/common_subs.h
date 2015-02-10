#ifndef __COMMON_SUBS__
#define __COMMON_SUBS__

#include "common_structs.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void swap_ints(int& i, int& j);

int min (int, int);
int min (vector<int>); // returns min of integer list

int max (int, int);
int max (vector<int>); // returns max of integer list

bool overlap (struct coordset& a, struct coordset& b); // deterimines if coordsets overlap
  
int** twoDarray (int x, int y);
bool** twoDarray (int x, int y, bool init);

void free2Darray (int** twoD);
void free2Darray (bool** twoD);

#endif

