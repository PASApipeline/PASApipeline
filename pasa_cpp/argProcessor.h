#ifndef __argProcessor__
#define __argProcessor__

#include <map>
#include <string>
using namespace std;

extern map<string,bool> argSet;
extern map<string,string> argVal;

void processArgs(int argc, char* argv[]);


#endif

