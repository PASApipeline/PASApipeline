#ifndef __argProcessor__
#define __argProcessor__

#include <map>
#include <string>
using namespace std;

extern map<string,bool> argSet;
extern map<string,string> argVal;

bool processArgs(int argc, char* argv[]);


#endif

