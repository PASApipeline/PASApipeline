#include "argProcessor.h"

map<string,bool> argSet;
map<string,string> argVal;

bool processArgs (int argc, char* argv[]) {
  for (int i=1; i < argc; i++) {
    char* arg = argv[i];
    string myArg (arg);
    if (arg[0] == '-') { 
      
      if (i == argc-1 || argv[i+1][0] == '-') { //last value or followed by another -opt
	// argument specified with no value
	argSet[myArg] = true;
      }
      
      // see if next argument is another arg or value:
      if (i != argc-1) {
	// value:
	string nextArg (argv[i+1]);
	argVal[myArg] = nextArg;
      }
    }
  }
}

