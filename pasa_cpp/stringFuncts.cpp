#include "stringFuncts.h"
#include <stdlib.h>

int stringToInt (string s) {
  return (atoi(s.c_str()));
}

float stringToFloat (string s) {
  return (atof(s.c_str()));
}

vector<string> stringSplitter (string s, string delimeter) {
  int begin = 0;
  int length = s.length();
  int delimeter_length = delimeter.length();
  
  int pos = s.find(delimeter.c_str());
  vector<string> wordList;
  while (pos != -1) {
    int sublength = pos - begin;
    string word = s.substr(begin, sublength);
    wordList.push_back(word);
    begin = pos + delimeter_length;
    pos = s.find(delimeter.c_str(), begin);
  }
  if (begin != 0) {
    string word = s.substr(begin);
    wordList.push_back(word);
  } else if (begin == 0) {
    // no delimeter located.
    wordList.push_back(s);
  }

  return (wordList);
}


