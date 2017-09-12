// Copyright 2017 Naruki Yohshikawa
// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith


#ifndef LAST_SPLIT_PE_HH
#define LAST_SPLIT_PE_HH

#include <string>
#include <vector>
#include <set>

struct LastPairProbsOptions {
  bool rna;
  bool estdist;
  double mismap;
  bool isFraglen;
  double fraglen;
  bool isSdev;
  double sdev;
  bool isDisjoint;
  double disjoint;
  std::set<std::string> circular;
  std::vector<std::string> inputFileNames;
  double outer;
  double inner;
  double disjointScore;
  double maxMissingScore1;
  double maxMissingScore2;
  bool isSamFormat;  
};

void lastSplitPe(LastPairProbsOptions& opts);

#endif
