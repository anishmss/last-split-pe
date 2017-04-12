// Copyright 2017 Naruki Yoshikawa
// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

#include "last-pair-probs.hh"

#include "io.hh"
#include "stringify.hh"

#include <algorithm>
#include <cctype>  // isalpha
#include <cerrno>
#include <cmath>
#include <cstdlib>  // atof
#include <cstring>  // strncmp
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <limits.h>
#include <cfloat>
#include <stddef.h>  // size_t

typedef const char *String;

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

struct Alignment {
  double score;
  std::string qName, rName;
  long rStart, qStart, size, rSrcSize, qSrcSize;
  char rStrand, qStrand;
  std::string qSeq, rSeq;
  std::vector<double> prob;
};

class AlignmentParameters {
  // Parses the score scale factor, minimum score, and genome size.
  double t;  // score scale factor
  double e;  // minimum score
  long   g;  // genome size

 public:
  AlignmentParameters() : t(-1), e(-1), g(-1) {}  // dummy values

  void update(const std::string& line) {
    std::istringstream ss(line);
    std::string i;
    while (ss >> i) {
      const char *c = i.c_str();
      if (t == -1.0 && i.substr(0,2) == "t=") {
	cbrc::unstringify(t, c + 2);
        if (t <= 0) err("t must be positive");
      }
      if (e == -1.0 && i.substr(0,2) == "e=") {
	cbrc::unstringify(e, c + 2);
        if (e <= 0) err("e must be positive");
      }
      if (g == -1.0 && i.substr(0,8) == "letters=") {
	cbrc::unstringify(g, c + 8);
        if (g <= 0) err("letters must be positive");
      }
    }
  }

  bool isValid() const { return t != -1 && e != -1 && g != -1; }

  void validate() const {
    if (t == -1) err("I need a header line with t=");
    if (e == -1) err("I need a header line with e=");
    if (g == -1) err("I need a header line with letters=");
  }

  double tGet() const { return t; }
  double eGet() const { return e; }
  long   gGet() const { return g; }
};

static AlignmentParameters readHeaderOrDie(std::istream& lines) {
  std::string line;
  AlignmentParameters params;
  while (getline(lines, line)) {
    if (line.substr(0,1) == "#") {
      params.update(line);
      if (params.isValid())
        return params;
    } else if (line.find_first_not_of(" ") != std::string::npos) {
      break;
    }
  }
  params.validate();  // die
  return params;  // dummy
}

static std::vector<Alignment> readAlignmentsSet(std::istream& input) {
    std::vector<Alignment> A;
    std::string line;
    while(true) {
        int n = 0;
        double score;
        std::string str_score, rName, qName, rSeq, qSeq, prob;
        char rStrand, qStrand;
        long rStart, qStart, rSrcSize, qSrcSize, size;
        char type;
        bool initial = true;
        while (std::getline(input, line)) {
            std::string head = line.substr(0,1);
            std::stringstream ss(line);
            if(initial) {
                if(head == "#") continue;
                else initial = false;
            }
            if(head == "#") {
                break;
            } else if(head == "a") {
                std::string str_score;
                ss >> type >> str_score;
                std::stringstream ss2(str_score);
                ss2.ignore(6);
                ss2 >> score;
            } else if(head == "s") {
                if(n==0) {
                    ss >> type >> rName >> rStart >> size >> rStrand >> rSrcSize >> rSeq;
                    n = 1;
                } else if(n==1) {
                    ss >> type >> qName >> qStart >> size >> qStrand >> qSrcSize >> qSeq;
                    n = 2;
                } else {
                    err("bad Maf format!");
                }
            } else if(head == "p") {
                ss >> type >> prob;
                // std::cout << score << std::endl << rName << std::endl << qName << std::endl;
                Alignment aln;
                aln.size = size;
                aln.score = score;
                aln.qName = qName;
                aln.rName = rName;
                aln.qStart = qStart;
                aln.rStart = rStart;
                aln.rSrcSize = rSrcSize;
                aln.qSrcSize = qSrcSize;
                aln.rStrand = rStrand;
                aln.qStrand = qStrand;
                aln.qSeq = qSeq;
                aln.rSeq = rSeq;
                aln.prob.resize(prob.size());
                for(size_t i=0; i<aln.prob.size(); ++i) {
                    aln.prob[i] = prob[i];
                }
                A.push_back(aln);
                n = 0;
            }
        }
        return A;
    }
}


std::vector<std::vector<double>> UpdateProbability(std::vector<Alignment>& X, std::vector<Alignment>& Y, AlignmentParameters& params) {
    const double prob_disjoint = 0.1;
    size_t sizeX = X[0].qSrcSize;
    std::vector<std::vector<double>> p_Hij(sizeX, std::vector<double>(X.size(), 0.0));
    for(size_t i=0; i<sizeX; ++i) {
        for(size_t j=0; j<X.size(); ++j) {
            if(i >= X[j].qStart && i < X[j].qStart+X[j].size) {
                char c = X[j].prob[sizeX - X[j].qStart];
                p_Hij[i][j] = pow(10.0, -((c - 33) / 10.0));
            }
        }
    }
    std::vector<std::vector<double>> p_y_Hij(sizeX, std::vector<double>(X.size()));
    for(size_t i=0; i<sizeX; ++i) {
        for(size_t j=0; j<X.size(); ++j) {
            for(size_t k=0; k<Y.size(); ++k) {
                double inferredf = 0.1;     // TODO: infer f correctly
                p_y_Hij[i][j] += inferredf * exp(Y[k].score/params.tGet()) * (1-prob_disjoint);
               p_y_Hij[i][j] +=  exp(Y[k].score/params.tGet()) * prob_disjoint / 2.0 / params.gGet();
            }
        }
    }
    std::vector<std::vector<double>> p_Hij_y(sizeX, std::vector<double>(X.size()));
    double Z = 0;
    for(size_t i=0; i<sizeX; ++i) {
        for(size_t l=0; l<Y.size(); ++l) {
            Z += p_y_Hij[i][l] * p_Hij[i][l];
        }
        for(size_t j=0; j<Y.size(); ++j) {
            p_Hij_y[i][j] = p_y_Hij[i][j] * p_Hij[i][j] / Z;
        }
    }
    
    for(size_t j=0; j<X.size(); ++j) {
        for(size_t i=0; i<sizeX; ++i) {
            char c = static_cast<char>(std::max(0.0, -log10(1.0-p_y_Hij[i][j])*10.0) + 33);
            std::cout << c;
        }
        std::cout << std::endl;
    }
    return p_Hij_y;
}

void lastSplitPe(LastPairProbsOptions& opts) {
    const std::vector<std::string>& inputs = opts.inputFileNames;
    size_t n = inputs.size();
    std::ifstream inFile1;
    std::istream& input = (n > 0) ? cbrc::openIn(inputs[0], inFile1) : std::cin;
    AlignmentParameters params = readHeaderOrDie(input);
    std::vector<Alignment> X = readAlignmentsSet(input);
    std::vector<Alignment> Y = readAlignmentsSet(input);
    auto result = UpdateProbability(X, Y, params);
}
