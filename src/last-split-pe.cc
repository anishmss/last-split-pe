// Copyright 2017 Naruki Yoshikawa
// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

#include "last-pair-probs.hh"

#include "io.hh"
#include "stringify.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static double logSumExp(const double a, const double b) {
  // Adds numbers, in log space, to avoid overflow.
  const double m = std::max(a, b);
  return std::log(std::exp(a-m)+std::exp(b-m)) + m;
}

struct Alignment {
  double score;
  std::string qName, rName;
  long rStart, qStart, size, rSrcSize, qSrcSize;
  char rStrand, qStrand;
  std::string qSeq, rSeq, prob;
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

class AlignmentPair {
public:
    std::string refName;
    long refIndex;
    char refStrand;
    double probability;
    char qBase;
    int best_pair;
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

static std::vector<Alignment> readAlignmentSet(std::istream& input) {
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
                aln.prob = prob;
                A.push_back(aln);
                n = 0;
            }
        }
        return A;
    }
}


void calcProbAndOutput(std::vector<Alignment>& X, std::vector<Alignment>& Y, AlignmentParameters& params, LastPairProbsOptions& opts) {
    // p(I=1) i.e. probability of disjoint
    const double prob_disjoint = 0.01;

    size_t sizeX = X[0].qSrcSize;
    std::vector<long> i_j(X.size(), -1);
    std::vector<double> log_p_Hj(X.size(), 0.0);
    std::vector<double> log_p_y_Hj(X.size(), 0.0);
    std::vector<double> p_Hj_y(X.size(), 0.0);
    
    std::vector<AlignmentPair> alnpair;
    for(size_t i=0; i<sizeX; ++i) {
        std::fill(i_j.begin(), i_j.end(), -1);
        std::fill(log_p_Hj.begin(), log_p_Hj.end(), -1.0e99);
        std::fill(log_p_y_Hj.begin(), log_p_y_Hj.end(), -1.0e99);
        std::fill(p_Hj_y.begin(), p_Hj_y.end(), 0.0);
        double p_R = 0.0;
        for(size_t j=0; j<X.size(); ++j) {
            // prior probability: p(H_j)
            if(i >= X[j].qStart && i < X[j].qStart+X[j].size) {
                i_j[j] = X[j].rStart + i - X[j].qStart;
                char c = X[j].prob[i - X[j].qStart];
                double p = 1.0 - std::pow(10.0, -((c - 33) / 10.0));
                log_p_Hj[j] = std::log(p);
                p_R += p;
            }
            // conditional probability p(Y_k | H_j)
            for(size_t k=0; k<Y.size(); ++k) {
                // estimate flagment length (TODO: check again!)
                int flag_len = -1;
                if(X[j].qStrand=='+' && Y[k].qStrand=='-') {
                    flag_len = (i - X[j].qStart - 1) + (Y[k].rStart - i_j[j] + 1) + Y[k].qSrcSize;
                } else if(X[j].qStrand=='+' && Y[k].qStrand=='-') {
                    flag_len = i_j[j] - Y[k].rStart + i_j[j] + X[j].qSrcSize - (i - X[j].qStart - 1);
                }
                // p(inferred |f|)
                double log_pF = -1e99;
                double pF = 0.0;
                if(flag_len > 0) {
                    pF = (1.0/opts.sdev/std::sqrt(2.0 * M_PI)) * std::exp(-pow(flag_len - opts.fraglen, 2.0)/2.0/pow(opts.sdev, 2.0));
                    log_pF = std::log(pF);
                }
                // p(Y_k, I=0 | H_j)
                log_p_y_Hj[j] = logSumExp(log_p_y_Hj[j], log_pF + Y[k].score/params.tGet() + std::log(1.0-prob_disjoint));
                // p(Y_k, I=1 | H_j)
                log_p_y_Hj[j] = logSumExp(log_p_y_Hj[j], -std::log(2.0*params.gGet()) + Y[k].score/params.tGet()/2.0 + std::log(prob_disjoint));
            }
        }
        // denominator
        double Z = -1.0e99;
        for(size_t l=0; l<X.size(); ++l) {
            if(i_j[l] != -1) {
                Z = logSumExp(log_p_y_Hj[l] + log_p_Hj[l], Z);
            }
        }

        for(size_t j=0; j<X.size(); ++j) {
            p_Hj_y[j] = std::exp(log_p_y_Hj[j] + log_p_Hj[j] - Z);
        }
        // output
        /*std::cout << '\t';
        double best_prob = 0.0;
        int best_pair = -1;
        for(size_t j=0; j<X.size(); ++j) {
            if(p_Hj_y[j] > best_prob) {
                best_prob = p_Hj_y[j];
                best_pair = j;
            }
        }
        if(best_pair == -1 || i_j[best_pair] == -1) {
            std::cout << "-";
        } else {
            std::cout << "(" << X[best_pair].rName << "," << i_j[best_pair] << "," << X[best_pair].qStrand << "," << p_R * p_Hj_y[best_pair] << ")";
        }*/
        double best_prob = 0.0;
        int best_pair = -1;
        for(size_t j=0; j<X.size(); ++j) {
            if(p_Hj_y[j] > best_prob) {
                best_prob = p_Hj_y[j];
                best_pair = j;
            }
        }
        if(best_pair != -1 && i_j[best_pair] != -1) {
            alnpair.push_back({X[best_pair].rName, i_j[best_pair], X[best_pair].qStrand, p_R * p_Hj_y[best_pair], X[best_pair].rSeq[i-X[best_pair].qStart], best_pair});
        }
    }
    std::string seq;
    int alnPos = -1;
    bool isContinuousAlignment = false;
    for(size_t i=1; i<sizeX+1; ++i) {
        if(i < sizeX && alnpair[i-1].refName == alnpair[i].refName &&
           alnpair[i-1].refIndex + 1 == alnpair[i].refIndex &&
           alnpair[i-1].refStrand == alnpair[i].refStrand) {
            if(!isContinuousAlignment) {
                isContinuousAlignment = true;
                alnPos = alnpair[i-1].refIndex;
                seq += alnpair[i-1].qBase;
                seq += alnpair[i].qBase;
            } else {
                seq += alnpair[i].qBase;
            }
        } else {
            if(isContinuousAlignment) {
                std::cout << X[0].qName << '\t'            // QNAME
                          << seq.size() << 'M' << '\t'     // FLAG (temporary)
                          << alnpair[i-1].refName << '\t'  // RNAME
                          << alnPos << '\t'                // POS
                          << 255 << '\t'                   // MAPQ (unavailable)
                          << '*' << '\t'                   // CIGAR (temporary)
                          << '*' << '\t'                   // RNEXT (unavailable)
                          << 0 << '\t'                     // PNEXT (unavailable)
                          << 0 << '\t'                     // TLEN (unavailable)
                          << seq << '\t'                   // SEQ
                          << '*' << std::endl;             // QUAL (unavailable)
                std::cout << std::endl;
                seq.clear();
                isContinuousAlignment = false;
            }
        }
    }
}

void lastSplitPe(LastPairProbsOptions& opts) {
    const std::vector<std::string>& inputs = opts.inputFileNames;
    size_t n = inputs.size();
    std::ifstream inFile1;
    std::istream& input = (n > 0) ? cbrc::openIn(inputs[0], inFile1) : std::cin;
    AlignmentParameters params = readHeaderOrDie(input);
    int time = 0;
    while(input.good()) {
        std::cout << "execution: " << ++time << std::endl;
        std::vector<Alignment> X = readAlignmentSet(input);
        std::vector<Alignment> Y = readAlignmentSet(input);
        if(X.size() > 1) calcProbAndOutput(X, Y, params, opts);
        if(Y.size() > 1) calcProbAndOutput(Y, X, params, opts);
    }
}
