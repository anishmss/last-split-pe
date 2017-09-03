// Copyright 2017 Naruki Yoshikawa
// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

#include "last-split-pe.hh"

#include "io.hh"
#include "stringify.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <bitset>
#include <array>
#include <cassert>

static void err(const std::string &s)
{
    throw std::runtime_error(s);
}

static double logSumExp(const double a, const double b)
{
    // Adds numbers, in log space, to avoid overflow.
    const double m = std::max(a, b);
    return std::log(std::exp(a - m) + std::exp(b - m)) + m;
}

struct Alignment
{
    double score;
    std::string rName, qName, qNameNoPair, qPairName;
    long rStart, qStart, rAlnSize, qAlnSize, rSrcSize, qSrcSize;
    char rStrand, qStrand;
    std::string qSeq, rSeq, prob;
};

class AlignmentParameters
{
    // Parses the score scale factor, minimum score, and genome size.
    double t; // score scale factor
    double e; // minimum score
    long g;   // genome size
    long a;   // gap open
    long b;   // gap extention
    long x;   // drop parameter

  public:
    AlignmentParameters() : t(-1), e(-1), g(-1), a(-1), b(-1), x(-1) {} // dummy values

    void update(const std::string &line)
    {
        std::istringstream ss(line);
        std::string i;
        while (ss >> i)
        {
            const char *c = i.c_str();
            if (t == -1.0 && i.substr(0, 2) == "t=")
            {
                cbrc::unstringify(t, c + 2);
                if (t <= 0)
                    err("t must be positive");
            }
            if (e == -1.0 && i.substr(0, 2) == "e=")
            {
                cbrc::unstringify(e, c + 2);
                if (e <= 0)
                    err("e must be positive");
            }
            if (g == -1.0 && i.substr(0, 8) == "letters=")
            {
                cbrc::unstringify(g, c + 8);
                if (g <= 0)
                    err("letters must be positive");
            }
            if (a == -1 && i.substr(0, 2) == "a=")
            {
                cbrc::unstringify(a, c + 2);
                if (a <= 0)
                    err("a must be positive");
            }
            if (b == -1 && i.substr(0, 2) == "b=")
            {
                cbrc::unstringify(b, c + 2);
                if (b <= 0)
                    err("b must be positive");
            }
            if (x == -1 && i.substr(0, 2) == "x=")
            {
                cbrc::unstringify(x, c + 2);
                if (x <= 0)
                    err("x must be positive");
            }
        }
    }

    bool isValid() const { return t != -1 && e != -1 && g != -1; }

    void validate() const
    {
        if (t == -1)
            err("I need a header line with t=");
        if (e == -1)
            err("I need a header line with e=");
        if (g == -1)
            err("I need a header line with letters=");
    }

    double tGet() const { return t; }
    double eGet() const { return e; }
    long gGet() const { return g; }
    long aGet() const { return a; }
    long bGet() const { return b; }
    long xGet() const { return x; }
};

struct AlignmentPair
{
    std::string rName, qName;
    long rIndex, qIndex;
    char rStrand, qStrand;
    char rBase, qBase;
    double probability;
};

static AlignmentParameters readHeaderOrDie(std::istream &lines)
{
    std::string line;
    AlignmentParameters params;
    while (getline(lines, line))
    {
        if (line.substr(0, 1) == "#")
        {
            params.update(line);
            if (params.isValid())
                return params;
        }
        else if (line.find_first_not_of(" ") != std::string::npos)
        {
            break;
        }
    }
    params.validate(); // die
    return params;     // dummy
}

static Alignment readSingleAlignment(std::istream &input)
{
    // variables to store alignment information
    double score;
    std::string rName, qName, pairName, qNameNoPair, rSeq, qSeq, prob;
    char rStrand, qStrand;
    long rStart, qStart, rAlnSize, qAlnSize, rSrcSize, qSrcSize;

    std::string line; // placeholder for storing a line
    int n = 0;        // counter to distinguish ref/query
    while (std::getline(input, line))
    {
        if (line == "")
            break; // blank line

        char begin = line[0];
        std::stringstream ss(line);
        if (begin == 'a')
        {
            if (line.substr(0, 8) != "a score=")
            {
                std::cerr << "Bad File Format. Exit." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            ss.ignore(8); // ignore "a score="
            ss >> score;
        }
        else if (begin == 's')
        {
            char type;
            if (n == 0)
            {
                ss >> type >> rName >> rStart >> rAlnSize >> rStrand >> rSrcSize >> rSeq;
                n = 1;
            }
            else if (n == 1)
            {
                ss >> type >> qName >> qStart >> qAlnSize >> qStrand >> qSrcSize >> qSeq;
                //extracting pair membership info
                std::string delim = "/";
                size_t found = qName.find(delim);
                if (found == std::string::npos)
                {
                    std::cerr << "Read names not in ReadID/1  ReadID/2 format" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                qNameNoPair = qName.substr(0, found);
                pairName = qName.substr(found + delim.length());
                //std::cout << pairName << std::endl;
                if (pairName != "1" && pairName != "2")
                {
                    std::cerr << "Reads names should end in /1 or /2" << std::endl;
                    std::exit(EXIT_FAILURE);
                }

                n = 2;
            }
            else
            {
                std::cerr << "Reads names should end in /1 or /2" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else if (begin == 'p')
        {
            char type;
            ss >> type >> prob;
        }
    }
    Alignment aln;
    aln.score = score;
    aln.qName = qName;
    aln.qNameNoPair = qNameNoPair;
    aln.qPairName = pairName;
    aln.rName = rName;
    aln.qStart = qStart;
    aln.rStart = rStart;
    aln.rAlnSize = rAlnSize;
    aln.qAlnSize = qAlnSize;
    aln.rSrcSize = rSrcSize;
    aln.qSrcSize = qSrcSize;
    aln.rStrand = rStrand;
    aln.qStrand = qStrand;
    aln.qSeq = qSeq;
    aln.rSeq = rSeq;
    aln.prob = prob;

    return aln;
}

void outputInSAM(const std::vector<Alignment> &X, const std::vector<Alignment> &Y,
                 const std::string &refName, const long alnPos, const size_t totalLength,
                 const std::string &seq, std::stringstream &cigar, bool isSupplementary,
                 bool isFirst)
{
    // output in SAM format
    std::bitset<12> flag;
    flag.set(0, true);  // template having multiple segments in sequencing
    flag.set(1, false); // each segment properly aligned according to the aligner
    flag.set(2, false); // segment unmapped
    flag.set(3, false); // next segment in the template unmapped

    flag.set(4, X[0].qStrand == '-'); // SEQ being reverse complemented
    // What if some strand is '+'??
    // SEQ of the next segment in the template being reverse complemented
    // TODO: What if there are several alignments in Y?
    flag.set(5, Y[0].qStrand == '-' ? true : false);
    flag.set(6, isFirst);          // the first segment in the template
    flag.set(7, !isFirst);         // the last segment in the template
    flag.set(8, false);            // secondary alignment
    flag.set(9, false);            // not passing filters, such as platform/vendor quality controls
    flag.set(10, false);           // PCR or optical duplicate
    flag.set(11, isSupplementary); // supplementary alignment

    auto LastTwo = X[0].qName.substr(X[0].qName.length() - 2);
    std::string qName;
    if (LastTwo == "/1" || LastTwo == "/2")
    {
        qName = X[0].qName.substr(0, X[0].qName.length() - 2);
    }
    else
    {
        qName = X[0].qName;
    }
    if (totalLength < X[0].qSrcSize)
        cigar << X[0].qSrcSize - totalLength << "H";
    std::cout << qName << '\t'           // QNAME
              << flag.to_ulong() << '\t' // FLAG
              << refName << '\t'         // RNAME
              << alnPos + 1 << '\t'      // RPOS
              << 255 << '\t'             // MAPQ (unavailable)
              << cigar.str() << '\t'     // CIGAR
              << '*' << '\t'             // RNEXT (unavailable)
              << 0 << '\t'               // PNEXT (unavailable)
              << 0 << '\t'               // TLEN (unavailable)
              << seq << '\t'             // SEQ
              << '*' << std::endl;       // QUAL (unavailable)
}
void outputAlignmentSam(const std::vector<Alignment> &X, const std::vector<Alignment> &Y,
                        bool isFirst)
{
    assert(X.size() == 1);
    std::bitset<12> flag;
    // template having multiple segments in sequencing
    flag.set(0, true);
    // each segment properly aligned according to the aligner
    flag.set(1, false);
    // segment unmapped
    flag.set(2, false);
    // next segment in the template unmapped
    flag.set(3, false);
    // SEQ being reverse complemented
    flag.set(4, X[0].qStrand == '-');
    // SEQ of the next segment in the template being reverse complemented
    flag.set(5, Y[0].qStrand == '-' ? true : false);
    // the first segment in the template
    flag.set(6, isFirst);
    // the last segment in the template
    flag.set(7, !isFirst);
    // secondary alignment
    flag.set(8, false);
    // not passing filters, such as platform/vendor quality controls
    flag.set(9, false);
    // PCR or optical duplicate
    flag.set(10, false);
    // supplementary alignment
    flag.set(11, false);

    double prob = 0.0;
    for (size_t i = 0; i < X.size(); ++i)
    {
        char c = X[0].prob[i];

        // http://last.cbrc.jp/doc/last-split.html
        prob = std::max(prob, 1.0 - std::pow(10.0, -((c - 33 + 1) / 10.0))); // +1 is for lower bound
    }

    auto LastTwo = X[0].qName.substr(X[0].qName.length() - 2);
    std::string qName;
    if (LastTwo == "/1" || LastTwo == "/2")
    {
        qName = X[0].qName.substr(0, X[0].qName.length() - 2);
    }
    else
    {
        qName = X[0].qName;
    }

    std::cout << qName << '\t'                                    // QNAME
              << flag.to_ulong() << '\t'                          // FLAG
              << X[0].rName << '\t'                               // RNAME
              << X[0].rStart + 1 << '\t'                          // RPOS (since sam is 1-indexed, we need +1)
              << static_cast<int>(-10 * std::log10(prob)) << '\t' // MAPQ (TODO log of the probability of most reliable pair)
              << X[0].qSeq.size() << 'M' << '\t'                  // CIGAR
              << '=' << '\t'                                      // RNEXT (same, TODO: implement appropriately)
              << Y[0].rStart + 1 << '\t'                          // PNEXT (unavailable)
              << 0 << '\t'                                        // TLEN (unavailable)
              << X[0].qSeq << '\t'                                // SEQ
              << "*" << std::endl;                                // QUAL (unavailable)
}

std::vector<std::vector<AlignmentPair>> calcProb(std::vector<Alignment> &X, std::vector<Alignment> &Y,
                                                 AlignmentParameters &params, LastPairProbsOptions &opts)
{

    if (X.size() == 0 || Y.size() == 0)
    {
        std::cerr << "Invalid input: Not pair" << std::endl;
        //std::exit(EXIT_FAILURE);
        return std::vector<std::vector<AlignmentPair>>();
    }
    // vector to store result
    std::vector<std::vector<AlignmentPair>> alnprobs(X[0].qSrcSize, std::vector<AlignmentPair>(X.size()));

    // p(I=1) i.e. probability of disjoint
    const double prob_disjoint = 0.01;
    // the base x[i] aligns to position i_j in reference
    std::vector<long> i_j(X.size());
    // reference position
    std::vector<long> rPos(X.size(), 0);
    // reference position including gap
    std::vector<long> rPosWithGap(X.size(), 0);
    // query position including gap
    std::vector<long> qPosWithGap(X.size(), 0);
    // p(H_j | R): prior probability that x[i] has been sequenced from position i_j in reference
    std::vector<double> log_p_Hj(X.size());
    // p(y | H_j, R): likelihood
    std::vector<double> log_p_y_Hj(X.size());
    // p(H_j | y, R): posterior probability
    std::vector<double> p_Hj_y(X.size());

    const long NOALIGN = -1;
    const long GAP = -2;
    const double LOG0 = -1.0e99;

    for (size_t qPos = 0; qPos < alnprobs.size(); ++qPos)
    {
        // initialization
        std::fill(i_j.begin(), i_j.end(), NOALIGN);
        std::fill(log_p_Hj.begin(), log_p_Hj.end(), LOG0);
        std::fill(log_p_y_Hj.begin(), log_p_y_Hj.end(), LOG0);
        std::fill(p_Hj_y.begin(), p_Hj_y.end(), 0);

        double p_R = 0.0;
        for (size_t aln = 0; aln < X.size(); ++aln)
        {
            if (qPos == X[aln].qStart)
            {
                rPos[aln] = rPosWithGap[aln] = X[aln].rStart;
            }
            while (X[aln].qSeq[qPosWithGap[aln] - X[aln].qStart] == '-')
            {
                qPosWithGap[aln]++;
            }
        }
        for (size_t aln = 0; aln < X.size(); ++aln)
        {
            // if current query position is reported in X[aln]
            if (qPos >= X[aln].qStart && qPos < X[aln].qStart + X[aln].qAlnSize)
            {
                if (X[aln].rSeq[rPosWithGap[aln] - X[aln].rStart] == '-')
                {
                    i_j[aln] = GAP;
                }
                else
                {
                    i_j[aln] = rPos[aln];
                }
                char c = X[aln].prob[rPosWithGap[aln] - X[aln].rStart];
                double p = 1.0 - std::pow(10.0, -((c - 33) / 10.0));
                log_p_Hj[aln] = std::log(p);
                p_R += p;
            }
            // conditional probability p(Y_k | H_j)
            for (size_t k = 0; k < Y.size(); ++k)
            {
                // estimate flagment length
                int flag_len = -1;
                if (X[aln].qStrand == '+' && Y[k].qStrand == '-')
                {
                    flag_len = (qPos - X[aln].qStart - 1) + (Y[k].rStart - i_j[aln] + 1) + Y[k].qSrcSize;
                }
                else if (X[aln].qStrand == '+' && Y[k].qStrand == '-')
                {
                    flag_len = i_j[aln] - Y[k].rStart + i_j[aln] + X[aln].qSrcSize - (qPos - X[aln].qStart - 1);
                }
                // TODO: What if other cases??

                // p(inferred |f|)
                double log_pF = -1e99;
                double pF = 0.0;
                if (flag_len > 0)
                {
                    pF = (1.0 / opts.sdev / std::sqrt(2.0 * M_PI)) * std::exp(-pow(flag_len - opts.fraglen, 2.0) / 2.0 / pow(opts.sdev, 2.0));
                    log_pF = std::log(pF);
                }
                // p(Y_k, I=0 | H_j)
                log_p_y_Hj[aln] = logSumExp(log_p_y_Hj[aln], log_pF + Y[k].score / params.tGet() + std::log(1.0 - prob_disjoint));
                // p(Y_k, I=1 | H_j)
                log_p_y_Hj[aln] = logSumExp(log_p_y_Hj[aln],
                                            -std::log(2.0 * params.gGet()) + Y[k].score / params.tGet() / 2.0 + std::log(prob_disjoint));
            }
        }
        // denominator
        double Z = -1.0e99;
        for (size_t l = 0; l < X.size(); ++l)
        {
            if (i_j[l] != -1)
            {
                Z = logSumExp(log_p_y_Hj[l] + log_p_Hj[l], Z);
            }
        }
        for (size_t j = 0; j < X.size(); ++j)
        {
            p_Hj_y[j] = std::exp(log_p_y_Hj[j] + log_p_Hj[j] - Z);
        }
        for (size_t n = 0; n < X.size(); ++n)
        {
            alnprobs[qPos][n].rName = X[n].rName;
            alnprobs[qPos][n].qName = X[n].qName;
            alnprobs[qPos][n].rIndex = i_j[n];
            alnprobs[qPos][n].qIndex = qPos;
            alnprobs[qPos][n].qStrand = X[n].qStrand;
            alnprobs[qPos][n].rStrand = X[n].rStrand;
            if (qPos >= X[n].qStart && qPos < X[n].qStart + X[n].qAlnSize)
            {
                alnprobs[qPos][n].rBase = X[n].rSeq[rPosWithGap[n] - X[n].rStart];
                alnprobs[qPos][n].qBase = X[n].qSeq[qPosWithGap[n] - X[n].qStart];
            }
            else
            {
                alnprobs[qPos][n].rBase = 'X';
                alnprobs[qPos][n].qBase = 'X';
            }
            alnprobs[qPos][n].probability = p_R * p_Hj_y[n];
        }
        for (size_t aln = 0; aln < X.size(); ++aln)
        {
            qPosWithGap[aln]++;
            // when rSize and rAlnSize is different (e.g. 123/2) Need to check.
            if (rPos[aln] < X[aln].rStart || rPos[aln] >= X[aln].rStart + X[aln].rAlnSize
                 || X[aln].rSeq.at(rPosWithGap[aln] - X[aln].rStart) != '-')
            {
                rPos[aln]++;
            }
            rPosWithGap[aln]++;
        }
    }
    return alnprobs;
}

std::vector<AlignmentPair> chooseBestPair(std::vector<std::vector<AlignmentPair>> readProbs)
{
    std::vector<AlignmentPair> bestPairs;
    bestPairs.reserve(readProbs.size());

    for (auto &pairs : readProbs)
    {
        double bestProb = 0;
        int bestPair = -1;
        for (size_t i = 0; i < pairs.size(); ++i)
        {
            if (pairs[i].probability > bestProb)
            {
                bestProb = pairs[i].probability;
                bestPair = i;
            }
        }
        if(bestPair != -1) {
            bestPairs.push_back(pairs[bestPair]);
        } else { // if every pair probability is zero
            AlignmentPair e;
            e.probability = -1;
            e.rName = "NO ALIGN";
            bestPairs.push_back(e);
        }
    }
    return bestPairs;
}

void outputSAM(std::vector<AlignmentPair> read1Aln, std::vector<AlignmentPair> read2Aln)
{
    /*if(opts.isSamFormat) {
        std::string seq;
        long a = params.aGet();
        long b = params.bGet();
        long x = params.xGet();
        size_t idx = 1;
        bool isSupplementary = false;
        while(idx < alnpair.size()) {
            bool isEndAlignment = false;
            if(idx == alnpair.size()-1) isEndAlignment = true;
            // if match, start one alignment
            if(alnpair[idx-1].refName == alnpair[idx].refName &&
               alnpair[idx-1].refIndex == alnpair[idx].refIndex - 1 &&
               alnpair[idx-1].refStrand == alnpair[idx].refStrand) {
                size_t totalLength = 0;
                std::stringstream cigar;
                if(idx > 1) {
                    cigar << idx-1 << "H";
                    totalLength += idx - 1;
                }
                // count match length
                size_t matchLength = 2;
                seq = "";
                seq += alnpair[idx-1].qBase;
                seq += alnpair[idx].qBase;
                for(size_t i=idx+1; ; ++i) {
                    if(alnpair[i-1].refName == alnpair[i].refName &&
                       alnpair[i-1].refStrand == alnpair[i].refStrand) {
                        auto refIndexDiff = alnpair[i].refIndex-alnpair[i-1].refIndex;
                        if(refIndexDiff == 1) {
                            matchLength++;
                            seq += alnpair[i].qBase;
                        } else {
                            cigar << matchLength << "M";
                            totalLength += matchLength;
                            if(a + b*refIndexDiff < x) {
                                matchLength = 0;
                                cigar << refIndexDiff << "D";
                                seq += alnpair[i].qBase;
                            } else {
                                isEndAlignment = true;
                            }
                        }
                    } else {
                        bool isInsertion = false;
                        for(size_t j=i+1; a+b*(j-i)<x; ++j) {
                            auto refIndexDiff = alnpair[i].refIndex-alnpair[i-1].refIndex;
                            if(alnpair[j-1].refName == alnpair[j].refName &&
                               alnpair[j-1].refStrand == alnpair[j].refStrand &&
                               refIndexDiff == 1) 
                            {
                                cigar << j-i << "D";
                                totalLength += j-i;
                                i = j;
                                isInsertion = true;
                                break;
                            }
                        }
                        if(!isInsertion) {
                            isEndAlignment = true;
                        }
                    }
                    if(i == alnpair.size()-1) {
                        isEndAlignment = true;
                        cigar << matchLength << "M";
                        totalLength += matchLength;
                    }
                    if(isEndAlignment) {
                        outputInSAM(X, Y, alnpair[idx].refName, alnpair[idx-1].refIndex, totalLength, seq, cigar, isSupplementary, isFirst);
                        idx = i + 1;
                        isSupplementary = true;
                        break;
                    }
                }
            } else {
                idx++;
            }
        }
    } else {
        std::cout << std::endl;
    }*/
}

void outputNative(std::vector<AlignmentPair> readAln)
{
    // invalid input
    if(readAln.size() == 0) {
        return;
    }
    std::cout << readAln[0].qName;
    for (auto &p : readAln)
    {
        std::cout << "\t(" << p.rName << "," << p.rIndex << "," << p.rStrand << ")";
    }
    std::cout << std::endl;
}

void startSplitPEProcess(std::vector<Alignment> &alns1, std::vector<Alignment> &alns2,
                         AlignmentParameters &params, LastPairProbsOptions &opts)
{

    //std::cout << "start calcProb" << std::endl;
    std::vector<std::vector<AlignmentPair>> read1Probs = calcProb(alns1, alns2, params, opts);
    std::vector<std::vector<AlignmentPair>> read2Probs = calcProb(alns2, alns1, params, opts);
    //std::cout << "finish calcProb" << std::endl;
    //std::cout << "start chooseBestPair" << std::endl;
    std::vector<AlignmentPair> read1FinalAln = chooseBestPair(read1Probs);
    std::vector<AlignmentPair> read2FinalAln = chooseBestPair(read2Probs);
    //std::cout << "end chooseBestPair" << std::endl;

    //std::cout << "start output" << std::endl;
    if (opts.isSamFormat)
    {
        outputSAM(read1FinalAln, read2FinalAln);
    }
    else
    {
        outputNative(read1FinalAln);
        outputNative(read2FinalAln);
    }
    //std::cout << "end output" << std::endl;
}

void lastSplitPe(LastPairProbsOptions &opts)
{
    //std::cout << "start lastSplitPe" << std::endl;
    const std::vector<std::string> &inputs = opts.inputFileNames;
    size_t n = inputs.size();
    std::ifstream inFile1;
    std::istream &input = (n > 0) ? cbrc::openIn(inputs[0], inFile1) : std::cin;
    AlignmentParameters params = readHeaderOrDie(input);

    std::vector<Alignment> X, Y; // X will hold alignments of read/1, and Y of read/2

    Alignment currentAln = readSingleAlignment(input);
    if (currentAln.qName == "")
    {
        std::cerr << "Error: maf file is empty" << std::endl;
        return; // empty maf file
    }
    while (true)
    {
        //std::cout << currentAln.qName << " " << currentAln.qNameNoPair << " " << currentAln.qPairName << std::endl;
        if (currentAln.qPairName == "1") {
            X.push_back(currentAln);
        } else {
            Y.push_back(currentAln);
        }

        Alignment nextAln = readSingleAlignment(input);

        if (nextAln.qName != "")
        {  //some alignment was read
            if (currentAln.qNameNoPair == nextAln.qNameNoPair)
            {   //either same read OR same pair
                //std::cout << "same pair" << std::endl;
                currentAln = nextAln;
                continue;
            }
            else
            {   //different pair
                //std::cout << "different pair" << std::endl;
                startSplitPEProcess(X, Y, params, opts);
                currentAln = nextAln;
                X.clear();
                Y.clear();
                continue;
            }
        }
        else
        {   // readSingleAlignment didn't read in alignment because of EOF or other reasons.
            //std::cout << "start!" << std::endl;
            startSplitPEProcess(X, Y, params, opts);
            //std::cout << "EOF" << std::endl;
            break;
        }
    }
}