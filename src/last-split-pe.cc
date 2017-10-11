// Copyright 2017 Naruki Yoshikawa, Anish MS Shrestha
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

const long NOALIGN = -1;
const long GAP = -2;
const double LOG0 = -1.0e99;
const double probZero = 1.0e-10;

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
    std::string qSeq, rSeq, prob, qual;
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
    std::string rName, qName, qNameNoPair, qPairName;
    long rIndex, qIndex;
    char rStrand, qStrand;
    char rBase, qBase;
    double probability;
};

struct SAMaln //stores information for one line of SAM of a read, minus alignment info of its mate 
{
    std::string qNameNoPair, refName, querySeq, cigar;
    long refStartPos, queryStartPos;
    int tlen; // template length = fragment length
    char queryStrand;
    bool isSupplementary,isFirst;
    double errProb;
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
    std::string rName, qName, pairName, qNameNoPair, rSeq, qSeq, prob, qual;
    char rStrand, qStrand;
    long rStart, qStart, rAlnSize, qAlnSize, rSrcSize, qSrcSize;
    bool read_a, read_s1, read_s2, read_q, read_p ;
    read_a = read_s1 = read_s2 = read_q = read_p = false;

    std::string line; // placeholder for storing a line
    int n = 0;        // counter to distinguish ref/query
    while (std::getline(input, line)) 
    {
        if (line == "")
            continue; // allow multiple blank lines between records

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
            read_a = true;
        }
        else if (begin == 's')
        {
            char type;
            if (n == 0)
            {
                ss >> type >> rName >> rStart >> rAlnSize >> rStrand >> rSrcSize >> rSeq;
                read_s1 = true;
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
                read_s2 = true;
                n = 2;
            }
            else
            {
                std::cerr << "Malformed maf file" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else if (begin == 'q')
        {
            char type;
            ss >> type >> qual ;
            read_q = true;
        }
        else if (begin == 'p')
        {
            char type;
            ss >> type >> prob;
            if (!(read_a && read_s1 && read_s2))
            {
                std::cerr << "Malformed maf file" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
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
    if (read_q) aln.qual = qual ;
    return aln;
    
}

std::vector<std::vector<AlignmentPair>> calcProb(std::vector<Alignment> &X, std::vector<Alignment> &Y,
                                                 AlignmentParameters &params, LastPairProbsOptions &opts, bool verbose)
{
    
    /*
    if (X.size() == 0 || Y.size() == 0)
    {
        std::cerr << "Invalid input: Not pair" << std::endl;
        //std::exit(EXIT_FAILURE);
        return std::vector<std::vector<AlignmentPair>>();
    }
    */
    // vector to store result
    std::vector<std::vector<AlignmentPair>> alnprobs(X[0].qSrcSize, std::vector<AlignmentPair>(X.size()));

    
    //construct query position to alignment column map. 
    //query position refers to position in actual read. for -ve strand alignments, we need to be careful.
    std::vector<std::vector<long>> qPos2AlnCol(X.size(),std::vector<long>(X[0].qSrcSize,NOALIGN)); //long is overkill, but for consistency
    for (size_t aln = 0; aln < X.size(); ++aln)
    {
        auto queryStrand = X[aln].qStrand;
        auto queryPosInAln = X[aln].qStart; //caution: if qStrand is '-', this will not reflect the actual position on the query.
        auto querySeq = X[aln].qSeq;
       
        for(long alnColPos=0; alnColPos < querySeq.length(); alnColPos++) // querySeq.length() gives alignment size. 
        {
            if (querySeq[alnColPos] != '-')
            {
                if (queryStrand == '+') qPos2AlnCol[aln][queryPosInAln] = alnColPos; //queryPosInAln is same as queryPosition
                else qPos2AlnCol[aln][X[0].qSrcSize-1-queryPosInAln] = alnColPos; //queryPos needs to be computed from queryPosInAln 
                queryPosInAln++;
            }
        }

    }
    
    /*
    std::cout << "qPos2AlnCol" << std::endl ;
    for (size_t aln = 0; aln < X.size(); ++aln)
    {
        std::cout << "aln " << aln << ":" ;
        for ( int i=0; i < qPos2AlnCol[aln].size(); i++)
        {
            std::cout << qPos2AlnCol[aln][i] << "," ;
        }
        std::cout << std::endl ;
    }
    */
    
    
    //construct alignment column to reference position map
    std::vector<std::vector<long>> alnCol2rPos; 
    for (size_t aln = 0; aln < X.size(); ++aln)
    {
        auto refPos = X[aln].rStart;
        auto refSeq = X[aln].rSeq;
        std::vector<long> tempMap;
        for(long alnColPos=0; alnColPos < refSeq.length(); alnColPos++)
        {
            if (refSeq[alnColPos] != '-')
            {
                tempMap.push_back(refPos);
                refPos++;
            }
            else
            {
                tempMap.push_back(GAP);
            }
            
        }
        alnCol2rPos.push_back(tempMap);
    }
    
    /*
    std::cout << "alnCol2rPos" << std::endl ;
    for (size_t aln = 0; aln < X.size(); ++aln)
    {
        std::cout << "aln " << aln << ":" ;
        for ( int i=0; i < alnCol2rPos[aln].size(); i++)
        {
            std::cout << alnCol2rPos[aln][i] << "," ;
        }
        std::cout << std::endl ;
    }
    */
    
    //the following is the core part of last-split-pe where probabilities are updated for each position in query

    // p(I=1) i.e. probability of disjoint
    const double prob_disjoint = 0.01;
    // the base x[i] aligns to position i_j in reference
    std::vector<long> i_j(X.size());
    
    /*
    // reference position
    std::vector<long> rPos(X.size(), 0);
    // reference position including gap
    std::vector<long> rPosWithGap(X.size(), 0);
    // query position including gap
    std::vector<long> qPosWithGap(X.size(), 0);
    */
    
    // p(H_j | R): prior probability that x[i] has been sequenced from position i_j in reference
    std::vector<double> log_p_Hj(X.size());
    // p(y | H_j, R): likelihood
    std::vector<double> log_p_y_Hj(X.size());
    // p(H_j | y, R): posterior probability
    std::vector<double> p_Hj_y(X.size());

    for (int qPos = 0; qPos < alnprobs.size(); ++qPos) 
    {   
        if (verbose) std::cout << "readPosition:" << qPos << ": " << std::endl;
        // initialization
        std::fill(i_j.begin(), i_j.end(), NOALIGN);
        std::fill(log_p_Hj.begin(), log_p_Hj.end(), LOG0);
        std::fill(log_p_y_Hj.begin(), log_p_y_Hj.end(), LOG0);
        std::fill(p_Hj_y.begin(), p_Hj_y.end(), 0);

        
        //set up prior probability
        double p_R = 0.0; // sum of column probabilities across candidate alignments. updated below
        for (size_t aln = 0; aln < X.size(); ++aln)
        {
            // if current query position is reported in X[aln]
            //if (qPos >= X[aln].qStart && qPos < X[aln].qStart + X[aln].qAlnSize)
            if (qPos2AlnCol[aln][qPos] != NOALIGN)
            {
                /*
                if (X[aln].rSeq[rPosWithGap[aln] - X[aln].rStart] == '-')
                {
                    i_j[aln] = GAP;
                }
                else
                {
                    i_j[aln] = rPos[aln];
                }
                */
                i_j[aln] = alnCol2rPos[aln][qPos2AlnCol[aln][qPos]] ;              
                
                //char c = X[aln].prob[rPosWithGap[aln] - X[aln].rStart]; 
                char c = X[aln].prob[qPos2AlnCol[aln][qPos]];
                double p = 1.0 - std::pow(10.0, -((c - 33) / 10.0)); 
                log_p_Hj[aln] = std::log(p);
                p_Hj_y[aln] = p ; //initialize posterior = prior.
                p_R += p;
                
            }
        }
        
        if (verbose)
        {
        
            //checking prior probs
            std::cout << "Prior probs" << std::endl;
            for (size_t aln = 0; aln < X.size(); ++aln)
            {
                std::cout << std::exp(log_p_Hj[aln]) << "\t" ; 
            }
            std::cout << std::endl;
            std::cout << "p_R:" << p_R << std::endl;

        }
        
        if ( p_R > 1 ) 
        {   
            //TODO: first rescale the priors so that they add up to 1?
            p_R = 1.0; //not sure how LAST discretizes its probability
        }
        //if i_j[aln] == GAP for one of the alns, we don't know what to do, yet.
        bool oneGAP = false;
        for (size_t aln = 0; aln < X.size(); ++aln)
        {
            if (i_j[aln] == GAP) 
            {
                oneGAP = true;
                break;
            }
        }
        
        if (!(X.size() == 1 || Y.empty() || oneGAP)) //update probabilities only if need be. makes code faster(?). last condition is TODO
        {   
            //std::cout << "updating probs" << std::endl;
            for (size_t aln = 0; aln < X.size(); ++aln)
            {
                // conditional probability p(Y_k | H_j)
                for (size_t k = 0; k < Y.size(); ++k)
                {
                    
                    // estimate fragment length
                    int frag_len = -1;
                    if (X[aln].qStrand == '+' && Y[k].qStrand == '-')
                    {
                        //rlag_len = (qPos - X[aln].qStart - 1) + (Y[k].rStart - i_j[aln] + 1) + Y[k].qSrcSize;
                        frag_len =    qPos + (Y[k].rStart - i_j[aln] + 1) + (Y[k].qSrcSize - Y[k].qStart-1) ; 
                        //            (a)  +           (b)                +  (c)  
                        //  (a): length of 5'-end portion before qPos,  
                        //  (b): distance on reference  between mapped position of qPos and the start of local alignment Y[k] of mate
                        //  (c): length of remaining portion of mate, starting from Y[k].qStart
                        //is this correct? or the one above?
                        if (verbose) std::cout << "fragment length:" << frag_len << std::endl ;
                        
                    }
                    else if (X[aln].qStrand == '-' && Y[k].qStrand == '+')
                    {
                        //frag_len = i_j[aln] - Y[k].rStart + i_j[aln] + X[aln].qSrcSize - (qPos - X[aln].qStart - 1); // this looks wrong because i_j appears twice 
                        frag_len = qPos + (i_j[aln] - Y[k].rStart+1) + Y[k].qStart  ;
                        if (verbose) std::cout << "fragment length:" << frag_len << std::endl ;
                    }
                    // TODO: What if other cases??
    
                    // p(inferred |f|)
                    double log_pF = -1e99;
                    double pF = 0.0;
                    if (frag_len > 0)
                    {
                        pF = (1.0 / opts.sdev / std::sqrt(2.0 * M_PI)) * std::exp(-pow(frag_len - opts.fraglen, 2.0) / 2.0 / pow(opts.sdev, 2.0));
                        log_pF = std::log(pF);
                        
                        if (verbose) std::cout << "fragment prob: " << pF << std::endl;
                    }
                    // p(Y_k, I=0 | H_j)
                    log_p_y_Hj[aln] = logSumExp(log_p_y_Hj[aln], log_pF + Y[k].score / params.tGet() + std::log(1.0 - prob_disjoint));
                    // p(Y_k, I=1 | H_j)
                    /*log_p_y_Hj[aln] = logSumExp(log_p_y_Hj[aln],
                                                -std::log(2.0 * params.gGet()) + Y[k].score / params.tGet() / 2.0 + std::log(prob_disjoint));*/
                    log_p_y_Hj[aln] = logSumExp(log_p_y_Hj[aln],
                                                -std::log(2.0 * params.gGet()) + Y[k].score / params.tGet() + std::log(prob_disjoint)); // removing 2.0 in the second expression . why is it there?
                }
            }
            
            if (verbose)
            {
                std::cout <<"Checking likelihood" << std::endl;
                for (size_t aln = 0; aln < X.size(); ++aln)
                {
                    std::cout <<  log_p_y_Hj[aln] << "\t" ; 
                }
                std::cout << std::endl;
            }
            
            
            //denominator
            double Z = -1.0e99;
            for (size_t l = 0; l < X.size(); ++l)
            {
                if (i_j[l] != NOALIGN)
                {
                    Z = logSumExp(log_p_y_Hj[l] + log_p_Hj[l], Z);
                }
            }
            for (size_t j = 0; j < X.size(); ++j)
            {
                //why is there no if (i_j[l] != NOALIGN) check like above?
                p_Hj_y[j] = std::exp(log_p_y_Hj[j] + log_p_Hj[j] - Z);
            }
        }
        
        
        if (verbose)
        {
            std::cout <<"Checking unscaled posterior" << std::endl;
            for (size_t aln = 0; aln < X.size(); ++aln)
            {
                std::cout <<  p_Hj_y[aln] << "\t" ; 
            }
            std::cout << std::endl;
        }
        for (size_t aln = 0; aln < X.size(); ++aln)
        {
            alnprobs[qPos][aln].rName = X[aln].rName;
            alnprobs[qPos][aln].qName = X[aln].qName;
            alnprobs[qPos][aln].qNameNoPair = X[aln].qNameNoPair;
            alnprobs[qPos][aln].qPairName = X[aln].qPairName;
            alnprobs[qPos][aln].rIndex = i_j[aln];
            alnprobs[qPos][aln].qIndex = qPos;
            alnprobs[qPos][aln].qStrand = X[aln].qStrand;
            alnprobs[qPos][aln].rStrand = X[aln].rStrand;
            alnprobs[qPos][aln].qBase = X[aln].qSeq[qPos2AlnCol[aln][qPos]];
            
            if (qPos2AlnCol[aln][qPos] != NOALIGN)
            {
                //alnprobs[qPos][n].rBase = X[n].rSeq[rPosWithGap[n] - X[n].rStart];
                alnprobs[qPos][aln].rBase = X[aln].rSeq[qPos2AlnCol[aln][qPos]];
            }
            else
            {
                alnprobs[qPos][aln].rBase = 'X';
                //alnprobs[qPos][n].qBase = 'X';
            }
            alnprobs[qPos][aln].probability = p_R * p_Hj_y[aln];
            
            
        }
    }
    return alnprobs;
}

//std::vector<AlignmentPair>
void chooseBestPair(const std::vector<std::vector<AlignmentPair>> &readProbs, std::vector<AlignmentPair> &bestPairs)
{
    //std::vector<AlignmentPair> bestPairs;
    
    if (readProbs.empty()) return;
    
    bestPairs.reserve(readProbs.size());
    

    for (auto &pairs : readProbs)
    {
        double bestProb = 0;
        int bestPair = -1;
        for (size_t i = 0; i < pairs.size(); ++i)
        {
            if ((pairs[i].rIndex != NOALIGN) && (pairs[i].probability > bestProb))
            {
                bestProb = pairs[i].probability;
                bestPair = i;
            }
        }
        if(bestPair != -1) {
            bestPairs.push_back(pairs[bestPair]);
        } 
        else if (bestProb == 0){
           AlignmentPair noAlignPair; 
           noAlignPair.rIndex = NOALIGN;
           bestPairs.push_back(noAlignPair);
        }
        
        else { // if every pair was NOALIGN, any one will do
            bestPairs.push_back(pairs[0]); 
            /*
            AlignmentPair e;
            e.probability = -1;
            e.rName = "-";
            e.qName = readProbs[0][0].qName ;
            e.rIndex = NOALIGN;
            e.rStrand='x';
            e.qBase = pairs[i].qBase;
            bestPairs.push_back(e);
            */
        }
    }
    
}

/*
void outputInSAM(const std::string qNameNoPair,
                 const std::string &refName, const long alnPos, const char qStrand, 
                 const std::string &seq, std::stringstream &cigar, bool isSupplementary,
                 bool isFirst , const std::vector<AlignmentPair> &Y)
{
    // output in SAM format
    std::bitset<12> flag;
    flag.set(0, true);  // template having multiple segments in sequencing
    flag.set(1, false); // each segment properly aligned according to the aligner
    flag.set(2, false); // segment unmapped
    flag.set(3, false); // next segment in the template unmapped

    flag.set(4, qStrand == '-'); // SEQ being reverse complemented
    // What if some strand is '+'??
    // SEQ of the next segment in the template being reverse complemented
    // TODO: What if there are several alignments in Y?
    //flag.set(5, Y[0].qStrand == '-' ? true : false);
    flag.set(5,false); //TODO
    flag.set(6, isFirst);          // the first segment in the template
    flag.set(7, !isFirst);         // the last segment in the template
    flag.set(8, false);            // secondary alignment
    flag.set(9, false);            // not passing filters, such as platform/vendor quality controls
    flag.set(10, false);           // PCR or optical duplicate
    flag.set(11, isSupplementary); // supplementary alignment

   
    std::cout << qNameNoPair << '\t'     // QNAME
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
*/



void writeSAMoutput(const SAMaln &readSAM,  
                    const std::string &rnext, const long &pnext, const bool &nextSegmentUnmapped, const bool &nextSegmentRevComp, const bool &isProperPair)
{
    // output in SAM format
    std::bitset<12> flag;
    flag.set(0, true);  // template having multiple segments in sequencing
    flag.set(1, isProperPair); // each segment properly aligned according to the aligner
    flag.set(2, false); // segment unmapped
    flag.set(3, nextSegmentUnmapped); // next segment in the template unmapped
    flag.set(4, readSAM.queryStrand == '-'); // SEQ being reverse complemented
    flag.set(5, nextSegmentRevComp); // mate reverse complemented
    flag.set(6, readSAM.isFirst);          // the first segment in the template
    flag.set(7, !readSAM.isFirst);         // the last segment in the template
    flag.set(8, false);            // secondary alignment
    flag.set(9, false);            // not passing filters, such as platform/vendor quality controls
    flag.set(10, false);           // PCR or optical duplicate
    flag.set(11, readSAM.isSupplementary); // supplementary alignment

   
    std::cout << readSAM.qNameNoPair << '\t'     // QNAME
              << flag.to_ulong() << '\t' // FLAG
              << readSAM.refName << '\t'         // RNAME
              << readSAM.refStartPos << '\t'      // RPOS
              << static_cast<int>(-10 * std::log10(readSAM.errProb)) << '\t'             // MAPQ -10 log_10 (prob alignment is wrong)
              //<< readSAM.bestProb << '\t' 
              << readSAM.cigar << '\t'     // CIGAR
              << rnext << '\t'             // RNEXT (unavailable)
              << pnext << '\t'               // PNEXT (unavailable)
              << readSAM.tlen << '\t'               // TLEN (unavailable)
              << readSAM.querySeq << '\t'             // SEQ
              << '*' << std::endl;       // QUAL (unavailable)
}

void buildSAMsegmentPlus(size_t &i, const std::string &qNameNoPair, std::vector<SAMaln> &samStructVec, 
                         const std::vector<AlignmentPair> &alnPair,
                         const long &a, const long &b, const long &x, 
                         const bool &isFirst, const bool &isSupplementary, const bool &verbose)
{
    if (verbose) std::cout << "inside buildSAMsegmentPlus" << std::endl;
    
    //entering actual segment
    long refStartPos = alnPair[i].rIndex; 
    //long refEndPos=refStartPos; //the last ref position in this SAM line. (needed for - strand alignments)
    long queryStartPos = i;
    long queryEndPos = queryStartPos;
    std::string refName = alnPair[i].rName;
    char queryStrand = alnPair[i].qStrand;
    std::string querySeq="";
    std::vector<std::string> preCigar; // e.g. cigar 4M2D3M is represented as  preCigar [M,M,M,M,D,D,M,M,M]. this allows cleaner looping
    while(i < alnPair.size()-1) // looping inside this SAM line. only the final exit condition is stated here. follow the explicit breaks and continues.
    {   
        
        preCigar.push_back("M") ; 
        // The starting point will always be an M position. Because:
        // If we are entering first time as new segment, we can't be at NOALIGN or GAP because of \cite{getToFirstM} \.
        // If we are here because of looping inside the same segment, we can't be at GAP because of \label{conditionGAP} below, 
        // and we don't allow same segment if we encountered NOALIGN.
        //refEndPos = alnPair[i].rIndex;
        querySeq += alnPair[i].qBase;
        queryEndPos = i;
        
        if (alnPair[i+1].rIndex != GAP && alnPair[i+1].rIndex != NOALIGN )
        {
            if (verbose) std::cout << "at pos" << i+1 << ": ALIGNED" << std::endl;
            if (alnPair[i].rName == alnPair[i+1].rName &&
                alnPair[i+1].qStrand == '+') 
            {
                if (verbose) std::cout << "sub-case: same chr, same strand" << std::endl;
                auto refIndexDiff = alnPair[i+1].rIndex-alnPair[i].rIndex;
                if (verbose) std::cout << "refIndexDiff: " << refIndexDiff << std::endl;
                if(refIndexDiff == 1) 
                {
                    i++;
                    continue; //same segment
                }
                else if((refIndexDiff > 0) && (a + b*refIndexDiff < x))
                {
                    for(int j=1; j<refIndexDiff; j++)
                        preCigar.push_back("D");
                    i++;
                    continue; //same segment
                }
                else // very large deletion, or some arrangment that alters left-to-right order
                {
                    i++;
                    break; //new segment
                }
               
            }
            else{ //different chr, or different strand
                if (verbose) std::cout << "sub-case: different chr , different strand " << std::endl;
                i++;
               break; //new segment
            }
        }
        else if (alnPair[i+1].rIndex == GAP) //\label{conditionGAP}
        {
            if(verbose) std::cout << "at pos" << i+1 << ": GAP" << std::endl;
            auto gapStart = i+1; //remember GAP start position
            while( alnPair[i+1].rIndex == GAP) // contiguous GAPs?  can't overshoot alnpair because of NOALIGN dummy alignment at the end
            {
                i++;
            }
            // i  now points to the last GAP position
            if (verbose) std::cout << "last GAP pos: " << i << std::endl;
            if (alnPair[i+1].rIndex == NOALIGN) 
            {
                if (verbose) std::cout << "Gap ended by NOALIGN" << std::endl;
                i++;
                break; 
            }         
            if (alnPair[i+1].rIndex-alnPair[gapStart-1].rIndex==1) //i is an M position. this is an insertion.
            {
                for (int m=gapStart; m<=i ; m++) //i=gapEnd
                {
                    if (verbose) std::cout << "Insert an I" << std::endl;
                    preCigar.push_back("I");
                    querySeq += alnPair[m].qBase;
                }
                i++;
                continue; //same segment. 
            } 
            else
            { //not an insertion. maybe those infamous insertion plus deletion.
                if (verbose) std::cout << "Not an insertion. Start new segment" << std::endl; 
                i++;
                break; //new segment. 
            }
        }
        else if (alnPair[i+1].rIndex == NOALIGN) // condition checking not required, but explicit is better than implicit.
        {
            if (verbose) std::cout << "at pos" << i+1 << ": NOALIGN" << std::endl;
            i++;
            break;
        }
        
    }
    //exiting SAMSegment. i points to start of next segment 
    
    //gather what we have so far, create a SAMaln struct,  and add to samStructVec
    SAMaln tempSAM;     
    if (verbose) 
    {
        std::cout << "preCigar:" ;
        for (auto pC: preCigar)
        {
            std::cout << pC ;
        }
        std::cout << std::endl;
    }
    
    //error probability
    //get error prob among alignment columns of this segment
    double bestProb = 0.0 ;
    if (verbose) std::cout << "Column probabilities" << std::endl;
    for (size_t r = queryStartPos; r <= queryEndPos ; r++ )
    {
        if (verbose) std::cout << alnPair[r].probability << std::endl; 
        if (alnPair[r].probability > bestProb) bestProb = alnPair[r].probability;
    }
    if (verbose) std::cout << std::endl;
    assert(bestProb > 0.0);
    double errProb = 1 - bestProb ; 
    if (errProb == 0) errProb = probZero;
    if (verbose) std::cout << "Errprob" << errProb << std::endl;
         
    
    //cigar string
    std::stringstream cigar;
    if (queryStartPos > 0) cigar << queryStartPos << "H" ;
    std::string currentState = preCigar[0];
    int length=0;
    for (auto &state:preCigar)
    {
        if (currentState == state)
        {
            length++;
        }
        else
        {
            cigar << length << currentState;
            length = 1;
            currentState = state;
        }
    }
    cigar << length << currentState;
    if (queryEndPos < alnPair.size()-2 )  //-2 because of dummy align at end
    {
        cigar << alnPair.size()-2-queryEndPos << "H";
    }

    tempSAM.cigar = cigar.str();
    tempSAM.errProb = errProb ; 
    tempSAM.isFirst = isFirst;
    tempSAM.isSupplementary = isSupplementary;
    tempSAM.qNameNoPair = qNameNoPair;
    tempSAM.querySeq = querySeq; 
    tempSAM.queryStartPos = queryStartPos; //on the read
    tempSAM.queryStrand = queryStrand;
    tempSAM.refName = refName;
    tempSAM.refStartPos = refStartPos+1; //SAM is 1-based 
    tempSAM.tlen = 0;
    
            
    samStructVec.push_back(tempSAM);
}


void buildSAMsegmentMinus(size_t &i, const std::string &qNameNoPair, std::vector<SAMaln> &samStructVec, 
                         const std::vector<AlignmentPair> &alnPair,
                         const long &a, const long &b, const long &x, 
                         const bool &isFirst, const bool &isSupplementary, const bool &verbose)
{
    if (verbose) std::cout << "inside buildSAMsegmentMinus" << std::endl;
    
    //entering actual segment
    long refStartPos = alnPair[i].rIndex; 
    long refEndPos=refStartPos; //the last ref position in this SAM line. (needed for - strand alignments)
    long queryStartPos = i;
    long queryEndPos = queryStartPos;
    std::string refName = alnPair[i].rName;
    char queryStrand = alnPair[i].qStrand;
    std::string querySeq="";
    std::vector<std::string> preCigar; // e.g. cigar 4M2D3M is represented as  preCigar [M,M,M,M,D,D,M,M,M]. this allows cleaner looping
    
 
    while(i < alnPair.size()-1) // looping inside this SAM line. only the final exit condition is stated here. to understand logic, follow the explicit breaks and continues .
    {   
        
        preCigar.push_back("M") ; 
        // The starting point will always be an M position. Because:
        // If we are entering first time as new segment, we can't be at NOALIGN or GAP because of \cite{getToFirstM} \.
        // If we are here because of looping inside the same segment, we can't be at GAP because of \label{conditionGAP} below, 
        // and we don't allow same segment if we encountered NOALIGN.
        refEndPos = alnPair[i].rIndex;
        queryEndPos = i;
        querySeq += alnPair[i].qBase;
        
        if (alnPair[i+1].rIndex != GAP && alnPair[i+1].rIndex != NOALIGN )
        {
            if (verbose) std::cout << "at pos" << i+1 << ": ALIGNED" << std::endl;
            if (alnPair[i].rName == alnPair[i+1].rName &&
                alnPair[i+1].qStrand == '-') 
            {
                if (verbose) std::cout << "sub-case: same chr, same strand" << std::endl;
                auto refIndexDiff = alnPair[i].rIndex-alnPair[i+1].rIndex; //i+1.rIndex is smaller because of - strand
                if (verbose) std::cout << "refIndexDiff: " << refIndexDiff << std::endl; 
                if(refIndexDiff == 1) 
                {
                    i++;
                    continue; //same segment
                }
                else if((refIndexDiff > 0) && (a + b*refIndexDiff < x))
                {
                    for(int j=1; j<refIndexDiff; j++)
                        preCigar.push_back("D");
                    i++;
                    continue; //same segment
                }
                else // very large deletion, or some arrangment that alters left-to-right order
                {
                    i++;
                    break; //new segment
                }
               
            }
            else{ //different chr, or different strand
                if (verbose) std::cout << "sub-case: different chr , different strand " << std::endl;
                i++;
               break; //new segment
            }
        }
        else if (alnPair[i+1].rIndex == GAP) //\label{conditionGAP}
        {
            //now we need to look forward until the end of the gap to decide what to do
            if(verbose) std::cout << "at pos" << i+1 << ": GAP" << std::endl;
            auto gapStart = i+1; //remember GAP start position
            while( alnPair[i+1].rIndex == GAP) // contiguous GAPs?  can't overshoot alnpair because of NOALIGN dummy alignment at the end
            {
                i++;
            }
            // i  now points to the last GAP position
            if (verbose) std::cout << "last GAP pos: " << i << std::endl;
            if (alnPair[i+1].rIndex == NOALIGN) 
            {
                if (verbose) std::cout << "Gap ended by NOALIGN" << std::endl;
                i++;
                break; 
            }         
            if (alnPair[gapStart-1].rIndex-alnPair[i+1].rIndex==1) //i is an M position. this is an insertion.
            {
                for (int m=gapStart; m<=i ; m++) //i=gapEnd
                {
                    if (verbose) std::cout << "Insert an I" << std::endl;
                    preCigar.push_back("I");
                    querySeq += alnPair[m].qBase;
                }
                i++;
                continue; //same segment. 
            } 
            else
            { //not an insertion. maybe those infamous insertion plus deletion.
                if (verbose) std::cout << "Not an insertion. Start new segment" << std::endl; 
                i++;
                break; //new segment. 
            }
        }
        else if (alnPair[i+1].rIndex == NOALIGN) // condition checking not required, but explicit is better than implicit.
        {
            if (verbose) std::cout << "at pos" << i+1 << ": NOALIGN" << std::endl;
            i++;
            break;
        }
        
    }
    //exiting SAMSegment. i points to start of next segment.
   
    //gather what we have so far, create a SAMaln struct,  and add to samStructVec
    SAMaln tempSAM;     
    if (verbose) 
    {
        std::cout << "preCigar:" ;
        for (auto pC: preCigar)
        {
            std::cout << pC ;
        }
        std::cout << std::endl;
    }
    
    //error probability
    //get error prob among alignment columns of this segment
    double bestProb = 0.0 ;
    if (verbose) std::cout << "Column probabilities" << std::endl;
    for (size_t r = queryStartPos; r <= queryEndPos ; r++ )
    {
        if (verbose) std::cout << alnPair[r].probability << std::endl; 
        if (alnPair[r].probability > bestProb) bestProb = alnPair[r].probability;
    }
    if (verbose) std::cout << std::endl;
    assert(bestProb > 0.0);
    double errProb = 1 - bestProb ; 
    if (errProb == 0) errProb = probZero;
    if (verbose) std::cout << "Errprob" << errProb << std::endl;
         
    
    //cigar string . because reverse strand. we have to reverse-iterate
    std::stringstream cigar;
    if (queryEndPos < alnPair.size()-2 )  //-2  because of dummy align at end
    {
        cigar << alnPair.size()-2-queryEndPos << "H";
    }
            
    std::string currentState = preCigar.back();
    int length=0;
    for (auto rev=preCigar.rbegin(); rev != preCigar.rend(); ++rev)
    {
        auto state = *rev;
        if (currentState == state)
        {
            length++;
        }
        else
        {
            cigar << length << currentState;
            length = 1;
            currentState = state;
        }
    }
    cigar << length << currentState;
    if (queryStartPos > 0) cigar << queryStartPos << "H" ;
    
    //querySeq also needs to be reversed. only order, no need for complementary sequence because lastal already reports that
    reverse(querySeq.begin(),querySeq.end());    
    
    tempSAM.cigar = cigar.str();
    tempSAM.errProb = errProb ; 
    tempSAM.isFirst = isFirst;
    tempSAM.isSupplementary = isSupplementary;
    tempSAM.qNameNoPair = qNameNoPair;
    tempSAM.querySeq = querySeq; 
    tempSAM.queryStartPos = alnPair.size()-queryEndPos-2; //pos on - strand  which aligns to refEndpos = tempSAM.refStartPos. -2 and not -1 because we added dummy element to alnPair.
    tempSAM.queryStrand = queryStrand;
    tempSAM.refName = refName;
    tempSAM.refStartPos = refEndPos+1; //SAM is 1-based 
    tempSAM.tlen = 0; 
            
    samStructVec.push_back(tempSAM);
}


/*
void buildSAMsegmentMinus()
{
    if (verbose) std::cout << "constructing SamAln struct for positive strand" << std::endl;
            
    //cigar
    std::stringstream cigar;
    if (i < alnPair.size()-1 )  //remember dummy align at end?
    {
        cigar << alnPair.size()-i-1 << "H";
    }
            
    std::string currentState = preCigar.back();
    int length=0;
    for (auto rev=preCigar.rbegin(); rev != preCigar.rend(); ++rev)
    {
        auto state = *rev;
        if (currentState == state)
        {
            length++;
        }
        else
        {
            cigar << length << currentState;
            length = 1;
            currentState = state;
        }
    }
    cigar << length << currentState;
    if (queryStartPos > 0) cigar << queryStartPos << "H" ;
    if (verbose) std::cout << "constructed cigar" << std::endl;
    
    //reverse querySeq.  maybe, not do it in-place?
    reverse(querySeq.begin(),querySeq.end());
    std::cout << querySeq << std::endl;
    
    //refStartPos   needs to start with what aligns with pos i, not queryStartPos
    tempSAM.cigar = cigar.str();
    tempSAM.errProb = errProb ; 
    tempSAM.isFirst = isFirst;
    tempSAM.isSupplementary = isSupplementary;
    tempSAM.qNameNoPair = qNameNoPair;
    tempSAM.querySeq =  querySeq;
    tempSAM.queryStrand = queryStrand;
    tempSAM.refName = refName;
    tempSAM.refStartPos = refEndPos+1; //SAM is 1-based 
    tempSAM.tlen = 0;


    samStructVec.push_back(tempSAM);
}

*/
void buildSAM(const std::string &qNameNoPair, std::vector<AlignmentPair> &alnPair, // no const as we add a dummy aln
               const AlignmentParameters &params,  std::vector<SAMaln> &samStructVec, const bool &isFirst, bool verbose)
{
    
    if (verbose) std::cout << "inside buildSAM" << std::endl;
    if (alnPair.empty()) return ;
    
    //alnPair[i].index can be one of followng: (1)some actual ref position or (2) NOALIGN or (3)GAP
    //logic in this function flows based on these cases
    long a = params.aGet();
    long b = params.bGet();
    long x = params.xGet();
    
    //since we need to look ahead in alnpair to make decisions, the loop is cleaner with a dummy alignment at the tail end 
    AlignmentPair dummy;
    dummy.rIndex = NOALIGN;
    dummy.probability = 0.0;
    alnPair.push_back(dummy);
    
    bool isSupplementary=false;

    size_t pos = 0; //position along alnpair = position along read
    while(pos < alnPair.size()-1) 
    {
        //start new alignment segment (i.e. alignment corresponding to one SAM line)
        if (verbose) std::cout << "starting new segment at qPos: " << pos << std::endl;
        
        //alignment segments must start with "M"
        if(alnPair[pos].rIndex == NOALIGN || alnPair[pos].rIndex == GAP ) //\label{getToFirstM}
        {
            if(verbose) std::cout << " Encountered GAP or NOALIGN at qPos:" << pos << std::endl;
            pos++ ;
            continue;
        }
        
        size_t j = pos; //paranoid. to check later that we did move forward
        if (alnPair[pos].qStrand == '+') buildSAMsegmentPlus(pos,qNameNoPair, samStructVec, alnPair, a,b,x,isFirst,isSupplementary, verbose);
        else if (alnPair[pos].qStrand == '-') buildSAMsegmentMinus(pos,qNameNoPair, samStructVec, alnPair, a,b,x,isFirst,isSupplementary, verbose);
        else 
        {
           std::cout << "Neither + or - strand. Illegal character" << std::endl;
           std::exit(EXIT_FAILURE);
        }
        assert(pos > j); //paranoid. check that we did move forward
        
        //prepare for next round
        isSupplementary = true; //any alignment segment that might come next is supplementary 
    }
}

void getMateInfo(const std::vector<SAMaln> &readSAMs, const std::vector<SAMaln> &mateSAMs, 
                  std::string &rnext, long &pnext, bool &nextSegmentUnmapped, bool &nextSegmentRevComp)
{
    //default values
    rnext="*"; 
    pnext=0;
    nextSegmentUnmapped = false;
    nextSegmentRevComp = false;
    
    if (mateSAMs.empty()) nextSegmentUnmapped = true;
    else if (mateSAMs.size() == 1) // not sure what to do if mate is split into multiple segments
    {
        if (readSAMs.size() == 1) //for now, insist that read also has single alignment only
        { 
            if (readSAMs[0].refName == mateSAMs[0].refName) rnext = "=";
            else rnext = mateSAMs[0].refName;
            pnext = mateSAMs[0].refStartPos;
        }
        if (mateSAMs[0].queryStrand == '-')
        {
            nextSegmentRevComp = true;
        }
    }
}

void outputSAM(const std::string qNameNoPair, std::vector<AlignmentPair> &read1AlnPair, std::vector<AlignmentPair> &read2AlnPair, const AlignmentParameters &params, LastPairProbsOptions &opts, bool verbose)
               
{   
    std::vector<SAMaln> read1SAMs;
    std::vector<SAMaln> read2SAMs;
    
    long read1Len = read1AlnPair.size();
    long read2Len = read2AlnPair.size();
    
    if (verbose) std::cout << "BuildSAM for read1" << std::endl;
    if (!(read1AlnPair.empty())) buildSAM(qNameNoPair, read1AlnPair, params, read1SAMs, true, verbose);
    if (verbose) std::cout << "BuildSAM for read2" << std::endl;
    if (!(read2AlnPair.empty())) buildSAM(qNameNoPair, read2AlnPair, params, read2SAMs, false, verbose);
    
    // to compute the "proper paired-end" flag. 
    int fragmentLen;
    bool isProperPair = false; //make it member of SAMaln struct?
    if (read1SAMs.size() == 1 && read2SAMs.size() == 1 ) 
    {
        if (read1SAMs[0].refName == read2SAMs[0].refName)
        {
            if (read1SAMs[0].queryStrand == '+' && read2SAMs[0].queryStrand == '-')
            {   
                fragmentLen = read2SAMs[0].refStartPos - read1SAMs[0].refStartPos    +   read1SAMs[0].queryStartPos + (read2Len - read2SAMs[0].queryStartPos);  
                //                  difference between start points                  +     read edges not reported in the alignment
                if ( fragmentLen <= opts.fraglen + 3*opts.sdev && fragmentLen >= opts.fraglen - 3*opts.sdev )
                {
                    isProperPair = true;
                }
                read1SAMs[0].tlen = fragmentLen;
                read2SAMs[0].tlen = -1 * fragmentLen;
            }
            else if (read1SAMs[0].queryStrand == '-' && read2SAMs[0].queryStrand == '+') 
            {
                fragmentLen = read1SAMs[0].refStartPos - read2SAMs[0].refStartPos + read2SAMs[0].queryStartPos + read1Len - read1SAMs[0].queryStartPos ;
                if ( fragmentLen <= opts.fraglen + 3*opts.sdev && fragmentLen >= opts.fraglen - 3*opts.sdev )
                {
                     isProperPair = true;
                 }
                read1SAMs[0].tlen = -1 * fragmentLen;
                read2SAMs[0].tlen = fragmentLen;
            }
        }
    }
    
    //mate info
    std::string rnext; 
    long pnext;
    bool nextSegmentUnmapped;
    bool nextSegmentRevComp;
    if (!(read1AlnPair.empty())) 
    {
        if (verbose) std::cout << "Get mate info for read1" << std::endl;
        getMateInfo(read1SAMs,read2SAMs,rnext,pnext,nextSegmentUnmapped,nextSegmentRevComp);
        if (verbose) std::cout << "Writing output" << std::endl;
        for (auto &samLine: read1SAMs)
        {
            writeSAMoutput(samLine,rnext,pnext,nextSegmentUnmapped,nextSegmentRevComp,isProperPair);
        }
    }
    
    if (!(read2AlnPair.empty())) 
    {
        if (verbose) std::cout << "Get mate info for read2" << std::endl;
        getMateInfo(read2SAMs,read1SAMs,rnext,pnext,nextSegmentUnmapped,nextSegmentRevComp);
        if (verbose) std::cout << "Writing output" << std::endl;
        for (auto &samLine: read2SAMs)
        {
            writeSAMoutput(samLine,rnext,pnext,nextSegmentUnmapped,nextSegmentRevComp,isProperPair);
        }
    }
       
}
               
               
/*
void outputSAMold(const std::string qNameNoPair, const std::vector<AlignmentPair> &alnpair, const std::vector<AlignmentPair> &mateAlnPair, 
               const AlignmentParameters &params, const bool isFirst)
{
   
    
    long a = params.aGet();
    long b = params.bGet();
    long x = params.xGet();
    
    std::string seq;
    size_t idx = 1;
    bool isSupplementary = false;
    while(idx < alnpair.size()) 
    {
        bool isEndAlignment = false;
        if(idx == alnpair.size()-1) isEndAlignment = true;
        // if match, start one alignment
        if(alnpair[idx-1].rName == alnpair[idx].rName &&
           alnpair[idx-1].rIndex == alnpair[idx].rIndex - 1 &&
           alnpair[idx-1].rStrand == alnpair[idx].rStrand) // else is weak
        {
            
            size_t totalLength = 0;
            std::stringstream cigar;
            if(idx > 1) 
            {
                cigar << idx-1 << "H";
                totalLength += idx - 1;
            }
            // count match length
            size_t matchLength = 2;
            seq = "";
            seq += alnpair[idx-1].qBase;
            seq += alnpair[idx].qBase;
            for(size_t i=idx+1; ; ++i) 
            {
                if(alnpair[i-1].rName == alnpair[i].rName &&
                   alnpair[i-1].rStrand == alnpair[i].rStrand) 
                {
                    
                    auto refIndexDiff = alnpair[i].rIndex-alnpair[i-1].rIndex;
                    if(refIndexDiff == 1) 
                    {
                        matchLength++;
                        seq += alnpair[i].qBase;
                    } 
                    else 
                    {
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
                    
                } 
                else 
                {
                    bool isInsertion = false;
                    for(size_t j=i+1; a+b*(j-i)<x; ++j) {
                        auto refIndexDiff = alnpair[i].rIndex-alnpair[i-1].rIndex;
                        if(alnpair[j-1].rName == alnpair[j].rName &&
                           alnpair[j-1].rStrand == alnpair[j].rStrand &&
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
                if(i == alnpair.size()-1) 
                {
                    isEndAlignment = true;
                    cigar << matchLength << "M";
                    totalLength += matchLength;
                }
                if(isEndAlignment) 
                {
                    outputInSAM(qNameNoPair,alnpair[idx].rName, alnpair[idx-1].rIndex, alnpair[idx].qStrand, alnpair.size() ,totalLength, seq, cigar, isSupplementary, isFirst, mateAlnPair);
                    idx = i + 1;
                    isSupplementary = true;
                    break;
                }
            }
        } 
        else 
        {
            idx++;
        }
    }
} 
*/
   
void outputNative(std::vector<AlignmentPair> readAln)
{
    // invalid input
    if(readAln.size() == 0) {
        return;
    }
    std::cout << readAln[0].qName;
    for (auto &p : readAln)
    {
        std::cout << "\t(" << p.rName << "," << p.rIndex << "," << p.qStrand << "," << p.probability << ")";
    }
    std::cout << std::endl;
}

void startSplitPEProcess(const std::string readNameNoPair, std::vector<Alignment> &alns1, std::vector<Alignment> &alns2, AlignmentParameters &params, LastPairProbsOptions &opts)
{
    
    bool verbose=false;
    
    if (verbose) std::cout << "Processing read: "<< readNameNoPair << std::endl;
    std::vector<AlignmentPair> read1FinalAln; 
    std::vector<AlignmentPair> read2FinalAln; 
        
    if (!alns1.empty())
    {
        if (verbose) std::cout << "start calcProb" << std::endl; 
        std::vector<std::vector<AlignmentPair>> read1Probs = calcProb(alns1, alns2, params, opts,verbose);
        /*
        for(auto& pos:read1Probs){
            for(auto& pair:pos){
                std::cout << pair.probability << ":" << pair.rIndex << ";\t" ; 
            }
            std::cout << std::endl;
        }
        */
        chooseBestPair(read1Probs,read1FinalAln);
        if (!opts.isSamFormat) outputNative(read1FinalAln);
    }
    
    if (!alns2.empty())
    {
        if (verbose)  std::cout << "start calcProb" << std::endl; 
        std::vector<std::vector<AlignmentPair>> read2Probs = calcProb(alns2, alns1, params, opts, verbose);
        chooseBestPair(read2Probs,read2FinalAln);
        if (!opts.isSamFormat) outputNative(read2FinalAln);
    }
    
    if (opts.isSamFormat)
    {
         
        if (verbose)  std::cout << "output SAM" << std::endl ; 
        outputSAM(readNameNoPair,read1FinalAln, read2FinalAln,params, opts, verbose);
    }
    /*
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
    */
   
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
                currentAln = nextAln;
                continue;
            }
            else
            {   //different pair
                startSplitPEProcess(currentAln.qNameNoPair,X, Y, params, opts);
                currentAln = nextAln;
                X.clear();
                Y.clear();
                continue;
            }
        }
        else
        {   // readSingleAlignment didn't read in alignment because of EOF or other reasons.
            startSplitPEProcess(currentAln.qNameNoPair,X, Y, params, opts);
            break;
        }
    }
}
