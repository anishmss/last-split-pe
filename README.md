# last-split-pe
last-split-pe is a method that can split-align a short DNA read to a reference genome. It achieves high accuracy by combining probabilistic alignments with information from paired-end reads.

Here's the basic idea:

![Alt text](last-split-pe-graphics.png)

# Download and Installation
You can download a copy of the source code from the `Downloads` link in the navigation bar, or you can clone the repository:  
~~~~
git clone https://bitbucket.org/splitpairedend/last-split-pe.git
~~~~

Once inside the `src` folder, run `make`.
You might want to place the `last-split-pe` binary in a standard bin directory of your OS.

You will also need to download and install [LAST](https://gitlab.com/mcfrith/last). Make sure you are able to run the following LAST commands: `last-pair-probs`, `lastdb`, `lastal`, `last-split`, `fastq-interleave`.

# Usage
The worflow consists of first using LAST to compute high-scoring local alignments of each read, ignoring pairing. Then `last-split-pe` computes a final (split-)alignment utilizing the pairing information.

Suppose we have a reference genome in fasta format in the file `refGenome.fa`, and paired-end reads in two files `reads1.fastq` and `reads2.fastq`.

Here are the commands:
~~~~
lastdb -uNEAR db refGenome.fa
fastq-interleave reads1.fastq reads2.fastq |
lastal -Q1 db   | 
last-split -m 0.9 -n -fMAF+ | last-split-pe -f MEAN -s STD_DEV > output.sam
~~~~
MEAN and STD_DEV stand for the mean and standard deviation of the distribution of fragment (from which the read pairs are sequenced) lengths. This can be estimated from a sample of the reads by using LAST's `last-pair-prob` as follows:

~~~~
fastq-interleave reads1.fastq reads2.fastq | head -n80000 > sample.fastq
lastal -Q1 -D1000 -i1 db sample.fastq | last-pair-probs -e
~~~~

## What are all the commands doing?
`lastdb` constructs a suffix array of the reference. The -uNEAR option specifies a seeding strategy to find short and strong similarities (e.g. when comparing human reads to human reference). 

`fastq-interleave` merges the two fastq files into one, such that reads of the same pair appear consecutively. It also appends "/1" and "/2" to the ends of the names of pair members, if that's not already the case.


`lastal` finds high-scoring local alignments of each read to the reference. 

`last-split` computes the probability of each column of each alignment. The -m option discards local alignments with error probability more than 0.9. Decrease it to get more confident local alignments. The -n option, which is necessary, tells `last-split` not to compute a final split-alignment, as this will be handled by the next command. 

`last-split-pe` updates the column probabilities based on alignments of the paired-end mates and reports a (split-)alignment of the reads.   By default, it reports alignments of error probability less than 0.01; this can be controlled using the `-m` option. 
 
## Where is the SAM header?
The `printSamHeader.py` script (inside the `scripts` folder) reads in the reference genome fasta file and generates SAM header lines containing information about the reference. This can be appended to the top of the SAM output. 
~~~~
python printSamHeader.py refGenome.fa > samHeader
cat samHeader output.sam > output_withHeader.sam
~~~~

The following achieves the same result without having to hold output.sam separately on disk.
~~~~
{
python printSamHeader.py refGenome.fa;
lastal -Q1 -i1 db reads.fastq | last-split -m 0.9 -n -fMAF+ | last-split-pe -f MEAN -s STD_DEV ;
} > output_withHeader.sam
~~~~

## Multi-threaded operation
Make sure  [GNU parallel](https://www.gnu.org/software/parallel) is installed. 

Assume we have already constructed the index of the reference and that reads have been interleaved.  

Then, suppose we want to launch 10 threads.
~~~~
parallel --gnu --pipe -L8 -j10 "lastal -Q1 -i1 db  | last-split -m 0.9 -n -fMAF+ | last-split-pe -f MEAN -s STD_DEV " < interleaved.fastq > output.sam
~~~~

## If reads are in fasta format
We need to modify only the `lastal` arguments as follows:
~~~~
lastal  -Q0 -i1 db reads.fasta 
~~~~
By default, LAST uses a different scoring if the query is in fasta format.
See LAST options for how to set match/mismatch scores and gap penalties.
See [last-train](https://gitlab.com/mcfrith/last/-/blob/main/doc/last-train.rst) to automatically find suitable scoring parameters for you data.

## Publication and Supplementary Data 
Details of the algorithm can be found in the following paper:

Anish M S Shrestha, Naruki Yoshikawa, and  Kiyoshi Asai :
Combining probabilistic alignments with read pair information improves accuracy of split-alignments, Bioinformatics, Volume 34, Issue 21, 1 November 2018, Pages 3631–3637

Scripts and data used for performance evaluation are available [here](https://drive.google.com/file/d/1fC83huEt2lXHo-BnFEnu0Pqcg1unBviZ/view?usp=sharing) .

## Contact
anish.shrestha --atmark-- dlsu.edu.ph or
anish.ms.shrestha --atmark-- gmail.com
