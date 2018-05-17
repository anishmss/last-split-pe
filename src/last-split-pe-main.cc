// Copyright 2017 Naruki Yoshikawa, 
// Copyright 2017 Anish M.S. Shrestha
// Copyright 2014 Toshiyuki Sato
// Copyright 2014 Martin C. Frith

// Parse the command line and run last-split-pec.

#include "last-split-pe.hh"
#include "stringify.hh"

#include <getopt.h>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>
#include <new>  // bad_alloc

static void run(int argc, char* argv[]) {
  LastPairProbsOptions opts;


  opts.mismap = 0.01;
  opts.isFraglen = false;
  opts.isSdev = false;
  opts.isDisjoint = false;
  opts.isNativeFormat = false;

  //const char *version = "last-split-pe "
//#include "version.hh"
//"\n";

  std::string help = "\
Usage:\n\
  " + std::string(argv[0]) + " --help\n\
  " + std::string(argv[0]) + " [options] maf-format-output-of-LAST-SPLIT\n\
\n\
Options:\n\
  -h, --help            show this help message and exit\n\
  -f BP, --fraglen=BP   mean fragment length \n\
  -s BP, --sdev=BP      standard deviation of fragment length\n\
  -m PROB, --mismap=PROB  do not report alignments with mismap (= error probability) > PROB  \n\
                          (default: " + cbrc::stringify(opts.mismap) + ")\n\
  -n , --native      native output format (for testing purposes only)\n\
";

  const char sOpts[] = "hf:s:m:n";

  static struct option lOpts[] = {
    { "help",     no_argument,       0, 'h' },
    { "fraglen",  required_argument, 0, 'f' },
    { "sdev",     required_argument, 0, 's' },
    { "mismap",   required_argument, 0, 'm' },
    { "native",   no_argument,    0, 'n' },
    { 0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return;
    case 'f':
      opts.isFraglen = true;
      cbrc::unstringify(opts.fraglen, optarg);
      break;
    case 'n':
      opts.isNativeFormat = true;
      break;
    case 's':
      opts.isSdev = true;
      cbrc::unstringify(opts.sdev, optarg);
      if (opts.sdev < 0.0) {
        throw std::runtime_error("option -s: should be >= 0");
      }
      break;
    case 'm':
      cbrc::unstringify(opts.mismap, optarg);
      break;
    case '?':
      throw std::runtime_error("");
    }
  }

  if (optind == argc && (!opts.isFraglen || !opts.isSdev)) {
    std::cerr << help;
    throw std::runtime_error("You must specify -f and -s");
  }

  opts.inputFileNames.assign(argv + optind, argv + argc);
  if (opts.inputFileNames.size() > 2) {
    throw std::runtime_error("too many file names");
  }

  std::ios_base::sync_with_stdio(false);  // makes std::cin much faster!!!

  lastSplitPe(opts);
}

int main(int argc, char* argv[]) {
  try {
    run(argc, argv);
    if (!flush(std::cout)) throw std::runtime_error("write error");
    return EXIT_SUCCESS;
  } catch (const std::bad_alloc& e) {  // bad_alloc::what() may be unfriendly
    std::cerr << argv[0] << ": out of memory\n";
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    const char *s = e.what();
    if (*s) std::cerr << argv[0] << ": " << s << '\n';
    return EXIT_FAILURE;
  }
}
