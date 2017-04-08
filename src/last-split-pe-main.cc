// Copyright 2017 Naruki Yoshikawa

// Parse the command line and run last-pair-probs.

#include "last-split-pe.hh"
#include "stringify.hh"

#include <getopt.h>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>
#include <new>  // bad_alloc

static void run(int argc, char* argv[]) {
  LastPairProbsOptions opts;

  opts.rna = false;
  opts.estdist = false;
  opts.mismap = 0.01;
  opts.isFraglen = false;
  opts.isSdev = false;
  opts.isDisjoint = false;

  const char *version = "last-split-pe "
#include "version.hh"
"\n";

  std::string help = "\
Usage:\n\
  " + std::string(argv[0]) + " --help\n\
\n\
Options:\n\
  -h, --help            show this help message and exit\n\
  -V, --version         show program's version number and exit\n\
";

  const char sOpts[] = "hrem:f:s:d:c:V";

  static struct option lOpts[] = {
    { "help",     no_argument,       0, 'h' },
    { "version",  no_argument,       0, 'V' },
    { 0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return;
    case 'V':
      std::cout << version;
      return;
    case '?':
      throw std::runtime_error("");
    }
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
