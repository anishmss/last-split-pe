CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -pedantic -g -std=c++11 -pthread -DHAS_CXX_THREADS -Wno-sign-compare -Wfatal-errors

SPOBJ = last-split-pe.o last-split-pe-main.o io.o

all: last-split-pe
last-split-pe: $(SPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SPOBJ)
#last-pair-probs: $(PPOBJ)
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(PPOBJ)
#last-pair-probs.o: last-pair-probs.cc last-pair-probs.hh io.hh stringify.hh
#last-pair-probs-main.o: last-pair-probs-main.cc last-pair-probs.hh stringify.hh version.hh
last-split-pe.o: last-split-pe.cc last-split-pe.hh io.hh stringify.hh
last-split-pe-main.o: last-split-pe-main.cc last-split-pe.hh stringify.hh version.hh
io.o: io.cc io.hh

clean: 
	rm -f last-split-pe *.o
