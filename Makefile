CC=g++
CFLAGS=-fopenmp -Wall -O3 -std=c++11 -g -fno-omit-frame-pointer -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -Wl,--no-as-needed -lprofiler -ltcmalloc -Wl,--as-needed #-parallel -xHost #-prof-gen -prof-dir=./  -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=1 

BASESRC=src/alphabet.cpp src/sequencefamily.cpp src/spacedworddb.cpp src/rswset.cpp src/representspw.cpp src/familyscore.cpp src/sequence.cpp src/patternset.cpp src/spacedword.cpp src/pattern.cpp
BASEHDR=src/alphabet.hpp src/sequencefamily.hpp src/spacedworddb.hpp src/rswset.hpp src/representspw.hpp src/familyscore.hpp src/sequence.hpp src/patternset.hpp src/spacedword.hpp src/pattern.hpp
BASEOBJ=$(BASESRC:.cpp=.o)

SWDBSRC=src/swdb.cpp src/dboptions.cpp src/createdatabase.cpp src/familyspacedword.cpp
SWDBHDR=src/dboptions.hpp src/createdatabase.hpp src/familyspacedword.hpp
SWDBOBJ=$(SWDBSRC:.cpp=.o)
SWDBPRG=swdb

SWDSSRC=src/swds.cpp src/dsoptions.cpp src/searchdatabase.cpp src/spacedhit.cpp
SWDSHDR=src/dsoptions.hpp src/searchdatabase.hpp src/spacedhit.hpp
SWDSOBJ=$(SWDSSRC:.cpp=.o)
SWDSPRG=swds

all: swdatabase swdatasearch

swdatabase: $(SWDBSRC) $(BASESRC) $(SWDBPRG)
swdatasearch: $(SWDSSRC) $(BASESRC) $(SWDSPRG)

$(SWDBPRG): $(SWDBOBJ) $(BASEOBJ)
	$(CC) $(CFLAGS) $(SWDBOBJ) $(BASEOBJ) -o $@
$(SWDSPRG): $(SWDSOBJ) $(BASEOBJ)
	$(CC) $(CFLAGS) $(SWDSOBJ) $(BASEOBJ) -o $@

.cpp.o: $(SWDBHDR) $(SWDSHDR) $(BASEHDR)
		$(CC) -c $(CFLAGS) $< -o $@

clean:
	find ./src/ -name "*.o" -delete
	find ./ -name $(SWDBPRG) -delete
	find ./ -name $(SWDSPRG) -delete
