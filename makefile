compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces -O3
profileflags=-std=gnu11 -march=native -Wall -fsanitize=address -g -pg

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: criticalityChecker.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_64_BIT -o criticalityChecker criticalityChecker.c utilities/readGraph6.c $(flags)

# There are two different implementations of the 128-bit version. The array version generally performs faster.
128bit: criticalityChecker.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_128_BIT -o criticalityChecker-128 criticalityChecker.c utilities/readGraph6.c $(flags)

128bitarray: criticalityChecker.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_128_BIT_ARRAY -o criticalityChecker-128a criticalityChecker.c utilities/readGraph6.c $(flags)	

profile: criticalityChecker.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_64_BIT -o criticalityChecker-pr criticalityChecker.c utilities/readGraph6.c $(profileflags)

all: 64bit 128bit 128bitarray

.PHONY: clean
clean:
	rm -f criticalityChecker criticalityChecker-128 criticalityChecker-128a criticalityChecker-pr

