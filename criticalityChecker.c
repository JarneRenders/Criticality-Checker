/**
 * criticalityChecker.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE "Usage: ./criticalityChecker [b|c|3] [v] [h]"
#define HELPTEXT \
"Filter cubic graphs satisfying certain criticality requirements.\n\
Can also be used to determine 3-edge-colourability of cubic graphs.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to stdout\n\
in graph6 format.\n\
\n\
Without any arguments this program output graphs which are critical.\n\
\n\
    -3, --colourability\n\
            outputs graphs which are 3-edge-colourable; does not work\n\
            with -b or -c\n\
    -b, --bicritical\n\
            outputs graphs which are bicritical; does not work with\n\
            -3 or -c\n\
    -c, --cocritical\n\
            outputs graphs which are cocritical; does not work with\n\
            -3 or -b\n\
    -h, --help\n\
            outputs this helptext\n\
    -v, --verbose\n\
            sends extra information to stderr\n\
"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include "utilities/readGraph6.h"
#include "utilities/bitset.h"

struct graph {
    int nv;
    bitset *adjacencyList;

    // 0 will be uncoloured, 1,2,3 will be the colours
    int colouring[MAXBITSETSIZE*MAXBITSETSIZE];

    bitset availableColours[MAXBITSETSIZE];
    bitset allColMask;

};

struct options {
    bool haveModResPair;
    bool bicriticalFlag;
    bool cocriticalFlag;
    bool colourabilityFlag;
    bool verboseFlag;
    int res;
    int mod;
};

struct counters {
    long long unsigned int skippedGraphs;
    long long unsigned int callsToIsColourable;
    long long unsigned int removedLeaf;
    long long unsigned int recursionSteps;
    long long unsigned int propagationSteps;
    long long unsigned int minRecursionSteps;
    long long unsigned int maxRecursionSteps;
};

void printBitset(bitset set) {
    forEach(el, set) {
        fprintf(stderr, "%d ", el);
    }
    fprintf(stderr, "\n");
}

void printGraph(bitset adjacencyList[], int nv) {
    for(int i = 0; i < nv; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

int readGraph(const char *graphString, struct graph *g,
 struct options *options, struct counters *counters) {
    g->nv = getNumberOfVertices(graphString);
    if(g->nv == -1 || g->nv > MAXBITSETSIZE) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        counters->skippedGraphs++;
        return 1;
    }
    g->adjacencyList = malloc(sizeof(bitset)*g->nv);
    if(loadGraph(graphString, g->nv, g->adjacencyList) == -1) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        counters->skippedGraphs++;
        return 1;
    }
    return 0;
}

void freeGraph(struct graph *g) {
    free(g->adjacencyList);
}

//**********************************************************************
//
//          
//
//**********************************************************************


bool canColour(struct graph *g, int v, int col) {
    return contains(g->availableColours[v], col);
}

#define GETCOLOUR(g, u, v) g->colouring[u * MAXBITSETSIZE + v]
#define SETCOLOUR(g, u, v, col) g->colouring[u * MAXBITSETSIZE + v] =\
                                 col;

// Colour edge uv with col. Only works properly if edge can be coloured
// in this colour. 
void colour(struct graph *g, int u, int v, int col) {
    SETCOLOUR(g, u, v, col);
    SETCOLOUR(g, v, u, col);
    removeElement(g->availableColours[u], col);
    removeElement(g->availableColours[v], col);
}

// Uncolour the edge uv. Only works properly if it was coloured before.
void uncolour(struct graph *g, int u, int v) {
    int col = GETCOLOUR(g, u, v);
    SETCOLOUR(g, u, v, 0);
    SETCOLOUR(g, v, u, 0);
    add(g->availableColours[u], col);
    add(g->availableColours[v], col);
}

bool colouringRecursion(struct graph *g, bitset removedVertices, 
 bitset toCheck, bitset checked, struct counters *counters);

// When colouring an edge, continue colouring if there are edges that
// should receive a fixed colour. Returns true if successful and false
// if this leads to a contradiction.
bool propagate(struct graph *g, bitset removedVertices,
 bitset propagateVertices, bitset toCheck, bitset checked,
 struct counters *counters) {

    counters->propagationSteps++;

    // Stop propagating if no vertices left
    if(isEmpty(propagateVertices)) { // Continue normal recursion
        if(colouringRecursion(g, removedVertices, toCheck,
         checked, counters)) {
            return true;
        }
        return false;
    }

    int v = next(propagateVertices, -1);
    removeElement(propagateVertices, v); // If bottleneck, do better

    if(size(g->availableColours[v]) != 1) {
        if(propagate(g, removedVertices, propagateVertices, toCheck,
         checked, counters)) {
            return true;
        }
        return false;
    }

    // Now there is one colour left.

    // Check if degree 2 vertex; They always have one colour left, but
    // no options
    bitset presentNbrs = 
     difference(g->adjacencyList[v], removedVertices);
    if(size(presentNbrs) == 2) {
        if(propagate(g, removedVertices, propagateVertices, toCheck,
         checked, counters)) {
            return true;
        }
        return false;
    }

    // Only one colour left for v. Check if only one uncoloured edge
    // left. Redundant check, but we need to find the specific edge.
    bitset unColNbrs = difference(presentNbrs, checked);
    if(size(unColNbrs) > 1) {
        fprintf(stderr, "Error: this should not happen. (Degree 3"
         "vertex with one colour left can only have one edge left.\n");
        exit(1);
    }

    // Only one uncoloured edge left. Let w be its other endpoint and
    // color it if possible.
    int w = next(unColNbrs, -1);
    int col = next(g->availableColours[v], -1);
    if(!canColour(g, w, col)) return false;
    colour(g, v, w, col);

    add(toCheck, w); 
    removeElement(toCheck, v);
    add(checked, v);
    add(propagateVertices, w);

    if(propagate(g, removedVertices, propagateVertices, toCheck,
     checked, counters)) {
        uncolour(g, v, w);
        return true;
    }

    uncolour(g, v, w);

    return false;
}

bool colouringRecursion(struct graph *g, bitset removedVertices, 
 bitset toCheck, bitset checked, struct counters *counters) {

    counters->recursionSteps++;

    // We now have a 3-colouring.
    if(isEmpty(toCheck)) { 

        // Exit the program if something went wrong (should not happen).
        // Does not impact computation time that much.
        for(int i = 0; i < g->nv; i++) {
            if(contains(removedVertices, i)) continue;
            bitset cols = EMPTY;
            forEach(el, g->adjacencyList[i]) {
                if(contains(removedVertices, el)) continue;
                int col = GETCOLOUR(g, i, el);
                if(col < 1 || col > 3) {
                    fprintf(stderr,
                     "Error: colour not filled in or incorrect\n");
                    exit(1);
                }
                if(contains(cols, col)) {
                    fprintf(stderr,
                     "Error: multiple same cols at vertex\n");
                    exit(1);
                }
                add(cols, col);
            }
        }

        return true;
    }

    // Colour the remaining edges incident with u in all possible ways
    int u = next(toCheck, -1);
    bitset availCols = g->allColMask;
    bitset nbrs = difference(g->adjacencyList[u], removedVertices);
    if(size(nbrs) == 1) exit(1);
    bitset edgesToColor = EMPTY;
    forEach(nbr, nbrs) {
        int col = GETCOLOUR(g, u, nbr);
        if(col == 0) {
            add(edgesToColor, nbr);
        }
        else {
            removeElement(availCols, col);
        }
    }

    if(size(edgesToColor) == 3) {
        fprintf(stderr, "Error: this should not happen.\n");
        exit(1);
    }

    // No uncoloured nbrs
    if(isEmpty(edgesToColor)) {
        if(colouringRecursion(g, removedVertices,
         difference(toCheck, singleton(u)),
         union(checked, singleton(u)), counters)) {
            return true;
        }
        return false;
    }

    // Should propagate all vertices having only one neighbour
    int nbr1 = next(edgesToColor, -1);        
    int nbr2 = next(edgesToColor, nbr1);        
    int col1 = next(availCols, -1);
    int col2 = next(availCols, col1);

    // Only one uncoloured edge
    if(nbr2 == -1) {
        if(canColour(g, nbr1, col1)) {
                
            colour(g, u, nbr1, col1);

            bitset newToCheck = toCheck;
            add(newToCheck, nbr1); 
            removeElement(newToCheck, u);
            bitset newChecked = union(checked, singleton(u));

            if(propagate(g, removedVertices, singleton(nbr1),
             newToCheck, newChecked, counters)) {
                uncolour(g, u, nbr1); 
                return true;
            }
            uncolour(g, u, nbr1); 
        }

        // Check if second colour is available
        if(col2 == -1) {
            return false;
        }

        if(canColour(g, nbr1, col2)) {

            colour(g, u, nbr1, col2);

            bitset newToCheck = toCheck;
            add(newToCheck, nbr1); 
            removeElement(newToCheck, u);
            bitset newChecked = union(checked, singleton(u));

            if(propagate(g, removedVertices, singleton(nbr1),
             newToCheck, newChecked, counters)) {
                uncolour(g, u, nbr1); 
                return true;
            }
            uncolour(g, u, nbr1);
        }

        return false;
    }

    // Two uncoloured edges

    if(canColour(g, nbr1, col1) && 
     canColour(g, nbr2, col2)) {
        colour(g, u, nbr1, col1);
        colour(g, u, nbr2, col2);

        bitset newToCheck = toCheck;
        add(newToCheck, nbr1);
        add(newToCheck, nbr2);
        removeElement(newToCheck, u);
        bitset newChecked = union(checked, singleton(u));

        if(propagate(g, removedVertices, union(singleton(nbr1),
         singleton(nbr2)), newToCheck, newChecked, counters)) {
            uncolour(g, u, nbr1);
            uncolour(g, u, nbr2);
            return true;
        }
        uncolour(g, u, nbr1);
        uncolour(g, u, nbr2);
    }

    if(canColour(g,nbr1, col2) && 
     canColour(g, nbr2, col1)) {
        colour(g, u, nbr1, col2);
        colour(g, u, nbr2, col1);

        bitset newToCheck = toCheck;
        add(newToCheck, nbr1);
        add(newToCheck, nbr2);
        removeElement(newToCheck, u);
        bitset newChecked = union(checked, singleton(u));

        if(propagate(g, removedVertices, union(singleton(nbr1),
         singleton(nbr2)), newToCheck, newChecked, counters)) {
            uncolour(g, u, nbr1);
            uncolour(g, u, nbr2);
            return true;
        }
        uncolour(g, u, nbr1);
        uncolour(g, u, nbr2);
    }

    return false;
}

// Determine if cubic graph after removing vertices is 3-edge-colourable
bool isColourable(struct graph *g, bitset removedVertices,
 struct counters *counters) {

    counters->callsToIsColourable++;

    bitset presentVertices = complement(removedVertices, g->nv);
    int u = next(presentVertices, -1);
    while(!isEmpty(intersection(g->adjacencyList[u], removedVertices))) {
        u = next(presentVertices, u);
    } 


    // Colour edges incident with u 
    int col = 1;
    forEach(nbr, g->adjacencyList[u]) {
        colour(g, u, nbr, col);
        col++;
    }

    bitset toCheck = g->adjacencyList[u];
    bitset checked = singleton(u);

    bool foundCol = colouringRecursion(g, removedVertices, toCheck,
     checked, counters);


    // Reset colouring to all zeroes
    forEach(nbr, g->adjacencyList[u]) {
        uncolour(g, u, nbr);
    }

    return foundCol;

}

void initAvailColors(struct graph *g) {
    for(int i = 0; i < g->nv; i++) {
        g->availableColours[i] = g->allColMask;
    }
}

void updateMinAndMax(struct counters *counters,
 long long unsigned int oldSteps) {

    long long unsigned int diff = counters->recursionSteps - oldSteps;
    if(diff > counters->maxRecursionSteps) {
        counters->maxRecursionSteps = diff;
    }
    if(diff < counters->minRecursionSteps) {
        counters->minRecursionSteps = diff;
    }
}

bool isCritical(struct graph *g, struct options *options,
 struct counters *counters) {

    initAvailColors(g);

    for(int u = 0; u < g->nv; u++) {
        forEachAfterIndex(v, g->adjacencyList[u], u) {
            long long unsigned int steps = counters->recursionSteps;
            bool isCol = isColourable(g,
             union(singleton(u), singleton(v)), counters);
            updateMinAndMax(counters, steps);
            if(!isCol) {
                return false;
            }
        }
    }

    return true;
}

// Assumes that g is critical.
bool isCocritical(struct graph *g, struct options *options,
 struct counters *counters) {

    initAvailColors(g);

    for(int u = 0; u < g->nv; u++) {
        for(int v = u+1; v < g->nv; v++) {
            if(contains(g->adjacencyList[u], v)) continue;

            // Vertices of degree 1 have no impact on colourability
            bitset removed =  union(singleton(u), singleton(v));
            for(int i = 0; i < g->nv; i++) {
                if(size(difference(g->adjacencyList[i], removed)) <= 1) {
                    add(removed, i);
                    counters->removedLeaf++;
                }
            }
            long long unsigned int steps = counters->recursionSteps;
            bool isCol = isColourable(g, removed, counters);
            updateMinAndMax(counters, steps);
            if(!isCol) { 
                return false;
            }
        }
    }
    return true;
}

void printExtraOutput(struct counters *counters,
 long long unsigned int nGraphs) {
    fprintf(stderr, "Calls to isColourable: \t %8llu (%.2f per graph)\n",
     counters->callsToIsColourable,
     1.00 * counters->callsToIsColourable / nGraphs);

    fprintf(stderr, "Removed leaves: \t %8llu (%.2f per subgraph)\n",
     counters->removedLeaf,
     1.00 * counters->removedLeaf / counters->callsToIsColourable);

    fprintf(stderr,
     "Recursion steps: \t %8llu (%.2f per subgraph), min: %8llu, max: %8llu\n",
     counters->recursionSteps,
     1.00 * counters->recursionSteps / counters->callsToIsColourable,
     counters->minRecursionSteps,
     counters->maxRecursionSteps);

    fprintf(stderr, "Propagation steps: \t %8llu (%.2f per subgraph)\n",
     counters->propagationSteps,
     1.00 * counters->propagationSteps/ counters->callsToIsColourable);

    long long unsigned int sum = counters->propagationSteps +
     counters-> recursionSteps;
    fprintf(stderr,
     "Sum of steps: \t\t %8llu (%.2f per subgraph) %.2f%% is propagation\n",
     sum,
     1.00 * sum / counters->callsToIsColourable,
     100.00 * counters->propagationSteps / sum);
}

int main(int argc, char ** argv) {
    struct counters counters = {0};
    counters.minRecursionSteps = INT_MAX;
    struct options options = {0};
    int opt;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"colourability", no_argument, NULL, '3'},
            {"bicritical", no_argument, NULL, 'b'},
            {"cocritical", no_argument, NULL, 'c'},
            {"help", no_argument, NULL, 'h'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "3bchv", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'b':
                options.bicriticalFlag = true;
                break;
            case '3':
                options.colourabilityFlag = true;
                break;
            case 'c':
                options.cocriticalFlag = true;
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'v':
                options.verboseFlag = true;
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./criticalityChecker --help for more detailed"
                 " instructions.\n");
                return 1;
        }
    }

    //  Loop over non-option arguments.
    while (optind < argc) {
        if(sscanf(argv[optind], "%d/%d", &options.res, &options.mod)
         == 2) {
            if(options.haveModResPair) {
                fprintf(stderr,
                 "Error: You can only add one mod/res pair as an"
                 " argument.\n");
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./hamiltonicityChecker --help for more detailed"
                 " instructions.\n");
                return 1;
            }
            options.haveModResPair = true;
        }
        else {
            fprintf(stderr,
             "Error: Unknown argument: %s\n", argv[optind]);
            fprintf(stderr, "%s\n", USAGE);
            fprintf(stderr,
             "Use ./criticalityChecker --help for more detailed"
             " instructions.\n");
            return 1;
        }
        optind++;
    }


    // Check only one option is present.
    if((options.bicriticalFlag && options.cocriticalFlag) ||
    (options.bicriticalFlag && options.colourabilityFlag) ||
    (options.colourabilityFlag && options.cocriticalFlag)) {
        fprintf(stderr, "Error: these options are not combinable.\n");
        return 1;
    }

    if(!options.bicriticalFlag && !options.colourabilityFlag &&
     !options.cocriticalFlag) {
        fprintf(stderr, "Outputting critical graphs.\n");
    }
    if(options.cocriticalFlag) {
        fprintf(stderr, "Outputting cocritical graphs.\n");
    }
    if(options.bicriticalFlag) {
        fprintf(stderr, "Outputting bicritical graphs.\n");
    }
    if(options.colourabilityFlag) {
        fprintf(stderr, "Outputting 3-edge-colourable graphs.\n");
    }

    if(!options.haveModResPair) {
        options.mod = 1;
        options.res = 0;
    }

    bitset allCols = union(union(singleton(1), singleton(2)),
     singleton(3));

    unsigned long long int counter = 0;
    unsigned long long int total = 0;
    unsigned long long int passedGraphs = 0;
    clock_t start = clock();

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {

        //  If for graph n: n % mod != res, skip the graph.
        if (total++ % options.mod != options.res) {
            continue;
        }
        struct graph g;
        if(readGraph(graphString, &g, &options, &counters) == 1) {
            fprintf(stderr, "Error: problem loading graph.\n");
            exit(1);
        }
        g.allColMask = allCols;

        counter++;

        if(options.colourabilityFlag) {
            initAvailColors(&g);
            if(isColourable(&g, EMPTY, &counters)) {
                passedGraphs++;
                printf("%s", graphString);
            }
            freeGraph(&g);
            continue;
        }

        if(options.cocriticalFlag) {
            if(isCocritical(&g, &options, &counters)) {
                passedGraphs++;
                printf("%s", graphString);
            }
            freeGraph(&g);
            continue;
        }

        if(options.bicriticalFlag) {
            if(isCritical(&g, &options, &counters) &&
             isCocritical(&g, &options, &counters)) {
                passedGraphs++;
                printf("%s", graphString);
            }
            freeGraph(&g);
            continue;
        }

        // CHECK GRAPH PROPERTY HERE
        if(isCritical(&g, &options, &counters)) {
            passedGraphs++;
            printf("%s", graphString);
        }

        freeGraph(&g);
    }
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    free(graphString);

    if(options.verboseFlag)
        printExtraOutput(&counters, counter);

    if(total == counter)
        fprintf(stderr,
         "\rChecked %llu graphs in %f seconds: %llu passed.\n",
         counter, time_spent, passedGraphs);
    else
        fprintf(stderr,
         "\rChecked %llu/%llu graphs in %f seconds: %llu passed.\n",
         counter, total, time_spent, passedGraphs);

    if(counters.skippedGraphs > 0) {
        fprintf(stderr,
         "Warning: %lld graphs were skipped.\n",
         counters.skippedGraphs);
    }

    return 0;
}