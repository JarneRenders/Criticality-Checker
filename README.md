# CriticalityChecker
This repository contains a filter program for checking properties involving the criticality of snarks. A snark G is *critical* if G - {u,v} is 3-edge-colourable for every edge uv; *cocritical* if G - {u,v} is 3-edge-colourable for every non-adjacent pair u,v and *bicritical* if G - {u,v} is 3-edge-colourable for any pair u,v.

The latest version of this program can be obtained from <https://github.com/JarneRenders/Criticality-Checker>.

This program can be used the verify whether a snark is critical, cocritical or bicritical. These properties can also be verified for arbitrary cubic graphs. The program can also be used to determine 3-edge-colourability of the input graphs.

### Installation

This requires a working shell and `make`. On Windows an easy way to simulate this is by using Windows Subsystem for Linux (WSL).

- Compile using: 
  * `make` to create a binary for the 64 bit version
  * `make 128bit` to create a binary for the 128 bit version
  * `make 128bitarray` to create a binary for the 128 bit version which implements bitsets using arrays
  * `make all` to create all of the above

The 64 bit version is much faster than the 128 bit version, hence it is recommended to use this one when running the program on graphs with 64 vertices or fewer.
There is a difference in implementation between both of the 128 bit versions, one might be faster than the other. Both can handle input graphs up to order 128.
Use `make clean` to remove all binaries created using `make`.

### Usage of CriticalityChecker

This helptext can be found by executing `./criticalityChecker -h`.

Usage: `./criticalityChecker [b|c|3] [v] [h]`

Filter cubic graphs satisfying certain criticality requirements. Can also be used to determine 3-edge-colourability of cubic graphs.

Graph are read from stdin in graph6 format. Graphs are sent to stdout in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>.

Without any arguments this program output graphs which are critical.
```
    -3, --colourability
            outputs graphs which are 3-edge-colourable; does not work
            with -b or -c
    -b, --bicritical
            outputs graphs which are bicritical; does not work with
            -3 or -c
    -c, --cocritical
            outputs graphs which are cocritical; does not work with
            -3 or -b
    -h, --help
            outputs this helptext
    -v, --verbose
            sends extra information to stderr
```
