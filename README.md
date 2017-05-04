# Comparison of pair matching algorithms

This is the source code for a paper published in JASSS by Nathan Geffen and
Stefan Scholz titled *Efficient and effective pair-matching algorithms for
microsimulations*.

To run the tests described in that paper, execute run_tests.sh. The results are
written to table.tex. Some graphs are also produced.

The C++ implementation of the pair matching algorithms is in partnermatch.cc.
To compile it for the STIMOD microsimulation, run *make release*. To compile it for
the ATTRACTREJECT microsimulation, run *make release_attract*.

## Command line options for the C++ implementation

-n POSITIVE_INTEGER : size of population

-k POSITIVE_INTEGER : number of neighbors for CSPM, RKPM, and WSPM

-c POSITIVE_INTEGER : number of clusters for CSPM

-r POSITIVE_INTEGER : number of simulations

-i POSITIVE_INTEGER : number of iterations for single simulation

-s POSITIVE_INTEGER : seed for random number generator

-f ZERO_OR_POSITIVE_INTEGER : number identifying first simulation

-a [RNWDCB] : which algorithms to run (R=RPM, N=RKPM, W=WSPM, D=DCPM, C=CSPM,
 B=BFPM)

-A 0..1 : Attractor value for ATTRACTREJECT simulation

-R 0..1 : Rejector value for ATTRACTREJECT simulation

-o filename : file to write partnership info to (hardly ever needed)

-d [A..Za..z1..9_]+ : Descriptive identifier for simulation

-b : run Blossom V algorithm

-t : only do timings

-h : write CSV header

-v : verbose output

-y POSITIVE_INTEGER : ages to go up and down for distribution matching (leave
for default of 4 which is sensible)
