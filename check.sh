#!/bin/bash

mkdir -p results

gprbuild -p -P check.gpr || exit

check 1
check 2
check 3
check 4
check 5
check 6
check 7
check 8

# https://www.math.utah.edu/~beebe/software/ndiff/

echo "diff >"

ndiff -relerr 1.0e-12 -quiet results/data-01.txt expected/data-01.txt
ndiff -relerr 1.0e-12 -quiet results/data-02.txt expected/data-02.txt
ndiff -relerr 1.0e-12 -quiet results/data-03.txt expected/data-03.txt
ndiff -relerr 1.0e-12 -quiet results/data-04.txt expected/data-04.txt
ndiff -relerr 1.0e-12 -quiet results/data-05.txt expected/data-05.txt
ndiff -relerr 1.0e-12 -quiet results/data-06.txt expected/data-06.txt
ndiff -relerr 1.0e-12 -quiet results/data-07.txt expected/data-07.txt
ndiff -relerr 1.0e-12 -quiet results/data-08.txt expected/data-08.txt
