# Regge - a tool to compute defects and their derivatives on a lattice

## Overview

This is an elementary code that computes the defects and their derivatives for a 4-dimensional simplicial lattice. The algorithm is described fully in the paper *"Fast Fast algorithms for computing defects and their derivatives in the Regge calculus"* [Class.Qantum.Grav (28) 2011, p.185005][1]. An earlier version of that paper is included here as paper.pdf.

The main data-structures are two arrays

    lsq    : Array1dReal (1..n_lsq_max)
    simp12 : Array (1..n_simp2_max) of Array1dIntg (1..n_simp12_max);

The first array, lsq, records the squared leg-lengths for every leg in the lattice. The second array, simp12, is an array of arrays, that records the set of legs that define a complete loop around each bone. The first index in simp12 selects the bone while the second index selects the leg. The data in simp12 is stored as follows. Consider a typical bone with vertices (ijk) and suppose that the loop that encloses that bone consists of the the vertices (abcdefa). In this example there are just six distinct vertices (a) through to (f). Each successive pair of vertices along with (ijk) defines a tetrahedron on the bone. The defect is computed from the angles between successive tetrahedra. Suppose (ijk) is the 27th bone in simp12. Then the data in simp12(27) would be

    simp (27)( 1..3)  = (leg_ij,leg_ik,leg_jk)
    simp (27)( 4..7)  = (leg_ia,leg_ja,leg_ka,leg_ab)
    simp (27)( 8..11) = (leg_ib,leg_jb,leg_kb,leg_bc)
    simp (27)(12..15) = (leg_ic,leg_jc,leg_kc,leg_cd)
    simp (27)(16..19) = (leg_id,leg_jd,leg_kd,leg_de)
    simp (27)(20..23) = (leg_ie,leg_je,leg_ke,leg_ef)
    simp (27)(24..27) = (leg_if,leg_jf,leg_kf,leg_fa)
    simp (27)(28..30) = (leg_ia,leg_ja,leg_ka)

The last row in simp(27) is included simply to make the coding a bit easier.

## Compiling

The code is written in Ada and should compile as it stands. A good free compiler for Ada can be found at

   https://www.adacore.com/download/more

The code can be compiled using either

    gnatmake check.adb

or

    gprbuild -p -P check.gpr

## Testing

There are 8 datasets that can be used to test the code (or compare with other codes). Each data set (found in "data/") consists of just one bone surrounded by 6 vertices. For each data set the expected results can be found in "expected/".

The script check.sh compiles and runs the code for each dataset and then uses "ndiff" to compare the actual output (in "results/") against the expected results in (in "expected/").

The [ndiff][2] program is a variation of the standard Unix diff program but with a specific focus on comparing numerical data.

## License

The files regge.adb and regge.ads are distributed under the terms of the [ISC License][3].

 [1]: http://dx.doi.org/10.1088/0264-9381/28/18/185005
 [2]: https://www.math.utah.edu/~beebe/software/ndiff/
 [3]: http://opensource.org/licenses/ISC
