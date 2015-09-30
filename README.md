# PseudoAnosov Mathematica Package

The PseudoAnosov package consists of Mathematica functions for
manipulating characteristic polynomials of pseudo-Anosov maps.

## Main functions

- `DehnTwist[i,{a,b}]` applies the Dehn twist `i` to the curve with homology `{a,b}` on a closed surface of genus `g`.  Here `a` and `b` are lists of length `g` and represent the coefficients in the standard homology basis.  The Dehn twists are the standard "Lickorish generators", numbers from `1` to `3g-1`, with the sign giving the direction of the twist.  `DehnTwist[{i1,i2,...},{a,b}]` applies successive generators starting from the first element of the list.

- `HomologyAction[{i1,i2,...}]` returns the matrix of the action on homology of a sequence of Dehn twists `{i1,i2,...}`.  (See `DehnTwist` for a description of the generators.)  `HomologyAction[{i1,i2,...},g]` specifies the genus `g` explictly, which is otherwise taken as small as possible.  The option `BasisOrder` can be set to `"abab"` or `"aabb"` to specify whether the standard basis for homology should be ordered by hole or by type.

- `LefschetzCombine[L1,L2,...]` adds lists of Lefschetz number.  If they are not the same length, then the blocks are repeated to the length of the longest list. `LefschetzCombine[L1,L2,...,Lm,n]` caps the total length at an integer `n`.

- `LefschetzNumbers[P,k]`, where `P` is the characteristic polynomial of some matrix `M`, returns the Lefschetz number `2-Tr[M^k]`.  `LefschetzNumbers[P,{k2}]` returns a list of Lefschetz numbers `2-Tr[M^k]` for `1 <= k <= k2`.  `LefschetzNumbers[P,{k1,k2}]` returns a list of Lefschetz numbers `2-Tr[M^k]` for `k1 <= k <= k2`.

- `LefschetzNumbersSingularities[k,m,prs]`, where `k` and `m` are integers, lists all possible Lefschetz number sequences for `m` singularities of degree `k` (`m` defaults to `1`).  The sign of the Perron root is given by `prs` (default `-1`).

- `LefschetzNumbersSingularitiesStratum[S,prs]`, where `S` is a list of the degrees of singularities in a stratum, returns a list of all possible Lefschetz number sequences corresponding to the singularities, without taking into account regular orbits.  `S` can be specified as an explicit list (i.e., `{4,2,2,2}`) or in tallied form (`{{4,1},{2,3}}`).  The sign of the Perron root is given by `prs` (default `-1`).

- `LefschetzNumbersTestQ[S,P]` returns `True` if the polynomial `P` is compatible with the stratum `S`.  Possible options are `GiveReasonForRejection` (default `False`), `MaxIterate` (default `1`), and `MaxLefschetz` (default `50`).

- `LefschetzRegularOrbits[L]`, where `L` is a list of Lefschetz numbers, returns the list of regular periodic orbits compatible with `L`, unless an incompatible orbit is detected, in which cases the function stops and returns what it found.  The sign of the Perron root can be specified by the option `PerronRootSign` (default `Automatic`).

- `OrientableStrataList[g]` gives the list of orientable strata for a closed surface of genus `g>1`.  Each stratum in the list is of the form `{k_1,...,k_m}`, where `k_i` is the (even) degree of each singularity, and the sum over the `k_i` gives -2(Euler Characteristic).  Use `Tally/@OrientableStrataList[g]` to group singularities by multiplicity.

- `PolynomialBoundedList[x,n,r,a[n]]` returns a list of polynomials `x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[n-2] x^2 + a[n-1] x + a[n]` with Perron root less than `r`.  For `n` even, only one of each polynomial pair `P(-x)=P(x)` is listed.  If not specified, `a[n]` (determinant) defaults to `1`.

- `PseudoAnosovPerronRootQ[P]` returns `True` if the largest root of the polynomial `P` is nondegenerate (in magnitude) and real.  `PseudoAnosovPerronRootQ[P,xmax]` also returns `False` unless the the Perron root is less than `xmax` (in magnitude).

- `ReciprocalPolynomialBoundedList[x,n,r]` returns a list of reciprocal polynomials `x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1` with Perron root less than `r`.  For `n` even, only one of each polynomial pair `P(-x)=P(x)` is listed.

- `StrataList[g]` gives the list of all strata for a closed surface of genus `g>1`.  Each stratum in the list is of the form `{k_1,...,k_m}`, where `k_i` is the degree of each singularity, and the sum over the `k_i` gives -2(Euler Characteristic).  Use `Tally/@StrataList[g]` to group singularities by multiplicity.

- `StratumDoubleCover[S]` gives the stratum corresponding to the orientating double-cover of the stratum `S={k_1,...,k_m}`.

- `StratumToGenus[S]` gives the genus of the surface containing a stratum `S={k_1,...,k_m}`.

- `StratumOrbits[S,P]` returns a list of possible orbit structure (singular and regular periodic orbits) for the polynomial `P` on stratum `S`.  Returns an empty list if this proves impossible.  `StratumOrbits[S,L]` does the same for a list of Lefschetz numbers `L`.

- `StratumOrbitsTable[so]` presents the output of `StratumOrbits` in a table. `StratumOrbitsTable[so,itmax]` displays at most `itmax` iterates.

## Utility functions

- `GuessPerronRootSign[L]` guesses the Perron root sign by checking for alternating Lefschetz numbers in `L`.

- `IrreducibleMatrixQ[M]` returns `True` if the matrix `M` is irreducible.

- `MahlerMeasure[P]` returns the Mahler measure of the polynomial `P`, which is the absolute value of the product of roots outside the unit circle.

- `PerronRoot[P]` returns the largest root (in magnitude) of the polynomial `P`.

- `PolynomialDegree[P]` returns the degree of the polynomial `P(x)`.

- `PolynomialFromTraces[x,T,det]` creates a polynomial of degree `Length[T]+1` from the determinant `det` (defaults to `1`) and a list `T` of traces of powers of its associated matrix.

- `PolynomialRoots[P]` returns the roots of the polynomial `P`, sorted in decreasing order of magnitude.

- `ReciprocalPolynomial[x,n]` returns a reciprocal polynomial `x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1`.  `ReciprocalPolynomial[x,n,c]` uses `c` as the base name for coefficients.  `ReciprocalPolynomial[x,n,{a_1,...,a_(n/2)}]` uses a list for the coefficient, where `(n/2)` denotes `Floor[n/2]`.

- `ReciprocalPolynomialFromTraces[x,n,T]` creates a reciprocal polynomial of degree `n` from a list of traces of powers of its associated matrix.  `ReciprocalPolynomialFromTraces[x,T]` creates a polynomial of degree `2 Length[T]`.

- `ReciprocalPolynomialQ[P]` returns `True` if `P` is a reciprocal polynomial, i.e. of the form `a[0] x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + a[0]`.

- `TracesPower[P,k]`, where `P` is the characteristic polynomial of some matrix `M`, returns the trace `Tr[M^k]`.  `TracesPower[P,{k2}]` returns a list of traces `Tr[M^k]` for `1 <= k <= k2`.  `TracesPower[P,{k1,k2}]` returns a list of traces `Tr[M^k]` for `k1 <= k <= k2`.

- `UnTally[L]` where `L` is a tallied list (see `Tally`) undoes `Tally`, or leaves `L` alone if already untallied.

## Internal helper functions

### ``Private`` context

- `iPolynomialTracesBounds`
- `iPolynomialVariable`
- `iReciprocalPolynomialTracesBounds`

### ``Lefschetz`Tests`` context

- `iAlmostPureStratumAQ`
- `iAlmostPureStratumBQ`
- `iLCMToString`
- `iLefschetzToString`
- `iMinimumSingularitiesQ`
- `iPureStratumAQ`
- `iPureStratumBQ`
- `iSingularityPermutationsAQ`
- `iSingularityPermutationsBQ1`
- `iSingularityPermutationsBQ`
- `iSingularityToString`
- `iStratumOrbitsTestQ`
- `iTestListQ`
- `iTestNegativeQ`
- `iTestPositiveQ`

### ``Lefschetz`Orbits`` context

- `iAllPossibilities`
- `iGroupByPartition`
- `iSingularCyclic`
- `iSumOrbits`

### ``Homology`` context

- `iDehnA`
- `iDehnB`
- `iDehnC`
