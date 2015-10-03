# PseudoAnosov Package

The PseudoAnosov package contains **Mathematica** and **C++** functions for extracting properties of characteristic polynomials of pseudo-Anosov maps of surfaces.  The goal is to determine whether a given polynomial is *allowable* on a given stratum.  A stratum is a collection of singularities for a closed surface of genus *g*.  A polynomial is allowable if it can (at least in principle) arise as the characteristic polynomial associated with a pseudo-Anosov map on that surface, where the map preserves the singularities up to permutations.

Most of the functions in the package assume that the pseudo-Anosov map stabilizes an *orientable* foliation, so the singularities must all have even degree.  If the pseudo-Anosov map exists, then its *dilatation* is given by the largest root of the polynomial, also called its *Perron root*.

This package is an extension of the code used in papers by [Erwan Lanneau][1] and [Jean-Luc Thiffeault][2] to prove lower bounds on the smallest dilatation on [closed surfaces][3] and [braids][4].  A minimal version was included as part of the [braids paper][4] as [PseudoAnosovLite](/lite/).

The PseudoAnosov package is maintained by [Jean-Luc Thiffeault][2].

See the [complete list of functions](/functions.md/) in the PseudoAnosov package.  There are also several sample notebooks in the repository.

### A sample Mathematica session

First we load the package:
```mathematica
In[1]:= <<PseudoAnosov.m
```
Let's start with the torus (genus `g=1`), where we are dealing with Anosov maps (no singularities).  We can list all the reciprocal monic polynomials of degree `2g=2` (with integer coefficients) with dilatation less than, say, 3:
```mathematica
In[2]:= P = ReciprocalPolynomialBoundedList[x,2,3]

                    2
Out[2]= {1 - 3 x + x }
```
There is only one such polynomial, and it gives the lowest dilatation on that surface.

Genus 2 surfaces are more interesting.  First we can list all the possible strata (singularity data) for the case were the foliation is orientable:
```mathematica
In[3]:= s = OrientableStrataList[2]

Out[3]= {{4}, {2, 2}}
```
There are two strata, one with a degree-4 singularity (6 prongs), the other with two degree-2 singularities (4 prongs).  The sum of the degrees in each case is `4=4g-4`=-2(Euler characteristic), where `g=2` is the genus.  Next we find all monic reciprocal polynomials of degree `2g=4` with dilatations less than 1.75:
```mathematica
In[4]:= P = ReciprocalPolynomialBoundedList[x,4,1.75]

                  2    3    4
Out[4]= {1 - x - x  - x  + x }
```
There is only one such polynomial (which is why we chose the value 1.75).  Its root with the largest magnitude (Perron root) is
```mathematica
In[5]:= PerronRoot /@ P

Out[5]= {1.72208}
```
Let's see if this polynomial is a good candidate for a pseudo-Anosov map.  The command `StratumOrbits` tries to find a set of regular periodic orbits that is compatible with a polynomial on a given stratum.  If we try the polynomial `P[[1]]` on the second stratum `s[[2]]`,
```mathematica
In[6]:= StratumOrbits[s[[2]],P[[1]]]

Out[6]= {}
```
we get nothing, suggesting that this polynomial is not associated with any pseudo-Anosov maps on the stratum `s[[2]]={2,2}`.  However, on the stratum `s[[1]]={4}` we get
```mathematica
In[7]:= so = StratumOrbits[s[[1]],P[[1]]]

                                 2    3    4
Out[7]= {{Polynomial -> 1 - x - x  - x  + x , Stratum -> {{4, 1}},
     SingularitiesPermutation -> {{{1}}},
     SeparatricesPermutation -> {{{2, 3, 1}}},
     SingularitiesLefschetzBlock -> {1, 1, -5},
     RegularOrbits ->
      {0, 1, 0, 1, 3, 3, 6, 9, 14, 21, 36, 54, 90, 141, 230, 369, 606, 977, 1608, 2619, 4312, 7074, 11682, 19248, 31872, 52731, 87514, 145260, 241644, 402137, 670380, 1118187, 1867560, 3121221, 5221938, 8742312, 14648958, 24562068, 41214696, 69199515, 116263056, 195445504, 328749954, 553264722, 931601482, 1569414123, 2645169030, 4460292930, 7524259626, 12698241600}}}
```
This is a candidate pseudo-Anosov, since there exists a collection of regular periodic orbits that together with the singularity respect the [Lefschetz fixed-point theorem][5].  We can present the information better with
```mathematica
In[8]:= StratumOrbitsTable[so[[1]]]

Out[8]=
```
![output of StratumOrbitsTable[so[[1]]]](/images/sotable.png/)

The table shows the polynomial and its Perron root.  The stratum is denoted `4^1` or `{{4,1}}`, which means degree `4` with multiplicity `1`, i.e., the same as `{4}`.  There is only one singularity, so the permutation of singularities by the map can only be `{1}`.  The singularity has six prongs or separatrices, 3 of which are labeled 'ingoing' and the other 3 'outgoing'.  These ingoing/outgoing separatrices are cyclically permuted as `{2,3,1}` by this (hypothetical) pseudo-Anosov map.  The table then gives Lefschetz number sequences for iterates `n` of the map.  The last row gives the number of regular orbits (`#ro`).  We see that there are no regular fixed points (`n=1`), one period-2 orbit (`n=2`), and so on.

It turns out that a pseudo-Anosov map having this polynomial does indeed exist, and we just deduced that it must have the minimum dilatation for a genus 2 surface, as was first shown by Zhirov (1995), since there are no candidate polynomials with a lower dilatation.

### License

PseudoAnosov is released under the [GNU General Public License v3][6].  See [COPYING](/COPYING/) and [LICENSE](/LICENSE/).

### Support

The development of PseudoAnosov was supported by the [US National Science Foundation][7], under grants [DMS-0806821][8] and [CMMI-1233935][9].

[1]: https://www-fourier.ujf-grenoble.fr/~lanneau/
[2]: http://www.math.wisc.edu/~jeanluc/
[3]: http://arxiv.org/abs/0905.1302 "On the minimum dilatation of pseudo-Anosov homeomorphisms on surfaces of small genus, Annales de l'Institut Fourier 61, 105–144, 2011"
[4]: http://arxiv.org/abs/1004.5344 "On the minimum dilatation of braids on the punctured disc, Geometriae Dedicata 152, 165–182, 2011."
[5]: https://en.wikipedia.org/wiki/Lefschetz_fixed-point_theorem
[6]: http://www.gnu.org/licenses/gpl-3.0.html
[7]: http://www.nsf.gov
[8]: http://www.nsf.gov/awardsearch/showAward?AWD_ID=0806821
[9]: http://www.nsf.gov/awardsearch/showAward?AWD_ID=1233935
