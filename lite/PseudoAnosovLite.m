(* PseudoAnosovLite Mathematica Package *)

(*
    Copyright 2010 Jean-Luc Thiffeault (jeanluc@mailaps.org)
                   Erwan Lanneau       (lanneau@cpt.univ-mrs.fr)

    This file is part of PseudoAnosov.

    PseudoAnosov is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    PseudoAnosov is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PseudoAnosov.  If not, see <http://www.gnu.org/licenses/>.
*)

BeginPackage["PseudoAnosovLite`"]


Quiet[Needs["Combinatorica`"],General::compat] (* For "ToCycles" *)


(*
   Command usage messages
*)

PseudoAnosov::usage = "Functions for manipulating characteristic polynomials of pseudo-Anosov maps."

PolynomialRoots::usage = "PolynomialRoots[P] returns the roots of the polynomial P, sorted in decreasing order of magnitude.";

PerronRoot::usage = "PerronRoot[P] returns the largest root (in magnitude) of the polynomial P.";

pseudoAnosovPerronRootQ::usage = "pseudoAnosovPerronRootQ[P] returns True if the largest root of the polynomial P is nondegenerate (in magnitude) and real.  pseudoAnosovPerronRootQ[P,xmax] also returns False unless the the Perron root is less than xmax (in magnitude).";

TracesPower::usage = "TracesPower[P,k], where P is the characteristic polynomial of some matrix M, returns the trace Tr[M^k].  TracesPower[P,{k2}] returns a list of traces Tr[M^k] for 1 <= k <= k2.  TracesPower[P,{k1,k2}] returns a list of traces Tr[M^k] for k1 <= k <= k2."

ReciprocalPolynomial::usage = "ReciprocalPolynomial[x,n] returns a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1.  ReciprocalPolynomial[x,n,c] uses c as the base name for coefficients.  ReciprocalPolynomial[x,n,{a_1,...,a_(n/2)}] uses a list for the coefficient, where (n/2) denotes Floor[n/2].";

ReciprocalPolynomialQ::usage = "ReciprocalPolynomialQ[P] returns true if P is a reciprocal polynomial, i.e. of the form a[0] x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + a[0].";

PolynomialFromTraces::usage = "PolynomialFromTraces[x,T,det] creates a polynomial of degree Length[T]+1 from the determinant det (defaults to 1) and a list T of traces of powers of its associated matrix."

ReciprocalPolynomialFromTraces::usage = "ReciprocalPolynomialFromTraces[x,n,T] creates a reciprocal polynomial of degree n from a list of traces of powers of its associated matrix.  ReciprocalPolynomialFromTraces[x,T] creates a polynomial of degree 2 Length[T]."

ReciprocalPolynomialBoundedList::usage = "ReciprocalPolynomialBoundedList[x,n,r] returns a list of reciprocal polynomials x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1 with Perron root less than r.  For n even, only one of each polynomial pair P(-x)=P(x) is listed.";

LefschetzNumbers::usage = "LefschetzNumbers[P,k], where P is the characteristic polynomial of some matrix M, returns the Lefschetz number 2-Tr[M^k].  LefschetzNumbers[P,{k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for 1 <= k <= k2.  LefschetzNumbers[P,{k1,k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for k1 <= k <= k2."

LefschetzCombine::usage = "LefschetzCombine[L1,L2,...] adds lists of Lefschetz number.  If they are not the same length, then the blocks are repeated to the length of the longest list.\nLefschetzCombine[L1,L2,...,Lm,n] caps the total length at an integer n."

StratumOrbits::usage = "StratumOrbits[S,P] returns a list of possible orbit structure (singular and regular periodic orbits) for the polynomial P on stratum S.  Returns an empty list if this proves impossible.\nStratumOrbits[S,L] does the same for a list of Lefschetz numbers L."

		    
(*
   Options
*)

MaxLefschetz::usage = "MaxLefschetz is an option to StratumOrbits, specifying how many Lefschetz numbers to compute."

PerronRootSign::usage = "PerronRootSign is an option to StratumOrbits (Lefschetz numbers form) to specify whether the Perron root is positive or negative.  Set to Automatic to try and guess by looking at the last two Lefschetz numbers (default Automatic)."

EqualityTolerance::usage = "EqualityTolerance is an option to pseudoAnosovPerronRootQ to decide whether two numbers are \"equal enough\"."

(* Rule labels for orbit structure of a stratum *)
Polynomial::usage = "Polynomial is a rule label to specify the polynomial in the output of StratumOrbits."
Stratum::usage = "Stratum is a rule label to specify the stratum (in tallied form) in the output of StratumOrbits."
SingularitiesPermutation::usage = "SingularitiesPermutation is a rule label in the output of StratumOrbits to specify the list of permutations on each singularity type of the stratum."
SeparatricesPermutation::usage = "SeparatricesPermutation is a rule label in the output of StratumOrbits to specify the list of permutations on the separatrices of a set of cyclically-permuted singularities.  This permutation applies only when all the singularities are fixed, and the iterate number is even."
SingularitiesLefschetzBlock::usage = "SingularitiesLefschetzBlock is a rule label in the output of StratumOrbits to specify the periodic block of Lefschetz numbers induced by the permutations of the singularities and their separatrices."
RegularOrbits::usage = "RegularOrbits is a rule label in the output of StratumOrbits to specify the list of regular periodic orbits.  Thus, a 3 in the fifth entry means 3 period-five orbits."

Protect[Polynomial,Stratum,SingularitiesPermutation,SeparatricesPermutation,SingularitiesLefschetzBlock,RegularOrbits]


Begin["`Private`"]


(* Helper function to get a polynomial's independent variable, with
   some checks *)
polynomialvariable[p_] := Module[{x = Variables[p]},
    If[x == {}, (* || Length[x] > 1, *)
        Message[polynomialvariable::notapolynomial,p,x]; Abort[]];
    Return[First[x]]
]
polynomialvariable::notapolynomial = "Input argument `1` is not a polynomial in `2`."


PolynomialRoots[p_,opts:OptionsPattern[]] := Module[
    {x = polynomialvariable[p]},
    Sort[x/.NSolve[p == 0, x, Sequence @@ FilterRules[{opts},Options[NSolve]]],
        Abs[#2] < Abs[#1] &]
]
(* PolynomialRoots inherits the options for NSolve *)
Options[PolynomialRoots] = Options[NSolve]


PerronRoot[p_,opts:OptionsPattern[]] := First[PolynomialRoots[p,opts]]
(* PerronRoot inherits the options for NSolve *)
Options[PerronRoot] = Options[NSolve]


(* Test for the Perron root *)
pseudoAnosovPerronRootQ[p_,lmax___:0,opts:OptionsPattern[]] := Module[
    {prl, pr, degen, testdegen},
    degen = OptionValue[EqualityTolerance];
    testdegen[diff_] := Module[{},
        If[diff < degen, Return[False]];
        If[degen/2 < diff < degen, Message[PseudoAnosov::tooclose]];
        Return[True]
    ];
    (* Take the first two roots (presorted by magnitude) *)
    (* The way options work in Mathematica is tedious: just because I
       set the default value to MachinePrecision below, it doesn't
       get passed on to PolynomialRoots automatically. *)
    prl = Take[PolynomialRoots[p,
        WorkingPrecision -> OptionValue[WorkingPrecision]]];
    pr = Abs[prl[[1]]];
    (* First criterion: largest eigenvalue is greater than 1 *)
    If[!testdegen[Abs[pr-1]], Return[False]];
    (* Second criterion: largest eigenvalue is separated from the second *)
    If[!testdegen[pr - Abs[prl[[2]]]], Return[False]];
    If[lmax != 0,
        (* Third criterion: largest eigenvalue is less than lmax *)
        If[!testdegen[Abs[lmax] - pr], Return[False]];
    ];
    Return[True];
]
Options[pseudoAnosovPerronRootQ] =
    {WorkingPrecision -> MachinePrecision, EqualityTolerance -> 10^-8}


TracesPower[p_,mm_List:{10}] := Module[
    {x = polynomialvariable[p], n, T, ml},
    (* Polynomial is x^n + c[1]x^(n-1) + ... + c[n-1]x + c[n] *)
    c = Reverse[CoefficientList[Collect[p,x],x]];
    (* Make sure leading coefficient is 1, then drop it. *)
    c = Drop[c/c[[1]],1];
    n = Length[c]; (* Degree of polynomial *)
    If[Length[mm] == 1, ml = {1,mm[[1]]}, ml = mm];
    (* Recursively compute the traces. *)
    (*   Note that I used to have Table instead of the Do, but
         Mathematica must sometimes construct tables out of order,
         because this sometimes caused a crash. *)
    (* Ideally, the function would not store the T's that are not requested. *)
    Do[
        T[k] = -Sum[c[[m]] T[k-m], {m,Min[k-1,n]}]
               - If[k > n, 0, k c[[k]]]
    ,{k,ml[[2]]}];
    Table[T[k],{k,ml[[1]],ml[[2]]}]
]


TracesPower[p_,m_Integer:1] := First[TracesPower[p,{m,m}]]


ReciprocalPolynomial[x_,n_,a_List] := Module[{aal},
    aal = Table[a[[m]], {m,n/2}];
    aal = If[OddQ[n], Join[aal,Reverse[aal]], Join[aal,Drop[Reverse[aal],1]]];
    aal = Append[Prepend[aal,1],1];
    x^n + Sum[aal[[k]] x^(k-1), {k,n}]
]


ReciprocalPolynomial[x_,n_,c_:Global`a] :=
    ReciprocalPolynomial[x,n,Table[c[k],{k,n/2}]]


ReciprocalPolynomialQ[p_] := Module[
    {x = polynomialvariable[p], c, n},
    c = CoefficientList[Collect[p,x],x];
    n = Length[c]-1;
    Return[Simplify[p - x^n (p/.x->1/x)] === 0]
]


PolynomialFromTraces[x_,T_List,det_:1] := Module[{c, n = Length[T]+1},
    Do[c[k] = (-T[[k]] - Sum[c[m] T[[k-m]],{m,k-1}])/k,{k,n-1}];
    x^n + Sum[c[k] x^(n-k), {k,n-1}] + det
]


ReciprocalPolynomialFromTraces[x_,T_List] :=
    ReciprocalPolynomialFromTraces[x,2 Length[T],T]


ReciprocalPolynomialFromTraces[x_,n_Integer,T_List] := Module[{a},
    If[n > 2 Length[T] + 1,
        Message[ReciprocalPolynomialFromTraces::toofewtraces]; Return[]];
    Do[a[k] = (-T[[k]] - Sum[a[m] T[[k-m]],{m,k-1}])/k,{k,Floor[n/2]}];
    ReciprocalPolynomial[x,n,a]
]
ReciprocalPolynomialFromTraces::toofewtraces = "Error: list of traces should be at least n/2, where n is the degree of the polynomial."


iReciprocalPolynomialTracesBounds[n_Integer,r_] :=
    Floor[n/2 (r^# + r^-#)] & /@ Range[n/2] /; EvenQ[n]


iReciprocalPolynomialTracesBounds[n_Integer,r_] :=
    Floor[(n-1)/2 (r^# + r^-#) + 1] & /@ Range[n/2] /; OddQ[n]


ReciprocalPolynomialBoundedList[x_,n_Integer,r_,opts:OptionsPattern[]] :=
Module[
    {p,pl,sl,T,Tm = iReciprocalPolynomialTracesBounds[n,r]},
    p = ReciprocalPolynomialFromTraces[x,n,Table[T[k],{k,Floor[n/2]}]];
    pl = Flatten[Fold[
        Table[#1,{T[#2],-Tm[[#2]],Tm[[#2]]}]&, p, Range[Floor[n/2]]
    ]];
    (* Discard the polynomials that don't have integer coefficients *)
    sl = And @@ # & /@  (IntegerQ /@ # & /@ (CoefficientList[#,x] & /@ pl));
    pl = Pick[pl,sl];
    (* Discard the ones without the proper Perron root (not real or
       too large) *)
    pl = Pick[pl,pseudoAnosovPerronRootQ[#,r,opts] & /@ pl];
    (* Make all roots positive, discard duplicates *)
    (* Only do this for even n: for odd n, the leading term changes sign. *)
    If[EvenQ[n],
      pl = Union[If[PerronRoot[#] > 0, #, #/.x->-x] & /@ pl]
    ];
    (* Sort by Perron root *)
    (* Note that we computed the Perron root many times: a waste *)
    Sort[pl, PerronRoot[#1] < PerronRoot[#2] &]
]
Options[ReciprocalPolynomialBoundedList] = Options[pseudoAnosovPerronRootQ]


End[(* "`Private`" *)]


Begin["`Lefschetz`"]

(*
   General helper functions and definitions
*)

LefschetzNumbers[p_,mm_List] := (2 - #) & /@ TracesPower[p,mm]


LefschetzNumbers[p_,m_Integer:1] := 2 - TracesPower[p,m]


LefschetzCombine[Ll__List, n_Integer:0] := Module[{L = List[Ll], len},
    If[n == 0, len = LCM @@ (Length /@ L), len = n];
    Plus @@ (PadRight[#,len,#] & /@ L)
]


Begin["`Orbits`"]


allpossibilities::usage = "allpossibilities[L1,L2,...] generates all possible ordered combinations of elements of the given lists with repetition."

allpossibilities[L__List, OptionsPattern[]] := Module[{f, pos},
    (* The function f helps generate all possible combinations of
       elements of given lists *)
    f[l_, a_] := Flatten[Table[
         Append[l[[i]], a[[j]]], {i, Length[l]}, {j, Length[a]}], 1];
    (* Use Fold to apply f to each list, giving all possibilities *)
    pos = Fold[f,{{}},List[L]];
    (* If Ordered is False, order doesn't matter *)
    If[OptionValue[Ordered], pos, Union[Sort/@pos]]
]
Options[allpossibilities] = {Ordered -> True}


sumorbits::usage = "sumorbits[p,l] with p an integer and l a list of nonnegative integers, returns the sum of (k l[[k]]), where k is a divisor of p."

sumorbits[p_Integer, rpo_List, OptionsPattern[]] := Module[
    {dl = Divisors[p]},
    If[!OptionValue[IncludeLast], dl = Most[dl]];
    Plus @@ (# rpo[[#]] & /@ dl)
]
Options[sumorbits] = {IncludeLast -> True}

IncludeLast::usage = "Option to sumorbits: set to True to include the final iterate's contribution to the sum from (default True)."


(* Groups the list l according to an integer partition of Length[l].
   Example: groupbypartition[{1,2,3,4},{1,3}] = {{1},{2,3,4}} *)
groupbypartition[l_, part_] := Module[{g = {}},
    Fold[(AppendTo[g, Take[#1, #2]]; Drop[#1, #2]) &, l, part];
    Return[g]
]


(*
   Construct regular orbits from Lefschetz numbers
*)

(* Must list this function before the next one to ensure that the more
   "specific" version is matched first *)
StratumOrbits[s_List,L_List, opts:OptionsPattern[]] := Module[
    {Ls, ro, prs = OptionValue[PerronRootSign], opts2, len},
    (* If the sign of the Perron root is unspecified, try and guess *)
    If[prs == Automatic, prs = Sign[Times @@ Take[L,-2]]];
    Ls = SingularStratum[s,prs];
    Off[Regular::badLefschetz];
    opts2 = Sequence @@ FilterRules[{opts},Options[Regular]];
    len = Min[Length[L],OptionValue[MaxLefschetz]];
    ro = Regular[
        LefschetzCombine[Take[L,len],
            -SingularitiesLefschetzBlock/.#,len],opts2] & /@ Ls;
    On[Regular::badLefschetz];
    (* Eliminate bad orbits: look for a negative or nonintegral last element *)
    ro = Transpose[{Ls,ro}];
    ro = Pick[ro, Last[#[[2]]] > 0 && IntegerQ[Last[#[[2]]]] & /@ ro];
    Join[#[[1]],{RegularOrbits -> #[[2]]}] & /@ ro
]
Options[StratumOrbits] = {PerronRootSign -> Automatic, MaxLefschetz -> 50}


StratumOrbits[s_List,p_, opts:OptionsPattern[]] := Module[
    {L, x = polynomialvariable[p], prs = opts2},
    (* Find the sign of the Perron root, add it to options *)
    prs = Sign[PerronRoot[p]];
    If[OptionValue[PerronRootSign] != Automatic,
        If[OptionValue[PerronRootSign] != prs,
            Message[StratumOrbits::noprs]; Abort[]
        ]
    ];
    opts2 = Sequence @@ Join[{opts},{PerronRootSign -> prs}];
    L = LefschetzNumbers[p,{OptionValue[MaxLefschetz]}];
    Prepend[#,Polynomial -> p]& /@ StratumOrbits[s,L,opts2]
] /; ReciprocalPolynomialQ[p]
StratumOrbits::noprs = "PerronRootSign option contradicts actual Perron root of polynomial.  Do not use PerronRootSign with a polynomial."


SingularStratum::usage = "PseudoAnosov`Lefschetz`Orbits`SingularStratum[S,prs], where S is a list of the degrees of singularities in a stratum, returns a list of all possible Lefschetz number sequences corresponding to the singularities, without taking into account regular orbits.  S can be specified as an explicit list (i.e., {4,2,2,2}) or in tallied form ({{4,1},{2,3}}).  The sign of the Perron root is given by prs (default -1)."

SingularStratum[s_List, prs_Integer:-1] := Module[
    {t, L, f},
    (* Accept s in either explicit or tallied form *)
    If[ListQ[s[[1]]], t = s, t = Tally[s]];
    (* Lists of Lefschetz numbers for each singularity type *)
    L = ((Singular[##,prs]&) @@ # &) /@ t;
    (* Take all possibilities *)
    L = allpossibilities @@ L;
    L = Transpose /@ L;
    (* Combine the lists of Lefschetz, repeating blocks as needed *)
    L = {Reverse[Transpose[#[[1]]]],LefschetzCombine @@ #[[2]]} & /@ L;
    (* Add some keywords *)
    {Stratum -> t,
     SingularitiesPermutation -> #[[1,1]],
     SeparatricesPermutation -> #[[1,2]],
     SingularitiesLefschetzBlock -> #[[2]]} & /@ L
]


Regular::usage = "PseudoAnosov`Lefschetz`Orbits`Regular[L], where L is a list of Lefschetz numbers, returns the list of regular periodic orbits compatible with L, unless an incompatible orbit is detected, in which cases the function stops and returns what it's found.  The sign of the Perron root can be specified by the option PerronRootSign (default Automatic)."

Regular[L_List, OptionsPattern[]] := Module[
    {rpo, def, prs = OptionValue[PerronRootSign]},
    (* If the sign of the Perron root is unspecified, try and guess *)
    If[prs == Automatic, prs = Sign[Times @@ Take[L,-2]]];
    rpo = {-prs L[[1]]};
    If[rpo[[1]] < 0, Message[Regular::badLefschetz]; Return[rpo]];
    Do[
        If[prs < 0,
            def = (-1)^(p+1) L[[p]] - sumorbits[p,rpo,IncludeLast->False]
        ,
            def = -L[[p]] - sumorbits[p,rpo,IncludeLast->False]
        ];
        AppendTo[rpo, def/p];
        If[!IntegerQ[def/p] || def < 0,
            Message[Regular::badLefschetz]; Break[]
        ]
    ,{p, 2, Length[L]}];
    rpo
]
Options[Regular] = {PerronRootSign -> Automatic}
Regular::badLefschetz = "Bad sequence of Lefschetz numbers."


Singular::usage = "PseudoAnosov`Lefschetz`Orbits`Singular[k,m,prs], where k and m are integers, lists all possible Lefschetz number sequences for m singularities of degree k (m defaults to 1).  The sign of the Perron root is given by prs (default -1)."

Singular[k_Integer, m_Integer:1, prs_Integer:-1] :=
Module[
    {pr = k/2+1, Pk, Pm, clen, blk, alpo, all},
    (* All possible cyclic permutations of separatrices *)
    Pk = RotateLeft[Range[pr], #] & /@ (Range[pr] - 1);
    (* Keep only the ones that have different cycle lengths *)
    clen[perm_] := Sort[Length/@ToCycles[perm]];
    Pk = Union[Pk, SameTest -> (clen[#1] == clen[#2] &)];
    (* Possible permutations of singularities *)
    Pm = groupbypartition[Range[m],#] & /@ IntegerPartitions[m];
    blk = Tally /@ ((Length/@#)&/@Pm);
    all = {};
    (* Loop over each possible Pm *)
    Do[
        (* Loop over each cyclic block length of Pm
             blk[[i,j,1]] gives the block length
             blk[[i,j,2]] gives the multiplicity (# of blocks of that length)
           Generate a list of combinations of permutations,
             allowing for symmetry do to multiplicity *)
        alpo = ((allpossibilities[##,Ordered->False]&) @@
                    Table[Pk,{#[[2]]}]& ) /@ blk[[i]];
        (* Now combine the multiplicities together *)
        alpo = allpossibilities[##]& @@ alpo;
        alpo = {#,Pm[[i]]}& /@ (Flatten[#,1]& /@ alpo);
        all = Join[all,alpo];
    ,{i,Length[Pm]}];
    (* Return a list of entries of the form
       {permutations,Lefschetznumbers} *)
    ({{##},SingularPermutations[##,prs]}&) @@ # & /@ all
]


SingularPermutations::usage = "PseudoAnosov`Lefschetz`Orbits`SingularPermutations[Pk,Pm,prs] returns a list of Lefschetz numbers corresponding to singularity of degree k with m-fold degeneracy.  Pk is a list of permutations on the (k+2)/2 in or outgoing separatrices of the singularities, and Pm a permutation on the m singularities.  Pm must be in cycles form.  Note that Pk must be a power of a cyclic permutation.  The sign of the Perron root is given by prs (default -1)."

SingularPermutations[Pk_List, Pm_List:{{1}}, prs_Integer:-1] := Module[
    {k = 2Length[Pk]-2, m = Length[Flatten[Pm]], Pmc = Length /@ Pm},
    If[Length[Pk] != Length[Pmc], Message[PseudoAnosov::bad]];
    (* Apply SingularCyclic to each subcycle, with Pk specifying a
       permutation on the separatrices for each subcycle *)
    LefschetzCombine @@
        ((SingularCyclic[##,prs]&) @@ # & /@ Transpose[{Pk,Pmc}])
]


(* Private helper function for Singular *)
(* Special case of Singular for the m singularities
   being permuted cyclically *)
SingularCyclic[Pk_List, m_Integer:1, prs_Integer:-1] := Module[
    {k = 2Length[Pk]-2, period, L, Pm, Pm1, Pk1 = Pk},
    (* Generate a cyclic permutation for the singularities *)
    Pm = RotateLeft[Range[m]]; Pm1 = Pm;
    (* The total period for the separatrices *)
    period = m Times @@ Union[Length /@ ToCycles[Pk]];
    If[prs < 0,
        (* If the period is odd, then double it to get back to an even
           iterate *)
        period = If[OddQ[m], 2period, period];
    ];
    L = {};
    Do[
        If[Pm1 == Range[m],
            If[Pk1 == Range[k/2+1],
	        If[prs < 0,
                    AppendTo[L,If[EvenQ[i], -m(k+1), m]]
		,
                    AppendTo[L,-m(k+1)]
		]
            ,
                AppendTo[L,m]
            ];
            (* Permute the separatrices (only on even iterates if prs < 0) *)
            If[(prs < 0 && EvenQ[i]) || prs > 0, Pk1 = Permute[Pk1,Pk]];
        ,
            AppendTo[L,0]
        ];
        Pm1 = Permute[Pm1,Pm];
    , {i, period}];
    Return[L]
]


End[(* "`Orbits`" *)]


End[(* "`Lefschetz`" *)]


EndPackage[(* "PseudoAnosovLite`" *)]
