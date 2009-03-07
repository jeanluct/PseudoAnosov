(* PseudoAnosov Mathematica Package *)

BeginPackage["PseudoAnosov`", {"Global`"}]

Needs["Combinatorica`"]


(*
   Command usage messages
*)

PseudoAnosov::usage = "Functions for manipulating characteristic polynomials of pseudo-Anosov maps."

PolynomialDegree::usage = "PolynomialDegree[P] returns the degree of the polynomial P(x)."

TracesPower::usage = "TracesPower[P,m], where P is the characteristic polynomial of a matrix M, lists the traces Tr[M^k] for 1 <= k <= m."

LefschetzNumbers::usage = "LefschetzNumbers[P,k], where P is the characteristic polynomial of some matrix M, returns the Lefschetz number 2-Tr[M^k].  LefschetzNumbers[P,{k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for 1 <= k <= k2.  LefschetzNumbers[P,{k1,k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for k1 <= k <= k2."

ReciprocalPolynomial::usage = "ReciprocalPolynomial[x,n] returns a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1.  ReciprocalPolynomial[x,n,c] uses c as the base name for coefficients.  ReciprocalPolynomial[x,n,{a_1,...,a_(n/2)}] uses a list for the coefficient, where (n/2) denotes Floor[n/2].";

ReciprocalPolynomialFromTraces::usage = "ReciprocalPolynomialFromTraces[x,n,T] creates a reciprocal polynomial of degree n from a list of traces of powers of its associated matrix.  ReciprocalPolynomialFromTraces[x,T] creates a polynomial of degree 2 Length[T]."

ReciprocalPolynomialQ::usage = "ReciprocalPolynomialQ[P] returns true if P is a reciprocal polynomial, i.e. of the form a[0] x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + a[0].";

ReciprocalPolynomialCoefficientBounds::usage = "ReciprocalPolynomialCoefficientBounds[n,t] lists the bound on the magnitude of the coefficients {a[1],...,a[Floor[n/2]]} of a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1, given that its largest eigenvalue is h with t = h + 1/h.";

ReciprocalPolynomialBoundedList::usage = "ReciprocalPolynomialBoundedList[x,n,amax] returns a list of reciprocal polynomials x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1 with coefficients bounded by |a[k]| <= amax[k].  For n even, only one of each polynomials pair P(-x)=P(x) is listed.";

PolynomialRoots::usage = "PolynomialRoots[P] returns the roots of the polynomial P, sorted in decreasing order of magnitude.";

PerronRoot::usage = "PerronRoot[P] returns the largest root (in magnitute) of the polynomial P.";

pseudoAnosovPerronRootQ::usage = "pseudoAnosovPerronRootQ[P] returns True if the largest root of the polynomial P is nondegenerate (in magnitute) and real.  pseudoAnosovPerronRootQ[P,xmax] also returns False unless the the Perron root is less than xmax (in magnitude).";

MahlerMeasure::usage = "MahlerMeasure[P] returns the Mahler measure of the polynomial P, which is the absolute value of the product of roots outside the unit circle.";

IrreducibleMatrixQ::usage = "IrreducibleMatrixQ[M] returns true if the matrix M is irreducible.";

(* StrataList::usage = "StrataList[g,n] gives the list of strata for a ...
 hyperbolic surface of genus g with n boundary components (default n=0).  If n>0, we need at least one singularity per boundary component.  Note that a p-pronged singularity on the boundary is really a (p-2)-prong when the boundary components are shrunk to punctures.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity."; *)

OrientableStrataList::usage = "OrientableStrataList[g] gives the list of orientable strata for a hyperbolic surface of genus g.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the (even) degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity.  Use Tally/@OrientableStrataList[g] to group by multiplicities.";

StratumToGenus::usage = "StratumToGenus[S] gives the genus of the surface containing a stratum S={k_1,...,k_m}."

SumOrbits::usage = "SumOrbits[p,l] with p an integer and l a list of nonnegative integers, returns the sum of (k l[[k]]), where k is a divisor of p."

AllPossibilities::usage = "AllPossibilities[L1,L2,...] generates all possible ordered combinations of elements of the given lists with repetition."

LefschetzRegularOrbits::usage = ""

LefschetzNumbersSingularity::usage = "LefschetzNumbersSingularity[k,m], where k and m are integers, lists all possible Lefschetz number sequences for m singularities of degree k (m defaults to 1)."

LefschetzNumbersSingularityPermutations::usage = ""

LefschetzNumbersStratum::usage = "LefschetzNumbersStratum[S], where S is a list of the degrees of singularities in a stratum, returns a list of all possible Lefschetz number sequences corresponding to the singularities, without taking into account regular orbits.  S can be specified as an explicit list (i.e., {4,2,2,2}) or in tallied form ({{4,1},{2,3}})."

LefschetzCombine::usage = "LefschetzCombine[L1,L2,...] adds lists of Lefschetz number.  If they are not the same length, then the blocks are repeated to the length of the longest list.\nLefschetzCombine[L1,L2,...,Lm,n] caps the total length at an integer n."

LefschetzNumbersTestQ::usage = ""

StratumRegularOrbits::usage = "StratumRegularOrbits[S,P] returns a list of regular periodic orbits for the polynomial P on stratum S; the list is empty if this proves impossible."

(* Test for everything: first call pseudoAnosovPerronRootQ, then
LefschetzNumbersTestQ by strata (for even power). *)
(* pseudoAnosovPolynomialQ::usage = "" *)


(*
   Options
*)

GiveReasonForRejection::usage = "Option to LefschetzNumbersTestQ: Set to True to return the reason for rejecting a stratum (default False)."

MaxIterate::usage = "Option to LefschetzNumbersTestQ: Set to an integer giving the largest power of the map to test (default 10)."

MaxLefschetz::usage = "Option to LefschetzNumbersTestQ: Set to an integer giving how many of Lefschetz numbers to compute for the test."

IncludeLast::usage = "Option to SumOrbits: set to True to include the final iterate's contribution to the sum from (default True)."

PerronRootSign::usage = "Option to LefschetzStratum, LefschetzSingularity and LefschetzRegularOrbits to specify whether the Perron root is positive or negative.  For LefschetzRegularOrbits, set to Automatic to try and guess by looking at the last two Lefschetz numbers (default Automatic for LefschetzRegularOrbits, -1 otherwise)."


(*
   Error messages and warnings
*)

PseudoAnosov::toofewtraces = "Error: list of traces should ne at least n/2, where n is the degree of the polynomial."

PseudoAnosov::notminimal = "Warning: Not a minimal polynomial."

PseudoAnosov::nottested = "Warning: This function is not well tested."

PseudoAnosov::needpositivePerron = "Error: This function only applies to positive Perron root."

PseudoAnosov::neednegativePerron = "Error: This function only applies to negative Perron root."

PseudoAnosov::moreLefschetz = "Need at least `1` Lefschetz numbers at `2`th power."

PseudoAnosov::notastring = "Function `1` did not return a proper string."

PseudoAnosov::badLefschetz = "Bad sequence of Lefschetz numbers."

PseudoAnosov::manyallowableperms = "More than one allowable permutation type on stratum `1`."

PseudoAnosov::notapolynomial = "Input argument `1` is not a polynomial in `2`."


Begin["`Private`"]


(* Helper function to get a polynomial's independent variable, with
   some checks *)
PolynomialVariable[p_] := Module[{x = Variables[p]},
    If[x == {} || Length[x] > 1,
        Message[PseudoAnosov::notapolynomial,p,x]; Abort[]];
    Return[First[x]]
]


PolynomialDegree[p_] := Module[{x = PolynomialVariable[p]},
    Length[CoefficientList[Collect[p,x],x]]-1
]


TracesPower[p_,mm_List:{10}] := Module[
    {x = PolynomialVariable[p], n, T, ml},
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


LefschetzNumbers[p_,mm_List] := (2 - #) & /@ TracesPower[p,mm]


LefschetzNumbers[p_,m_Integer:1] := 2 - TracesPower[p,m]


ReciprocalPolynomial[x_,n_,a_List] := Module[{aal},
    aal = Table[a[[m]], {m,n/2}];
    aal = If[OddQ[n], Join[aal,Reverse[aal]], Join[aal,Drop[Reverse[aal],1]]];
    aal = Append[Prepend[aal,1],1];
    x^n + Sum[aal[[k]] x^(k-1), {k,n}]
]


ReciprocalPolynomial[x_,n_,c_:Global`a] :=
    ReciprocalPolynomial[x,n,Table[c[k],{k,n/2}]]


ReciprocalPolynomialFromTraces[x_,T_List] :=
    ReciprocalPolynomialFromTraces[x,2 Length[T],T]


ReciprocalPolynomialFromTraces[x_,n_Integer,T_List] := Module[
    {tl, a},
    If[n > 2 Length[T],
        Message[PseudoAnosov::toofewtraces]; Return[]];
    tl = TracesPower[ReciprocalPolynomial[x,n,a],x,{n/2}];
    ReciprocalPolynomial[x,n,
        Table[a[k],{k,n/2}] /.
        Solve[Table[T[[k]] == tl[[k]],{k,n/2}], Table[a[k],{k,n/2}]][[1]]]
]


ReciprocalPolynomialQ[p_] := Module[
    {x = PolynomialVariable[p], c, n},
    c = CoefficientList[Collect[p,x],x];
    n = Length[c]-1;
    Return[Simplify[p - x^n (p/.x->1/x)] === 0]
]


ReciprocalPolynomialCoefficientBounds[n_,t_] := Module[
    {poly,h,hh,c,ssol},
    (* Make a reciprocal polynomial with eigenvalues h and 1/h,
       with Floor[n/2]-fold degeneracy each. *)
    poly = Collect[((x-h)(x-1/h))^Floor[n/2],x,Simplify];
    (* If odd, need an extra -1 eigenvalue to make it reciprocal. *)
    If[OddQ[n], poly = poly (x+1)];
    c = CoefficientList[poly,x];
    (* Solve for h, where h + 1/h = t.  Keep positive solution. *)
    hsub = Solve[t == h + 1/h,h][[2,1]];
    (* Finally, take coefficients and assume all eigenvalues are as
       large as possible. *)
    Simplify[Table[(-1)^k c[[k+1]]/.hsub,{k,n/2}]]
]


(* Very kludgy: eliminate 'by hand' the duplicate polynomials that are
   related by P(X)=P(-X), since these have the same dominant
   eigenvalue.  Only hand-coded for n=2,4,6,8 for now. *)
(* To do: implement the elimination of P(-x) for any even n. *)
ReciprocalPolynomialBoundedList[x_,n_,a_List] := Module[
    {p,c,l,l2},
    p = ReciprocalPolynomial[x,n,c];
    If[n == 2 || n == 4 || (n > 8 && EvenQ[n]),
        l = Flatten[Fold[
	    Table[#1,{c[#2],-a[[#2]],a[[#2]]}]&, p, Table[k,{k,2,n/2}]
        ]];
        Flatten[Table[l,{c[1],-a[[1]],0}]]
    ,
    If[n == 6 || n == 8,
        (* First list the cases with c[1]<0 *)
        l = Table[p,{c[1],-a[[1]],-1}];
        l = Flatten[Fold[
	    Table[#1,{c[#2],-a[[#2]],a[[#2]]}]&, l, Table[k,{k,2,n/2}]
        ]];
        (* Then list the cases with c[1]=0, c[3]<=0 *)
        l2 = Table[p/.c[1]->0,{c[3],-a[[3]],0},{c[2],-a[[2]],a[[2]]}];
        l2 = Flatten[Fold[
	    Table[#1,{c[#2],-a[[#2]],a[[#2]]}]&, l2, Table[k,{k,4,n/2}]]];
	Join[l,l2]
        ,
        (* If all else fails, just include everything. *)
	Flatten[Fold[
	    Table[#1,{c[#2],-a[[#2]],a[[#2]]}]&, p, Table[k,{k,n/2}]
        ]]
    ]
    ]
]

(* These versions of the function takes a dilatation as an argument
   rather than a coefficient list.  They are hand-coded for specific
   n, which is ugly but this function is only needed for specific,
   difficult cases. *)
ReciprocalPolynomialBoundedList[x_,6,h_Real] := Module[
    {p,n = 6,c,a,pl = {},hh,hmin = 1+10^-4, degen = 10^-5},
    a = Floor/@ReciprocalPolynomialCoefficientBounds[n,h+1/h];
    p = ReciprocalPolynomial[x,n,c];
    keep[hh_] := ((hmin < hh[[1]] < h) && Abs[hh[[1]]-hh[[2]]] > degen);
    Do[
        hh = Abs/@Take[PolynomialRoots[p],2];
        If[keep[hh], pl=Append[pl,p]; Print[p]]
    ,{c[1],-a[[1]],-1},{c[2],-a[[2]],a[[2]]},{c[3],-a[[3]],a[[3]]}];
    p = p/.{c[1]->0};
    Do[
        hh = Abs/@Take[PolynomialRoots[p],2];
        If[keep[hh], pl=Append[pl,p]; Print[p]]
    ,{c[2],-a[[2]],a[[2]]},{c[3],-a[[3]],0}];
    pl
]
ReciprocalPolynomialBoundedList[x_,8,h_Real] := Module[
    {p,n = 8,c,a,pl = {},hh,hmin = 1+10^-4, degen = 10^-5, cases = 0},
    a = Floor/@ReciprocalPolynomialCoefficientBounds[n,h+1/h];
    p = ReciprocalPolynomial[x,n,c];
    keep[hh_] := ((hmin < hh[[1]] < h) && Abs[hh[[1]]-hh[[2]]] > degen);
    Do[
        If[Mod[++cases,10000] == 0, Print[cases]];
        hh = Abs/@Take[PolynomialRoots[p],2];
        If[keep[hh], pl=Append[pl,p]; Print[p]]
    ,{c[1],-a[[1]],-1},{c[2],-a[[2]],a[[2]]}
    ,{c[3],-a[[3]],a[[3]]},{c[4],-a[[4]],a[[4]]}];
    p = p/.{c[1]->0};
    Do[
        If[Mod[++cases,10000] == 0, Print[cases]];
        hh = Abs/@Take[PolynomialRoots[p],2];
        If[keep[hh], pl=Append[pl,p]; Print[p]]
    ,{c[2],-a[[2]],a[[2]]},{c[3],-a[[3]],0},{c[4],-a[[4]],a[[4]]}];
    pl
]


PolynomialRoots[p_,opts:OptionsPattern[]] := Module[
    {x = PolynomialVariable[p]},
    Sort[x/.NSolve[p == 0, x, opts], Abs[#2] < Abs[#1] &]
]
(* PolynomialRoots inherits the options for NSolve *)
Options[PolynomialRoots] = Options[NSolve]


PerronRoot[p_,opts:OptionsPattern[]] := First[PolynomialRoots[p,opts]]
(* PerronRoot inherits the options for NSolve *)
Options[PerronRoot] = Options[NSolve]


MahlerMeasure[p_,opts:OptionsPattern[]] := Module[
    {x = PolynomialVariable[p]},
    Times @@ Select[Abs/@PolynomialRoots[p,opts],#>1&]
]
(* MahlerMeasure inherits the options for NSolve *)
Options[MahlerMeasure] = Options[NSolve]


IrreducibleMatrixQ[M_List] := Module[{n = Length[M], powmax},
    (* See Ham and Song paper (2007), p. 172; Seneta 73. Theorem 2.8 *)
    powmax = n^2 - 2n + 2;
    (* I'm not sure it's not better to take the Abs value of the
       elements.  Depends on what we want. *)
    Fold[#1 && #2 != 0 &, True, Flatten[MatrixPower[M, powmax]]]
]


StrataList[g_Integer,n_Integer:0] := Module[{},
    (* If n>0, we need at least one singularity per boundary
       component.  Note that a p-pronged singularity on the boundary
       is really a (p-2)-prong when the boundary components are shrunk
       to punctures. *)
    Message[PseudoAnosov::nottested];
    (* Need to separate out the singularities on punctures and those
       away from them.  The ones on punctures should get a -2. To
       support a pA, 1-prongs must be on the boundary. *)
    Sort /@ Select[IntegerPartitions[2(2g+n-2)],Length[#]>=n&]
    (* Convert n 3-prongs to 1-prongs on boundaries. *)
    (* If we don't have enough 3-prongs, convert 4-prongs to 2-prongs
       on boundaries. *)
    (* The other cases I care much less about. *)
    (* For g>0, possible to have orientable foliations with punctures,
       for instance by putting 2-prongs on each boundary. Use a
       convention to disinguish singularities on boundaries.  For
       instance, list a stratum as a pair of lists: the first list
       contains the order of the n singularities on the punctures; the
       second the other singularities. Change "boundaries" to
       "punctures" later on. *)
]
(* Maybe it's better to do this from the partition function directly, by subtracting 2 afterwards. *)

OrientableStrataList[g_Integer] := 2 IntegerPartitions[2g-2]


StratumToGenus[s_List] := (Plus @@ s + 4)/4


(*

  The pseudo-Anosov polynomial tests
  
  These are based on the Lefschetz numbers calculated from the
  polynomial.
  
  ToDo:

   - Give more detailed reason, customized for each test.
   - Find more harmonious names: StratumTest?  PolynomialStratumQ?  pAStratumQ?

 *)


(* Test for the Perron root *)
pseudoAnosovPerronRootQ[p_,lmax_:0,opts:OptionsPattern[]] := Module[
    {prl, pr, degen, testdegen},
    If[PolynomialDegree[p] < 2, Return[False]];
    degen = OptionValue[EqualityTolerance];
    testdegen[diff_] := Module[{},
        If[diff < degen, Return[False]];
        If[degen/2 < diff < degen, Message[PseudoAnosov::tooclose]];
        Return[True]
    ];
    (* Take the first two roots (presorted by magnitude) *)
    prl = Take[PolynomialRoots[p,opts], 2];
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
    {WorkingPrecision -> $MachinePrecision, EqualityTolerance -> 10^-8}


(* General helper functions and definitions for the tests *)

SingularityToString[k_Integer, m_Integer] :=
    ToString[k] <> "^" <> ToString[m]


LefschetzToString[L_Integer, n_Integer] :=
    "L[" <> ToString[n] <> "]=" <> ToString[L]


LCMToString[k_Integer] := "LCM(partitions of " <> ToString[k] <> ")"


AllPossibilities[L__List, OptionsPattern[]] := Module[{f, pos},
    (* The function f helps generate all possible combinations of
       elements of given lists *)
    f[l_, a_] := Flatten[Table[
         Append[l[[i]], a[[j]]], {i, Length[l]}, {j, Length[a]}], 1];
    (* Use Fold to apply f to each list, giving all possibilities *)
    pos = Fold[f,{{}},List[L]];
    (* If Ordered is False, order doesn't matter *)
    If[OptionValue[Ordered], pos, Union[Sort/@pos]]
]
Options[AllPossibilities] = {Ordered -> True}


Allowable = "Allowable"


(*
   Lefschetz tests for positive Perron root
*)

(* Test whether there are enough singularities on a stratum to support
   this pseudo-Anosov (when Perron root is positive). *)
LefschetzMinimumSingularitiesQ[s_List,L_List] := Module[{idx},
    idx = Flatten[Position[# >= 0 & /@ (Length[s] - L), False]];
    If[idx == {},
        Return[Allowable]
    ,
        idx = First[idx];
        Return["Not enough singularities for " <>
               LefschetzToString[L[[idx]],idx]]
    ]
]


LefschetzSingularityPermutationsAQ[s_List,L_List, OptionsPattern[]] := Module[
    {d, m, Nn = Length[s], pow, idx},
    (* Group singularities by multiplicity *)
    d = #[[1]]/2& /@ Tally[s];
    m = #[[2]]& /@ Tally[s];
    (* Check if the formula is satisfied for each singularity type,
       logical-And the results. *)
    pow = (#[[2]])! (#[[1]]/2+1) & /@ Tally[s];
    (* Don't go past the end of L for this test (no complaints --
       factorial is just too big).  Also check for the optional
       IterateTest criterion (used to select even iterates, say). *)
    idx = Pick[Range[Length[d]],
        (# <= Length[L] && OptionValue[IterateTest][#]) & /@ pow];
    If[idx != {},
        idx = Pick[idx, Not/@(L[[pow[[#]]]] <= Nn-2m[[#]](d[[#]]+1) & /@ idx)];
        If[idx == {},
            Return[Allowable]
        ,
            idx = First[idx];
            Return[SingularityToString[2d[[idx]],m[[idx]]] <>
                " incompatible with " <>
                LefschetzToString[L[[pow[[idx]]]],pow[[idx]]]]
        ]
    ];
    Return[Allowable]
]
Options[LefschetzSingularityPermutationsAQ] = {IterateTest -> (True &)}


LefschetzSingularityPermutationsBQ[s_List,L_List] := Module[
    {k, m, Nn = Length[s], idx},
    (* Group singularities by multiplicity *)
    k = #[[1]]& /@ Tally[s];
    m = #[[2]]& /@ Tally[s];
    (* Check if the formula is satisfied for each singularity type,
       logical-And the results. *)
    idx = Flatten[Position[Table[
            LefschetzSingularityPermutationsBQ1[k[[j]]/2,m[[j]],Nn,L]
        ,{j,Length[k]}], False]];
    If[idx == {},
        Return[Allowable]
    ,
        idx = First[idx];
        Return[SingularityToString[k[[idx]],m[[idx]]] <>
            " incompatible with " <> LefschetzToString[L[[idx]],idx] <>
            " and " <> LCMToString[m[[idx]]]]
    ]
]
(* Private helper function for LefschetzSingularityPermutationsBQ. *)
LefschetzSingularityPermutationsBQ1[d_Integer,m_Integer,Nn_,L_] := Module[
    (* Compute the (unique) LCMs of integer partitions of m *)
    {lcm = Union[LCM @@ # & /@ IntegerPartitions[m]]},
    (* Test for each LCM in the list and Or the result, since at least
       one of them has to be true. *)
    If[Max[lcm](d+1) > Length[L], Throw[Max[lcm](d+1), oor]];
    Or @@ (L[[#(d+1)]] <= Nn - 2m(d+1) & /@ lcm)
]


LefschetzPureStratumAQ[s_List,L_List] := Module[
    {d, m, t = Tally[s]},
    (* If there is more than one singularity type, return Allowable. *)
    If[Length[t] > 1, Return[Allowable]];
    (* Need positive Lefschetz number L[[1]] *)
    If[L[[1]] <= 0, Return[Allowable]];
    d = t[[1,1]]/2;
    m = t[[1,2]];
    If[(d+1) > Length[L], Throw[d+1, oor]];
    If[L[[d+1]] <= m - 2 (d+1) L[[1]],
        Return[Allowable]
    ,
        Return[SingularityToString[2d,m] <>
            " incompatible with " <> LefschetzToString[L[[1]],1] <> " and " <>
            LefschetzToString[L[[d+1]],d+1]]
    ]
]


LefschetzPureStratumBQ[s_List,L_List] := Module[
    {d, m, t = Tally[s], k},
    (* If there is more than one singularity type, return Allowable. *)
    If[Length[t] > 1, Return[Allowable]];
    d = t[[1,1]]/2;
    m = t[[1,2]];
    (* Need positive Lefschetz number L[[1]] and enough singularities *)
    If[m-L[[1]] < 1 || L[[1]] <= 0, Return[Allowable]];
    k = Union[LCM @@ # & /@ IntegerPartitions[m-L[[1]]]];
    If[Max[k](d+1) > Length[L], Throw[Max[k](d+1), oor]];
    If[Or @@ (L[[#(d+1)]] <= m - 2m(d+1) & /@ k),
        Return[Allowable]
    ,
        Return[SingularityToString[2d,m] <>
            " incompatible with " <> LefschetzToString[L[[1]],1] <>
            " and " <> LCMToString[m-L[[1]]]]
    ]
]


LefschetzAlmostPureStratumAQ[s_List,L_List] := Module[
    {d2, m2, f, Nn = Length[s],
     (* Group and sort strata by increasing multiplicity *)
     t = Sort[Tally[s], #1[[2]] < #2[[2]] &]},
    (* If there aren't exactly two singularity types, return Allowable. *)
    If[Length[t] != 2, Return[Allowable]];
    (* If the first singularity type doesn't have multiplicity 1,
       return Allowable. *)
    If[t[[1,2]] > 1, Return[Allowable]];
    f = L[[1]]-1;
    If[f <= 0, Return[Allowable]];
    (* The repeated singularity is in the second slot, since we sorted. *)
    d2 = t[[2,1]]/2;
    m2 = t[[2,2]];
    If[(d2+1) > Length[L], Throw[d2+1, oor]];
    If[L[[d2+1]] <= Nn - 2f (d2+1),
        Return[Allowable]
    ,
        Return[SingularityToString[2d2,m2] <>
            " incompatible with " <> LefschetzToString[L[[1]],1] <> " and " <>
            LefschetzToString[L[[d2+1]],d2+1]]
    ]
]


LefschetzAlmostPureStratumBQ[s_List,L_List] := Module[
    {d2, m2, f, Nn = Length[s], k,
     (* Group and sort strata by increasing multiplicity *)
     t = Sort[Tally[s], #1[[2]] < #2[[2]] &]},
    (* If there aren't exactly two singularity types, return Allowable. *)
    If[Length[t] != 2, Return[Allowable]];
    (* If the first singularity type doesn't have multiplicity 1,
       return Allowable. *)
    If[t[[1,2]] > 1, Return[Allowable]];
    (* The repeated singularity is in the second slot, since we sorted. *)
    d2 = t[[2,1]]/2;
    m2 = t[[2,2]];
    f = L[[1]]-1;
    If[f <= 0 || m2-f < 1, Return[Allowable]];
    k = Union[LCM @@ # & /@ IntegerPartitions[m2-f]];
    If[Max[k](d2+1) > Length[L], Throw[Max[k](d2+1), oor]];
    If[Or @@ (L[[#(d2+1)]] <= Nn - 2m2 (d2+1) & /@ k),
        Return[Allowable]
    ,
        Return[SingularityToString[2d2,m] <>
            " incompatible with " <> LefschetzToString[L[[1]],1] <>
            " and " <> LCMToString[m2-L[[1]]+1]]
    ]
]


(*
   Lefschetz tests for negative Perron root
*)

LefschetzSingularityPermutationsNegativeQ[s_List,L_List] :=
    (* Call the normal test but only apply to even iterates *)
    LefschetzSingularityPermutationsAQ[s,L, IterateTest->EvenQ]


(*
   Construct regular orbits from Lefschetz numbers
*)

SumOrbits[p_Integer, rpo_List, OptionsPattern[]] := Module[
    {dl = Divisors[p]},
    If[!OptionValue[IncludeLast], dl = Most[dl]];
    Plus @@ (# rpo[[#]] & /@ dl)
]
Options[SumOrbits] = {IncludeLast -> True}


LefschetzRegularOrbits[L_List, OptionsPattern[]] := Module[
    {rpo, def, prs = OptionValue[PerronRootSign]},
    (* If the sign of the Perron root is unspecified, try and guess *)
    If[prs == Automatic, prs = Sign[Times @@ Take[L,-2]]];
    rpo = {-prs L[[1]]};
    If[rpo[[1]] < 0, Message[PseudoAnosov::badLefschetz]; Return[rpo]];
    If[prs < 0 && ((Times @@ Take[L,-2]) > 0),
        Message[PseudoAnosov::neednegativePerron]; Abort[]
    ];
    If[prs > 0 && ((Times @@ Take[L,-2] < 0) || Last[L] > 0),
        Message[PseudoAnosov::needpositivePerron]; Abort[]
    ];
    Do[
        If[prs < 0,
            def = (-1)^(p+1) L[[p]] - SumOrbits[p,rpo,IncludeLast->False]
        ,
            def = -L[[p]] - SumOrbits[p,rpo,IncludeLast->False]
        ];
        AppendTo[rpo, def/p];
        If[!IntegerQ[def/p] || def < 0,
            Message[PseudoAnosov::badLefschetz]; Break[]
        ]
    ,{p, 2, Length[L]}];
    rpo
]
Options[LefschetzRegularOrbits] = {PerronRootSign -> Automatic}


StratumRegularOrbits[s_List,p_, opts:OptionsPattern[]] := Module[
    {L, x = PolynomialVariable[p]},
    L = LefschetzNumbers[p,{OptionValue[MaxLefschetz]}];
    LefschetzStratumRegularOrbits[s,L,opts]
]
Options[StratumRegularOrbits] =
    {MaxLefschetz -> 50, PerronRootSign -> Automatic}


LefschetzStratumRegularOrbits[s_List,L_List, opts:OptionsPattern[]] := Module[
    {Ls, ro},
    Ls = LefschetzNumbersStratum[s,opts];
    Off[PseudoAnosov::badLefschetz];
    ro = LefschetzRegularOrbits[LefschetzCombine[L,-#,Length[L]],opts] & /@ Ls;
    On[PseudoAnosov::badLefschetz];
    (* Eliminate bad orbits: look for a negative or nonintegral last element *)
    ro = Pick[ro, Last[#] > 0 && IntegerQ[Last[#]] & /@ ro];
    If[Length[ro] > 1, Message[PseudoAnosov::manyallowableperms, s]];
    Return[ro]
]
Options[LefschetzStratumRegularOrbits] = Options[StratumRegularOrbits]


StratumRegularOrbitsTestQ[s_List,L_List, opts:OptionsPattern[]] := Module[{},
    Off[PseudoAnosov::manyallowableperms];
    If[LefschetzStratumRegularOrbits[s,L,opts] == {},
        Return["No permutation works"]
    ,
        Return[Allowable]
    ];
    On[PseudoAnosov::manyallowableperms];
]
Options[StratumRegularOrbitsTestQ] = {PerronRootSign -> Automatic}


(* Groups the list l according to an integer partition of Length[l].
   Example: GroupByPartition[{1,2,3,4},{1,3}] = {{1},{2,3,4}} *)
GroupByPartition[l_, part_] := Module[{g = {}},
    Fold[(AppendTo[g, Take[#1, #2]]; Drop[#1, #2]) &, l, part];
    Return[g]
]


LefschetzNumbersSingularity[k_Integer, m_Integer:1, opts:OptionsPattern[]] :=
Module[
    {pr = k/2+1, Pk, Pm, clen, blk, alpo, all},
    (* All possible cyclic permutations of separatrices *)
    Pk = RotateLeft[Range[pr], #] & /@ (Range[pr] - 1);
    (* Keep only the ones that have different cycle lengths *)
    clen[perm_] := Sort[Length/@ToCycles[perm]];
    Pk = Union[Pk, SameTest -> (clen[#1] == clen[#2] &)];
    (* Possible permutations of singularities *)
    Pm = GroupByPartition[Range[m],#] & /@ IntegerPartitions[m];
    blk = Tally /@ ((Length/@#)&/@Pm);
    all = {};
    (* Loop over each possible Pm *)
    Do[
        alpo = {};
        (* Loop over each cyclic block length of Pm
             blk[[i,j,1]] gives the block length
             blk[[i,j,2]] gives the multiplicity (# of blocks of that length)
           Generate a list of combinations of permutations,
             allowing for symmetry do to multiplicity *)
        alpo = ((AllPossibilities[##,Ordered->False]&) @@
                    Table[Pk,{#[[2]]}]& ) /@ blk[[i]];
        (* Now combine the multiplicities together *)
        alpo = AllPossibilities[##]& @@ alpo;
        alpo = {#,Pm[[i]]}& /@ (Flatten[#,1]& /@ alpo);
        all = Join[all,alpo];
    ,{i,Length[Pm]}];
    Union[(LefschetzNumbersSingularityPermutations[##,opts]&) @@ # & /@ all]
]
Options[LefschetzNumbersSingularity] = {PerronRootSign -> -1}


(* Private helper function for LefschetzNumbersSingularity *)
LefschetzNumbersSingularityPermutations[Pk_List, Pm_List:{{1}},
                                        opts:OptionsPattern[]] :=
Module[
    {k = 2Length[Pk]-2, m = Length[Flatten[Pm]], Pmc,
     prs = OptionValue[PerronRootSign]
    },
    Pmc = Length /@ Pm;
    If[Length[Pk] != Length[Pmc], Message[PseudoAnosov::bad]];
    (* Apply SingularityCyclic to each subcycle, with Pk specifying a
       permutation on the separatrices for each subcycle *)
    LefschetzCombine @@
        ((LefschetzNumbersSingularityCyclic[##,opts]&) @@ # & /@
            Transpose[{Pk,Pmc}])
]
Options[LefschetzNumbersSingularityPermutations] =
    Options[LefschetzNumbersSingularity]
(* LefschetzNumbersSingularity::usage = "LefschetzNumbersSingularity[Pk,Pm] returns a list of Lefschetz numbers corresponding to singularity of degree k with m-fold degeneracy.  Pk is a list of permutations on the (k+2)/2 in or outgoing separatrices of the singularities, and Pm a permutation on the m singularities.  Pm must be in cycles form.  Note that Pk must be a power of a cyclic permutation." *)


(* Private helper function for LefschetzNumbersSingularity *)
(* Special case of LefschetzNumbersSingularity for the m singularities
   being permuted cyclically *)
LefschetzNumbersSingularityCyclic[Pk_List, m_Integer:1, OptionsPattern[]] :=
Module[
    {k = 2Length[Pk]-2, period, L, Pm, Pm1, Pk1 = Pk,
     prs = OptionValue[PerronRootSign]
    },
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
            If[EvenQ[i], Pk1 = Permute[Pk1,Pk]];
        ,
            AppendTo[L,0]
        ];
        Pm1 = Permute[Pm1,Pm];
    , {i, period}];
    Return[L]
]
Options[LefschetzNumbersSingularityCyclic] =
    Options[LefschetzNumbersSingularity]


LefschetzNumbersStratum[s_List, opts:OptionsPattern[]] := Module[
    {t, L, f},
    (* Accept s in either explicit or tallied form *)
    If[ListQ[s[[1]]], t = s, t = Tally[s]];
    (* Lists of Lefschetz numbers for each singularity type *)
    L = ((LefschetzNumbersSingularity[##,opts]&) @@ # &) /@ t;
    (* Take all possibilities, and then combine the lists *)
    (LefschetzCombine @@ # &) /@ (AllPossibilities @@ L)
]
Options[LefschetzNumbersStratum] = Options[LefschetzNumbersSingularity];


LefschetzCombine[Ll__List, n_Integer:0] := Module[{L = List[Ll], len},
    If[n == 0, len = LCM @@ (Length /@ L), len = n];
    Plus @@ (PadRight[#,len,#] & /@ L)
]


(*
   Test for everything together
*)

LefschetzNumbersTestQ[s_List,p_, opts:OptionsPattern[]] := Module[
    {t, opts2 = FilterRules[{opts},LefschetzNumbersTestListQ]},
    If[PerronRoot[p] > 0,
        L = LefschetzNumbers[p,{OptionValue[MaxLefschetz]}];
        t = LefschetzNumbersTestPositiveQ[s,L,opts2]
    ,
        (* Generate twice as many Lefschetz numbers, since we need to
           apply the even tests to half the list *)
        L = LefschetzNumbers[p,{2 OptionValue[MaxLefschetz]}];
        t = LefschetzNumbersTestNegativeQ[s,L,opts2]
    ];
    (* If all the tests have failed so far, try the ultimate one, but
       only as a last resort! *)
    If[t == Allowable,
        tests = {StratumRegularOrbitsTestQ};
        t = LefschetzNumbersTestListQ[s,L,tests,opts2,MaxIterate -> 1];
    ];
    If[OptionValue[GiveReasonForRejection],
        Return[{t == Allowable,t}]
    ,
        Return[t == Allowable]
    ]
]
Options[LefschetzNumbersTestQ] =
    {GiveReasonForRejection -> False, MaxIterate -> 5, MaxLefschetz -> 100}


(* Private helper function for LefschetzNumbersTestQ *)
LefschetzNumbersTestListQ[s_List,L_List,tests_List, OptionsPattern[]] :=
Module[
    {Lm, reason = Allowable, dm = 1},
    If[OptionValue[OnlyOddIterates], dm = 2];
    Catch[
        Do[
            Do[
                (* The testing functions throw an exception if there
                   are not enough Lefschetz numbers *)
                Catch[
                    (* Make a list of Lefschetz numbers fpr phi^m *)
                    Lm = L[[#]]& /@ Range[m,Length[L],m];
                    (* If fails once, exit loops to avoid the other tests *)
                    reason = tests[[k]][s,Lm];
                    If[!StringQ[reason], Message[PseudoAnosov::notastring,
                        SymbolName[tests[[k]]]]; Abort[]];
                    If[reason != Allowable,
                        (* For the reason string, append the reason
                           returned by the test to the test's function
                           name. *)
                        reason = StringReplace[SymbolName[tests[[k]]],
                            {"Lefschetz" -> "",
                             RegularExpression["Q$"] -> "",
                             RegularExpression["AQ$"] -> "(a)",
                             RegularExpression["BQ$"] -> "(b)"}]
                             <> ": " <> reason;
                        (* Throw exception to escape both loops *)
                        Throw[k, testfailed];
                    ];
                , oor, Message[PseudoAnosov::moreLefschetz,#1,m] &]
            , {m, 1, OptionValue[MaxIterate], dm}]
        , {k,Length[tests]}];
    , testfailed];
    Return[reason]
]
Options[LefschetzNumbersTestListQ] =
    Append[FilterRules[Options[LefschetzNumbersTestQ], MaxIterate],
           OnlyOddIterates -> False]


(* Private helper function for LefschetzNumbersTestQ *)
LefschetzNumbersTestPositiveQ[s_List,L_, opts:OptionsPattern[]] := Module[
    {tests =
        {LefschetzMinimumSingularitiesQ,
         LefschetzSingularityPermutationsAQ,
         LefschetzSingularityPermutationsBQ,
         LefschetzPureStratumAQ,
         LefschetzPureStratumBQ,
         LefschetzAlmostPureStratumAQ,
         LefschetzAlmostPureStratumBQ}},
    LefschetzNumbersTestListQ[s,L,tests,
        FilterRules[{opts},Options[LefschetzNumbersTestListQ]]]
]
Options[LefschetzNumbersTestPositiveQ] = Options[LefschetzNumbersTestListQ]


(* Private helper function for LefschetzNumbersTestQ *)
LefschetzNumbersTestNegativeQ[s_List,L_, opts:OptionsPattern[]] := Module[
    {L2, t,
    tests = {LefschetzSingularityPermutationsNegativeQ}},
    t = LefschetzNumbersTestListQ[s,L,tests,opts,OnlyOddIterates->True];
    If[t != Allowable, Return[t]];
    (* List of Lefschetz numbers of phi^2 *)
    L2 = L[[#]] & /@ Range[2,Length[L],2];
    (* Do the test with p2, which has positive Perron root. *)
    LefschetzNumbersTestPositiveQ[s,L2,opts]
]
Options[LefschetzNumbersTestNegativeQ] = Options[LefschetzNumbersTestListQ]


End[(* "`Private`" *)]

EndPackage[(* "PseudoAnosov`" *)]
