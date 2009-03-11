(* PseudoAnosov Mathematica Package *)

BeginPackage["PseudoAnosov`"]


Needs["Combinatorica`"] (* For "ToCycles" *)


(*
   Command usage messages
*)

PseudoAnosov::usage = "Functions for manipulating characteristic polynomials of pseudo-Anosov maps."

PolynomialDegree::usage = "PolynomialDegree[P] returns the degree of the polynomial P(x)."

PolynomialRoots::usage = "PolynomialRoots[P] returns the roots of the polynomial P, sorted in decreasing order of magnitude.";

PerronRoot::usage = "PerronRoot[P] returns the largest root (in magnitute) of the polynomial P.";

pseudoAnosovPerronRootQ::usage = "pseudoAnosovPerronRootQ[P] returns True if the largest root of the polynomial P is nondegenerate (in magnitute) and real.  pseudoAnosovPerronRootQ[P,xmax] also returns False unless the the Perron root is less than xmax (in magnitude).";

MahlerMeasure::usage = "MahlerMeasure[P] returns the Mahler measure of the polynomial P, which is the absolute value of the product of roots outside the unit circle.";

TracesPower::usage = "TracesPower[P,m], where P is the characteristic polynomial of a matrix M, lists the traces Tr[M^k] for 1 <= k <= m."

ReciprocalPolynomial::usage = "ReciprocalPolynomial[x,n] returns a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1.  ReciprocalPolynomial[x,n,c] uses c as the base name for coefficients.  ReciprocalPolynomial[x,n,{a_1,...,a_(n/2)}] uses a list for the coefficient, where (n/2) denotes Floor[n/2].";

ReciprocalPolynomialQ::usage = "ReciprocalPolynomialQ[P] returns true if P is a reciprocal polynomial, i.e. of the form a[0] x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + a[0].";

ReciprocalPolynomialCoefficientBounds::usage = "ReciprocalPolynomialCoefficientBounds[n,t] lists the bound on the magnitude of the coefficients {a[1],...,a[Floor[n/2]]} of a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1, given that its largest eigenvalue is h with t = h + 1/h.";

ReciprocalPolynomialBoundedList::usage = "ReciprocalPolynomialBoundedList[x,n,amax] returns a list of reciprocal polynomials x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1 with coefficients bounded by |a[k]| <= amax[k].  For n even, only one of each polynomials pair P(-x)=P(x) is listed.";

IrreducibleMatrixQ::usage = "IrreducibleMatrixQ[M] returns true if the matrix M is irreducible.";

OrientableStrataList::usage = "OrientableStrataList[g] gives the list of orientable strata for a hyperbolic surface of genus g.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the (even) degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity.  Use Tally/@OrientableStrataList[g] to group by multiplicities.";

StratumToGenus::usage = "StratumToGenus[S] gives the genus of the surface containing a stratum S={k_1,...,k_m}."

LefschetzNumbers::usage = "LefschetzNumbers[P,k], where P is the characteristic polynomial of some matrix M, returns the Lefschetz number 2-Tr[M^k].  LefschetzNumbers[P,{k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for 1 <= k <= k2.  LefschetzNumbers[P,{k1,k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for k1 <= k <= k2."

LefschetzCombine::usage = "LefschetzCombine[L1,L2,...] adds lists of Lefschetz number.  If they are not the same length, then the blocks are repeated to the length of the longest list.\nLefschetzCombine[L1,L2,...,Lm,n] caps the total length at an integer n."

LefschetzNumbersTestQ::usage = "LefschetzNumbersTestQ[S,P] returns true if the polynomial P is compatible with the stratum S.  Possible options are GiveReasonForRejection (default False), MaxIterate (default 1), and MaxLefschetz (default 50)."

StratumOrbits::usage = "StratumOrbits[S,P] returns a list of possible orbit structure (singular and regular periodic orbits) for the polynomial P on stratum S.  Returns an empty list if this proves impossible.\nStratumOrbits[S,L] does the same for a list of Lefschetz numbers L."

StratumOrbitsTable::usage = "StratumOrbitsTable[so] presents the output of StratumOrbits in a table.\nStratumOrbitsTable[so,itmax] displays at most itmax iterates."

(* Test for everything: first call pseudoAnosovPerronRootQ, then
LefschetzNumbersTestQ by strata (for even power). *)
(* pseudoAnosovPolynomialQ::usage = "" *)


(*
   Options
*)

GiveReasonForRejection::usage = "GiveReasonForRejection is an option to LefschetzNumbersTestQ: Set to True to return the reason for rejecting a stratum (default False)."

MaxIterate::usage = "MaxIterate is an option to LefschetzNumbersTestQ: Set to an integer giving the largest power of the map to test (default 1)."

MaxLefschetz::usage = "MaxLefschetz is an option to LefschetzNumbersTestQ and StratumOrbits, specifying how many Lefschetz numbers to compute."

PerronRootSign::usage = "PerronRootSign is an option to StratumOrbits (Lefschetz numbers form) to specify whether the Perron root is positive or negative.  Set to Automatic to try and guess by looking at the last two Lefschetz numbers (default Automatic)."

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
    If[x == {} || Length[x] > 1,
        Message[PseudoAnosov::notapolynomial,p,x]; Abort[]];
    Return[First[x]]
]
polynomialvariable::notapolynomial = "Input argument `1` is not a polynomial in `2`."


PolynomialDegree[p_] := Module[{x = polynomialvariable[p]},
    Length[CoefficientList[Collect[p,x],x]]-1
]


PolynomialRoots[p_,opts:OptionsPattern[]] := Module[
    {x = polynomialvariable[p]},
    Sort[x/.NSolve[p == 0, x, opts], Abs[#2] < Abs[#1] &]
]
(* PolynomialRoots inherits the options for NSolve *)
Options[PolynomialRoots] = Options[NSolve]


PerronRoot[p_,opts:OptionsPattern[]] := First[PolynomialRoots[p,opts]]
(* PerronRoot inherits the options for NSolve *)
Options[PerronRoot] = Options[NSolve]


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


MahlerMeasure[p_,opts:OptionsPattern[]] := Module[
    {x = polynomialvariable[p]},
    Times @@ Select[Abs/@PolynomialRoots[p,opts],#>1&]
]
(* MahlerMeasure inherits the options for NSolve *)
Options[MahlerMeasure] = Options[PolynomialRoots]


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


IrreducibleMatrixQ[M_List] := Module[{n = Length[M], powmax},
    (* See Ham and Song paper (2007), p. 172; Seneta 73. Theorem 2.8 *)
    powmax = n^2 - 2n + 2;
    (* I'm not sure it's not better to take the Abs value of the
       elements.  Depends on what we want. *)
    Fold[#1 && #2 != 0 &, True, Flatten[MatrixPower[M, powmax]]]
]


OrientableStrataList[g_Integer] := 2 IntegerPartitions[2g-2]


StratumToGenus[s_List] := (Plus @@ s + 4)/4


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


Begin["`Tests`"]


Allowable = "Allowable"


singularitytostring[k_Integer, m_Integer] :=
    ToString[k] <> "^" <> ToString[m]


lefschetztostring[L_Integer, n_Integer] :=
    "L[" <> ToString[n] <> "]=" <> ToString[L]


lcmtostring[k_Integer] := "LCM(partitions of " <> ToString[k] <> ")"


(*
   Lefschetz tests for positive Perron root
*)

(* Test whether there are enough singularities on a stratum to support
   this pseudo-Anosov (when Perron root is positive). *)
MinimumSingularitiesQ[s_List,L_List] := Module[{idx},
    idx = Flatten[Position[# >= 0 & /@ (Length[s] - L), False]];
    If[idx == {},
        Return[Allowable]
    ,
        idx = First[idx];
        Return["Not enough singularities for " <>
               lefschetztostring[L[[idx]],idx]]
    ]
]


SingularityPermutationsAQ[s_List,L_List, OptionsPattern[]] := Module[
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
            Return[singularitytostring[2d[[idx]],m[[idx]]] <>
                " incompatible with " <>
                lefschetztostring[L[[pow[[idx]]]],pow[[idx]]]]
        ]
    ];
    Return[Allowable]
]
Options[SingularityPermutationsAQ] = {IterateTest -> (True &)}


SingularityPermutationsBQ[s_List,L_List] := Module[
    {k, m, Nn = Length[s], idx},
    (* Group singularities by multiplicity *)
    k = #[[1]]& /@ Tally[s];
    m = #[[2]]& /@ Tally[s];
    (* Check if the formula is satisfied for each singularity type,
       logical-And the results. *)
    idx = Flatten[Position[Table[
            SingularityPermutationsBQ1[k[[j]]/2,m[[j]],Nn,L]
        ,{j,Length[k]}], False]];
    If[idx == {},
        Return[Allowable]
    ,
        idx = First[idx];
        Return[singularitytostring[k[[idx]],m[[idx]]] <>
            " incompatible with " <> lefschetztostring[L[[idx]],idx] <>
            " and " <> lcmtostring[m[[idx]]]]
    ]
]
(* Private helper function for SingularityPermutationsBQ. *)
SingularityPermutationsBQ1[d_Integer,m_Integer,Nn_,L_] := Module[
    (* Compute the (unique) LCMs of integer partitions of m *)
    {lcm = Union[LCM @@ # & /@ IntegerPartitions[m]]},
    (* Test for each LCM in the list and Or the result, since at least
       one of them has to be true. *)
    If[Max[lcm](d+1) > Length[L], Throw[Max[lcm](d+1), oor]];
    Or @@ (L[[#(d+1)]] <= Nn - 2m(d+1) & /@ lcm)
]


PureStratumAQ[s_List,L_List] := Module[
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
        Return[singularitytostring[2d,m] <>
            " incompatible with " <> lefschetztostring[L[[1]],1] <> " and " <>
            lefschetztostring[L[[d+1]],d+1]]
    ]
]


PureStratumBQ[s_List,L_List] := Module[
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
        Return[singularitytostring[2d,m] <>
            " incompatible with " <> lefschetztostring[L[[1]],1] <>
            " and " <> lcmtostring[m-L[[1]]]]
    ]
]


AlmostPureStratumAQ[s_List,L_List] := Module[
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
        Return[singularitytostring[2d2,m2] <>
            " incompatible with " <> lefschetztostring[L[[1]],1] <> " and " <>
            lefschetztostring[L[[d2+1]],d2+1]]
    ]
]


AlmostPureStratumBQ[s_List,L_List] := Module[
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
        Return[singularitytostring[2d2,m] <>
            " incompatible with " <> lefschetztostring[L[[1]],1] <>
            " and " <> lcmtostring[m2-L[[1]]+1]]
    ]
]


StratumOrbitsTestQ[s_List,L_List, opts:OptionsPattern[]] := Module[
    {prs = OptionValue[PerronRootSign], opts2 = opts},
    (* If the sign of the Perron root is unspecified, try and guess *)
    If[prs == Automatic,
        prs = Sign[Times @@ Take[L,-2]];
        opts2 = FilterRules[{opts},Except[PerronRootSign]];
        opts2 = Sequence @@ Join[opts2,{PerronRootSign -> prs}];
    ];
    Off[StratumOrbits::manyallowableperms];
    If[StratumOrbits[s,L,opts2] == {},
        Return["No permutation works"]
    ,
        Return[Allowable]
    ];
    On[StratumOrbits::manyallowableperms];
]
Options[StratumOrbits] = {PerronRootSign -> Automatic, MaxLefschetz -> 50}
Options[StratumOrbitsTestQ] = Options[StratumOrbits]


(*
   Test for everything together
*)

LefschetzNumbersTestQ[s_List,p_, opts:OptionsPattern[]] := Module[
    {t,
     opts2 = Sequence @@ FilterRules[{opts},Options[TestListQ]],
     opts3 = Sequence @@ FilterRules[{opts},Options[StratumOrbitsTestQ]]},
    If[PerronRoot[p] > 0,
        L = LefschetzNumbers[p,{OptionValue[MaxLefschetz]}];
        t = TestPositiveQ[s,L,opts2]
    ,
        (* Generate twice as many Lefschetz numbers, since we need to
           apply the even tests to half the list *)
        L = LefschetzNumbers[p,{2 OptionValue[MaxLefschetz]}];
        t = TestNegativeQ[s,L,opts2]
    ];
    (* If all the tests have failed so far, try the ultimate one, but
       only as a last resort! *)
    If[t == Allowable,
        tests = {StratumOrbitsTestQ[##,opts3]&};
        t = TestListQ[s,L,tests,opts2,MaxIterate -> 1];
    ];
    If[OptionValue[GiveReasonForRejection],
        Return[{t == Allowable,t}]
    ,
        Return[t == Allowable]
    ]
]
Options[LefschetzNumbersTestQ] =
    {GiveReasonForRejection -> False, MaxIterate -> 1, MaxLefschetz -> 50}


(* Private helper function for LefschetzNumbersTestQ *)
TestListQ[s_List,L_List,tests_List, OptionsPattern[]] :=
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
                    If[!StringQ[reason], Message[TestListQ::notastring,
                        SymbolName[tests[[k]]]]; Abort[]];
                    If[reason != Allowable,
                        (* For the reason string, append the reason
                           returned by the test to the test's function
                           name. *)
                        (* Had to modify this: SymbolName no longer
                           works, since StratumOrbitsTestQ is passed
                           with options as a pure function.  Hence,
                           use ToString, which does not remove the
                           context: we have to do it manually. *)
                        reason = StringReplace[ToString[tests[[k]]],
                            {"PseudoAnosov`Lefschetz`Tests`" -> "",
                             RegularExpression["Q$"] -> "",
                             RegularExpression["AQ$"] -> "(a)",
                             RegularExpression["BQ$"] -> "(b)"}]
                             <> ": " <> reason;
                        (* Throw exception to escape both loops *)
                        Throw[k, testfailed];
                    ];
                , oor, Message[TestListQ::moreLefschetz,#1,m] &]
            , {m, 1, OptionValue[MaxIterate], dm}]
        , {k,Length[tests]}];
    , testfailed];
    Return[reason]
]
Options[TestListQ] =
    Append[FilterRules[Options[LefschetzNumbersTestQ], MaxIterate],
           OnlyOddIterates -> False]
TestListQ::moreLefschetz = "Need at least `1` Lefschetz numbers at `2`th power."
TestListQ::notastring = "Function `1` did not return a proper string."


(* Private helper function for LefschetzNumbersTestQ *)
TestPositiveQ[s_List,L_, opts:OptionsPattern[]] := Module[
    {tests =
        {MinimumSingularitiesQ,
         SingularityPermutationsAQ,
         SingularityPermutationsBQ,
         PureStratumAQ,
         PureStratumBQ,
         AlmostPureStratumAQ,
         AlmostPureStratumBQ}},
    TestListQ[s,L,tests, Sequence @@ FilterRules[{opts},Options[TestListQ]]]
]
Options[TestPositiveQ] = Options[TestListQ]


(* Private helper function for LefschetzNumbersTestQ *)
TestNegativeQ[s_List,L_, opts:OptionsPattern[]] := Module[
    {L2},
    (* There are no separate tests for negative Perron root *)
    (* List of Lefschetz numbers of phi^2 *)
    L2 = L[[#]] & /@ Range[2,Length[L],2];
    (* Do the test with p2, which has positive Perron root *)
    TestPositiveQ[s,L2,opts]
]
Options[TestNegativeQ] = Options[TestListQ]


End[(* "`Tests`" *)]


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
    If[Length[ro] > 1, Message[StratumOrbits::manyallowableperms, s]];
    Join[#[[1]],{RegularOrbits -> #[[2]]}] & /@ ro
]
StratumOrbits::manyallowableperms = "More than one allowable permutation type on stratum `1`."


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
    If[prs < 0 && ((Times @@ Take[L,-2]) > 0),
        Message[Regular::neednegativePerron]; Abort[]
    ];
    If[prs > 0 && ((Times @@ Take[L,-2] < 0) || Last[L] > 0),
        Message[Regular::needpositivePerron]; Abort[]
    ];
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
Regular::needpositivePerron = "Error: This function only applies to positive Perron root."
Regular::neednegativePerron = "Error: This function only applies to negative Perron root."
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


StratumOrbitsTable[so_, itmax_: 20] := Module[
    {itmx = Min[itmax, Length[RegularOrbits /. so]], pr,
     ps, ro, Lro, top, mid, bot, sip, lsip, sep, k, sp, ln},
    (* If the input is a list of solutions rather than one solution,
       call the function on each element *)
    If[Count[Flatten[so],Rule[Polynomial,_]] > 1 || Length[so] == 1,
        Return[StratumOrbitsTable[#,itmax]& /@ so]
    ];
    (* Function to set the math text style *)
    ps[str_] := Style[str,Italic,12];
    (* Regular orbits *)
    ro = Take[RegularOrbits /. so, itmx];
    Lro = sumorbits[#,ro, IncludeLast->True] & /@ Range[itmx];
    pr = PerronRoot[Polynomial/.so];
    If[pr < 0,
        Lro = MapAt[-#&, Lro, {#}&/@Range[2,itmx,2]]
    ,
        Lro = -Lro;
    ];
    top = {
        Framed[Polynomial /. so, Background -> White, FrameStyle -> None],
        "  ",
        Framed[Row[{"Perron root ",PerronRoot[Polynomial/.so]}],
             Background -> White, FrameStyle -> None]
        };
    top = Row[top];

    sip = FromCycles[{#}]& /@ Flatten[SingularitiesPermutation/.so,1];
    lsip = Length /@ (SingularitiesPermutation/.so);
    sep = Flatten[SeparatricesPermutation/.so,1];
    k = Flatten[Table[Table[(First /@ (Stratum/.so))[[i]],{lsip[[i]]}],
        {i,Length[SingularitiesPermutation/.so]}]];
    sp = Superscript @@ # & /@
        Table[{k[[i]],Length[sip[[i]]]}, {i,Length[k]}];
    mid = Row[{
        Framed[Row[{"stratum ",(Superscript @@ #) & /@ (Stratum/.so)}],
             Background -> White, FrameStyle -> None],
        "   ",
        Framed[
        Grid[Join[{ps/@{"sing",Subscript["P","sing"],Subscript["P","sep"]}},
            Transpose[{sp,sip,sep}]]], Background -> White, FrameStyle -> None]
    }];

    sp = Row[{ps["L"],"(",Superscript @@ #,")"}] & /@
        Table[{k[[i]],Length[sip[[i]]]}, {i,Length[k]}];
    ln = Table[SingularCyclic[sep[[i]], Length[sip[[i]]], Sign[pr]],
            {i,Length[k]}];
    ln = PadRight[#,itmx,#]& /@ ln;
    ln = Flatten /@ Transpose[{sp,ln}];

    bot =
        Grid[{
            Join[{ps["n"]}, Range[itmx]],
            Join[{ps["L"]},
                LefschetzNumbers[Polynomial /. so, {itmx}]],
            (* Join[{ps[Subscript["L","so"]]},
                PadRight[SingularitiesLefschetzBlock /. so, itmx,
                    SingularitiesLefschetzBlock /. so]], *)
            Sequence @@ ln,
            Join[{ps[Subscript["L","ro"]]}, Lro],
            Join[{"#ro"}, ro]
        }, Dividers -> {{2 -> True}, {2 -> True, 3 -> True, -2 -> True}}];

    Framed[Column[{top,mid,bot}, Left, 1],
        Background -> Lighter[LightGray, .6], FrameStyle -> None]
]


End[(* "`Orbits`" *)]


End[(* "`Lefschetz`" *)]


EndPackage[(* "PseudoAnosov`" *)]
