(* PseudoAnosov Mathematica Package *)

BeginPackage["PseudoAnosov`", {"Global`"}]


(*
   Command usage messages
*)

PseudoAnosov::usage = "Functions for manipulating characteristic polynomials of pseudo-Anosov maps."

TracesPower::usage = "TracesPower[p,x,m], where p(x) is the characteristic polynomial of a matrix M, lists the traces Tr[M^k] for 1 <= k <= m."

LefschetzNumbers::usage = "LefschetzNumbers[p,x,m], where p(x) is the characteristic polynomial of a matrix M, lists the Lefschetz numbers 2-Tr[M^k] for 1 <= k <= m."

ReciprocalPolynomial::usage = "ReciprocalPolynomial[x,n] returns a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1.  ReciprocalPolynomial[x,n,c] uses c as the base name for coefficients.  ReciprocalPolynomial[x,n,{a_1,...,a_(n/2)}] uses a list for the coefficient, where (n/2) denotes Floor[n/2].";

ReciprocalPolynomialQ::usage = "ReciprocalPolynomialQ[p,x] returns true if p(x) is a reciprocal polynomial, i.e. of the form a[0] x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + a[0].";

ReciprocalPolynomialCoefficientBounds::usage = "ReciprocalPolynomialCoefficientBounds[n,t] lists the bound on the magnitude of the coefficients {a[1],...,a[Floor[n/2]]} of a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1, given that its largest eigenvalue is h with t = h + 1/h.";

ReciprocalPolynomialBoundedList::usage = "ReciprocalPolynomialBoundedList[x,n,amax] returns a list of reciprocal polynomials x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1 with coefficients bounded by |a[k]| <= amax[k].  For n even, only one of each polynomials pair P(-x)=P(x) is listed.";

PolynomialRoots::usage = "PolynomialRoots[p,x] returns the roots of the polynomial p(x), sorted in decreasing order of magnitude.";

PerronRoot::usage = "PerronRoot[p,x] returns the largest root (in magnitute) of the polynomial p(x).";

MahlerMeasure::usage = "MahlerMeasure[p,x] returns the Mahler measure of the polynomial p(x), which is the absolute value of the product of roots outside the unit circle.";

MinimalPolynomialQ::usage = "MinimalPolynomialQ[p] returns true if the polynomial p cannot be factored.";

IrreducibleMatrixQ::usage = "IrreducibleMatrixQ[M] returns true if the matrix M is irreducible.";

StrataList::usage = "StrataList[g,n] gives the list of strata for a hyperbolic surface of genus g with n boundary components (default n=0).  If n>0, we need at least one singularity per boundary component.  Note that a p-pronged singularity on the boundary is really a (p-2)-prong when the boundary components are shrunk to punctures.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity.";

OrientableStrataList::usage = "OrientableStrataList[g] gives the list of orientable strata for a hyperbolic surface of genus g.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the (even) degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity.  Use Tally/@OrientableStrataList[g] to group by multiplicities.";

StratumToGenus::usage = "StratumToGenus[S] gives the genus of the surface containing a stratum S={k_1,...,k_m}."

LefschetzMinimumSingularitiesQ::usage = ""

LefschetzLonelySingularitiesQ::usage = ""

LefschetzSingularityPairsQ::usage = ""

LefschetzNumbersTestQ::usage = ""


(*
   Error messages and warnings
*)

PseudoAnosov::notminimal = "Warning: Not a minimal polynomial."

PseudoAnosov::nottested = "Warning: This function is not well tested."

PseudoAnosov::needpositivePerron = "Error: This test only works for positive Perron root."

PseudoAnosov::neednegativePerron = "Error: This test only works for negative Perron root."


Begin["`Private`"]


TracesPower[p_,x_,mm_:1] := Module[
    (* Polynomial is x^n + c[1]x^(n-1) + ... + c[n-1]x + c[n] *)
    {c = Reverse[CoefficientList[Collect[p,x],x]], n, T},
    (* Make sure leading coefficient is 1, then drop it. *)
    c = Drop[c/c[[1]],1];
    n = Length[c]; (* Degree of polynomial *)
    Table[
        T[k] = -Sum[c[[m]] T[k-m], {m,Min[k-1,n]}]
               - If[k > n, 0, k c[[k]]]
    ,{k,mm}]
]


LefschetzNumbers[p_,x_,mm_:1] := (2 - #) & /@ TracesPower[p,x,mm]


ReciprocalPolynomial[x_,n_,a_List] := Module[{aal},
    aal = Table[a[[m]], {m,n/2}];
    aal = If[OddQ[n], Join[aal,Reverse[aal]], Join[aal,Drop[Reverse[aal],1]]];
    aal = Append[Prepend[aal,1],1];
    x^n + Sum[aal[[k]] x^(k-1), {k,n}]
]


ReciprocalPolynomial[x_,n_,c_:Global`a] :=
    ReciprocalPolynomial[x,n,Table[c[k],{k,n/2}]]


ReciprocalPolynomialQ[p_,x_] := Module[
    {c = CoefficientList[Collect[p,x],x], n},
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
        hh = Abs/@Take[PolynomialRoots[p,x],2];
        If[keep[hh], pl=Append[pl,p]; Print[p]]
    ,{c[1],-a[[1]],-1},{c[2],-a[[2]],a[[2]]},{c[3],-a[[3]],a[[3]]}];
    p = p/.{c[1]->0};
    Do[
        hh = Abs/@Take[PolynomialRoots[p,x],2];
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
        hh = Abs/@Take[PolynomialRoots[p,x],2];
        If[keep[hh], pl=Append[pl,p]; Print[p]]
    ,{c[1],-a[[1]],-1},{c[2],-a[[2]],a[[2]]}
    ,{c[3],-a[[3]],a[[3]]},{c[4],-a[[4]],a[[4]]}];
    p = p/.{c[1]->0};
    Do[
        If[Mod[++cases,10000] == 0, Print[cases]];
        hh = Abs/@Take[PolynomialRoots[p,x],2];
        If[keep[hh], pl=Append[pl,p]; Print[p]]
    ,{c[2],-a[[2]],a[[2]]},{c[3],-a[[3]],0},{c[4],-a[[4]],a[[4]]}];
    pl
]


PolynomialRoots[p_,x_,opts:OptionsPattern[]] :=
    Sort[x/.NSolve[p == 0, x, opts], Abs[#2] < Abs[#1] &]

(* PolynomialRoots inherits the options for NSolve *)
Options[PolynomialRoots] = Options[NSolve]


PerronRoot[p_,x_,opts:OptionsPattern[]] := First[PolynomialRoots[p,x,opts]]

(* PerronRoot inherits the options for NSolve *)
Options[PerronRoot] = Options[NSolve]


MahlerMeasure[p_,x_,opts:OptionsPattern[]] := Module[{},
    If[!MinimalPolynomialQ[p],Message[PseudoAnosov::notminimal]];
    Fold[Times,1,Select[Abs/@PolynomialRoots[p,x,opts],#>1&]]]

(* MahlerMeasure inherits the options for NSolve *)
Options[MahlerMeasure] = Options[NSolve]


MinimalPolynomialQ[p_] := Factor[p] === Expand[p]


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

OrientableStrataList[g_Integer] := Sort /@ 2 IntegerPartitions[2g-2]


StratumToGenus[s_List] := (Fold[Plus[#1,#2]&,0,s] + 4)/4


(*

  The pseudo-Anosov polynomial tests
  
  These are based on the Lefschetz numbers calculated from the
  polynomial.
  
  ToDo:

   - If Perronroot<0, then apply the test to phi^2.
   - Find a way to record the reason.
   - Find more harmonious names: StratumTest?  PolynomialStratumQ?  pAStratumQ?

 *)

(* Test whether there are enough singularities on a stratum to support
   this pseudo-Anosov (when Perron root is positive). *)
LefschetzMinimumSingularitiesQ[s_List,p_,x_] := Module[
    {nmax = 100},
    If[PerronRoot[p,x] < 0,
        Message[PseudoAnosov::needpositivePerron]; Return[]];
    Length[s] >= Max[LefschetzNumbers[p,x,nmax]]
]


(* Test whether a stratum with a "lonely" singularity can support this
   pseudo-Anosov (when the Perron root is positive. *)
LefschetzLonelySingularitiesQ[s_List,p_,x_] := Module[
    {k, Nn = Length[s]},
    If[PerronRoot[p,x] < 0,
        Message[PseudoAnosov::needpositivePerron]; Return[]];
    (* Select the "lonely" singularities *)
    k = First /@ Select[Tally[s], #[[2]] == 1 &];
    (* If there are no lonely singularities, return true. *)
    If[k == {},Return[True]];
    (* Check if the formula is satisfied for each singularity,
       logical-And the results. *)
    Fold[And,True,LefschetzLonelySingularityQ[#/2,Nn,p,x]&/@k]
]

(* Helper function for LefschetzLonelySingularitiesQ. *)
LefschetzLonelySingularityQ[d_Integer,Nn_,p_,x_] :=
    Last[LefschetzNumbers[p,x,d+1]] <= Nn - 2(d+1)


(* Test whether a stratum with a "lonely" singularity can support this
   pseudo-Anosov (when the Perron root is positive. *)
LefschetzSingularityPairsQ[s_List,p_,x_] := Module[
    {k, Nn = Length[s]},
    If[PerronRoot[p,x] < 0,
        Message[PseudoAnosov::needpositivePerron]; Return[]];
    (* Select the singularity pairs *)
    k = First /@ Select[Tally[s], #[[2]] == 2 &];
    (* If there are no singularity pairs, return true. *)
    If[k == {},Return[True]];
    (* Check if the formula is satisfied for each singularity pair,
       logical-And the results. *)
    Fold[And,True,LefschetzSingularityPairQ[#/2,Nn,p,x]&/@k]
]

(* Helper function for LefschetzSingularityPairsQ. *)
LefschetzSingularityPairQ[d_Integer,Nn_,p_,x_] :=
    Last[LefschetzNumbers[p,x,2d+2]] <= Nn - 4(d+1)


(* Test for everything. *)
LefschetzNumbersTestQ[s_,p_,x_] :=
    If[PerronRoot[p,x] > 0,
        LefschetzMinimumSingularitiesQ[s,p,q] &&
        LefschetzLonelySingularitiesQ[s,p,x] &&
        LefschetzSingularityPairsQ[s,p,x]
    ]

End[(* "`Private`" *)]

EndPackage[(* "PseudoAnosov`" *)]
