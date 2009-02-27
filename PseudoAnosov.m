(* PseudoAnosov Mathematica Package *)

BeginPackage["PseudoAnosov`", {"Global`"}]


(*
   Command usage messages
*)

PseudoAnosov::usage = "Functions for manipulating characteristic polynomials of pseudo-Anosov maps."

PolynomialDegree::usage = "PolynomialDegree[p,x] returns the degree of the polynomial p(x)."

TracesPower::usage = "TracesPower[p,x,m], where p(x) is the characteristic polynomial of a matrix M, lists the traces Tr[M^k] for 1 <= k <= m."

LefschetzNumbers::usage = "LefschetzNumbers[p,x,k], where p(x) is the characteristic polynomial of some matrix M, returns the Lefschetz number 2-Tr[M^k].  LefschetzNumbers[p,x,{k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for 1 <= k <= k2.  LefschetzNumbers[p,x,{k1,k2}] returns a list of Lefschetz numbers 2-Tr[M^k] for k1 <= k <= k2."

ReciprocalPolynomial::usage = "ReciprocalPolynomial[x,n] returns a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1.  ReciprocalPolynomial[x,n,c] uses c as the base name for coefficients.  ReciprocalPolynomial[x,n,{a_1,...,a_(n/2)}] uses a list for the coefficient, where (n/2) denotes Floor[n/2].";

ReciprocalPolynomialFromTraces::usage = "ReciprocalPolynomialFromTraces[x,n,T] creates a reciprocal polynomial of degree n from a list of traces of powers of its associated matrix.  ReciprocalPolynomialFromTraces[x,T] creates a polynomial of degree 2 Length[T]."

ReciprocalPolynomialQ::usage = "ReciprocalPolynomialQ[p,x] returns true if p(x) is a reciprocal polynomial, i.e. of the form a[0] x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + a[0].";

ReciprocalPolynomialCoefficientBounds::usage = "ReciprocalPolynomialCoefficientBounds[n,t] lists the bound on the magnitude of the coefficients {a[1],...,a[Floor[n/2]]} of a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1, given that its largest eigenvalue is h with t = h + 1/h.";

ReciprocalPolynomialBoundedList::usage = "ReciprocalPolynomialBoundedList[x,n,amax] returns a list of reciprocal polynomials x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1 with coefficients bounded by |a[k]| <= amax[k].  For n even, only one of each polynomials pair P(-x)=P(x) is listed.";

PolynomialRoots::usage = "PolynomialRoots[p,x] returns the roots of the polynomial p(x), sorted in decreasing order of magnitude.";

PerronRoot::usage = "PerronRoot[p,x] returns the largest root (in magnitute) of the polynomial p(x).";

MahlerMeasure::usage = "MahlerMeasure[p,x] returns the Mahler measure of the polynomial p(x), which is the absolute value of the product of roots outside the unit circle.";

MinimalPolynomialQ::usage = "MinimalPolynomialQ[p] returns true if the polynomial p cannot be factored.";

IrreducibleMatrixQ::usage = "IrreducibleMatrixQ[M] returns true if the matrix M is irreducible.";

(* StrataList::usage = "StrataList[g,n] gives the list of strata for a ...
 hyperbolic surface of genus g with n boundary components (default n=0).  If n>0, we need at least one singularity per boundary component.  Note that a p-pronged singularity on the boundary is really a (p-2)-prong when the boundary components are shrunk to punctures.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity."; *)

OrientableStrataList::usage = "OrientableStrataList[g] gives the list of orientable strata for a hyperbolic surface of genus g.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the (even) degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity.  Use Tally/@OrientableStrataList[g] to group by multiplicities.";

StratumToGenus::usage = "StratumToGenus[S] gives the genus of the surface containing a stratum S={k_1,...,k_m}."

LefschetzNumbersTestQ::usage = ""


(*
   Options
*)

GiveReasonForRejection::usage = "Option to LefschetzNumbersTestQ.  Set to True to return the reason for rejecting a stratum."


(*
   Error messages and warnings
*)

PseudoAnosov::toofewtraces = "Error: list of traces should ne at least n/2, where n is the degree of the polynomial."

PseudoAnosov::notminimal = "Warning: Not a minimal polynomial."

PseudoAnosov::nottested = "Warning: This function is not well tested."

PseudoAnosov::needpositivePerron = "Error: This test only works for positive Perron root."

PseudoAnosov::neednegativePerron = "Error: This test only works for negative Perron root."


Begin["`Private`"]


PolynomialDegree[p_,x_] := Length[CoefficientList[Collect[p,x],x]]-1


TracesPower[p_,x_,mm_List:{10}] := Module[
    (* Polynomial is x^n + c[1]x^(n-1) + ... + c[n-1]x + c[n] *)
    {c = Reverse[CoefficientList[Collect[p,x],x]], n, T, ml},
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


TracesPower[p_,x_,m_Integer:1] := First[TracesPower[p,x,{m,m}]]


LefschetzNumbers[p_,x_,mm_List] := (2 - #) & /@ TracesPower[p,x,mm]


LefschetzNumbers[p_,x_,m_Integer:1] := 2 - TracesPower[p,x,m]


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
    Times @@ Select[Abs/@PolynomialRoots[p,x,opts],#>1&]
]
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

OrientableStrataList[g_Integer] := 2 IntegerPartitions[2g-2]


StratumToGenus[s_List] := (Plus @@ s + 4)/4


(*

  The pseudo-Anosov polynomial tests
  
  These are based on the Lefschetz numbers calculated from the
  polynomial.
  
  ToDo:

   - Include the Perron root as a test: >1, nondegen, etc.
   - Find more harmonious names: StratumTest?  PolynomialStratumQ?  pAStratumQ?

 *)

(* Test whether there are enough singularities on a stratum to support
   this pseudo-Anosov (when Perron root is positive). *)
LefschetzMinimumSingularitiesQ[s_List,p_,x_] := Module[
    {nmax = 100},
    If[PerronRoot[p,x] < 0,
        Message[PseudoAnosov::needpositivePerron]; Return[]];
    Length[s] >= Max[LefschetzNumbers[p,x,{nmax}]]
]


LefschetzSingularityPermutationsQ[s_List,p_,x_] := Module[
    {k, m, Nn = Length[s]},
    If[PerronRoot[p,x] < 0,
        Message[PseudoAnosov::needpositivePerron]; Return[]];
    (* Group singularities by multiplicity *)
    k = #[[1]]& /@ Tally[s];
    m = #[[2]]& /@ Tally[s];
    (* Check if the formula is satisfied for each singularity type,
       logical-And the results. *)
    And @@ Table[
        LefschetzSingularityPermutationsQ1[k[[j]]/2,m[[j]],Nn,p,x]
    ,{j,Length[k]}]
]
(* Private helper function for LefschetzSingularityPermutationsQ. *)
LefschetzSingularityPermutationsQ1[d_Integer,m_Integer,Nn_,p_,x_] := Module[
    (* Compute the (unique) LCMs of integer partitions of m *)
    {lcm = Union[LCM @@ # & /@ IntegerPartitions[m]]},
    (* Test for each LCM in the list and Or the result, since at least
       one of them has to be true. *)
    Or @@ (LefschetzNumbers[p,x,#(d+1)] <= Nn - 2m(d+1) & /@ lcm)
]


LefschetzPureStratumQ[s_List,p_,x_] := Module[
    {d, m, t = Tally[s]},
    If[PerronRoot[p,x] < 0,
        Message[PseudoAnosov::needpositivePerron]; Return[]];
    (* If there is more than one singularity type, return True. *)
    If[Length[t] > 1, Return[True]];
    (* Group singularities by multiplicity *)
    d = t[[1,1]]/2;
    m = t[[1,2]];
    LefschetzNumbers[p,x,d+1] <= m - 2 (d+1) LefschetzNumbers[p,x,1]
]


LefschetzAlmostPureStratumQ[s_List,p_,x_] := Module[
    {d2, m, f, Nn = Length[s],
     (* Group and sort strata by increasing multiplicity *)
     t = Sort[Tally[s], #1[[2]] < #2[[2]] &]},
    If[PerronRoot[p,x] < 0,
        Message[PseudoAnosov::needpositivePerron]; Return[]];
    (* If there aren't exactly two singularity types, return True. *)
    If[Length[t] != 2, Return[True]];
    (* If the first singularity type doesn't have multiplicity 1,
       return True. *)
    If[t[[1,2]] > 1, Return[True]];
    f = LefschetzNumbers[p,x,1] - 1;
    If[f <= 0, Return[True]];
    (* The repeated singularity is in the second slot, since we sorted. *)
    d2 = t[[2,1]]/2;
    m2 = t[[2,2]];
    LefschetzNumbers[p,x,d2+1] <= Nn - 2f (d2+1) 
]


(* Test for everything. *)
LefschetzNumbersTestQ[s_List,p_,x_, opts:OptionsPattern[]] := Module[
    {tests, pass = True, reason = "Allowable", n, tl, p2},
    If[PerronRoot[p,x] > 0,
        tests =
            {LefschetzMinimumSingularitiesQ,
             LefschetzSingularityPermutationsQ,
             LefschetzPureStratumQ,
             LefschetzAlmostPureStratumQ};
        Do[
            (* If False once, exit the loop to avoid the other tests *)
            If[!(pass = tests[[k]][s,p,x]),
                If[OptionValue[GiveReasonForRejection],
                    reason = SymbolName[tests[[k]]]];
                Break[]
            ];
        , {k,Length[tests]}];
        If[OptionValue[GiveReasonForRejection],
            Return[{pass,reason}],
            Return[pass]]
    ,
        (* Compute the polynomial of phi^2 *)
        n = PolynomialDegree[p,x];
        tl = TracesPower[p,x,{n}];
        tl = Table[tl[[2k]],{k,n/2}]; (* Keep only even traces *)
        p2 = ReciprocalPolynomialFromTraces[x,tl];
        (* Do the test with p2, which has positive Perron root. *)
        LefschetzNumbersTestQ[s,p2,x,opts]
    ]
]
Options[LefschetzNumbersTestQ] = {GiveReasonForRejection -> False}


End[(* "`Private`" *)]

EndPackage[(* "PseudoAnosov`" *)]
