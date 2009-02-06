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


(*
   Error messages and warnings
*)

PseudoAnosov::notminimal = "Not a minimal polynomial."


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


End[(* "`Private`" *)]

EndPackage[(* "PseudoAnosov`" *)]
