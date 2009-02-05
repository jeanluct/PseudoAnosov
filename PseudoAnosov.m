(* PseudoAnosov Mathematica Package *)

BeginPackage["PseudoAnosov`", {"Global`"}]


PseudoAnosov::usage = "Functions for manipulating characteristic polynomials of pseudo-Anosov maps."

TracesPower::usage = "TracesPower[p,x,m], where p(x) is the characteristic polynomial of a matrix M, lists the traces Tr[M^k] for 1 <= k <= m."

LefschetzNumbers::usage = "LefschetzNumbers[p,x,m], where p(x) is the characteristic polynomial of a matrix M, lists the Lefschetz numbers 2-Tr[M^k] for 1 <= k <= m."

ReciprocalPolynomial::usage = "ReciprocalPolynomial[x,n] returns a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1.  ReciprocalPolynomial[x,n,c] uses c as the base name for coefficients.  ReciprocalPolynomial[x,n,{a_1,...,a_(n/2)}] uses a list for the coefficient, where (n/2) denotes Floor[n/2].";

ReciprocalPolynomialCoefficientBounds::usage = "ReciprocalPolynomialCoefficientBounds[n,t] lists the bound on the magnitude of the coefficients {a[1],...,a[Floor[n/2]]} of a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1, given that its largest eigenvalue is h with t = h + 1/h.";


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


End[(* "`Private`" *)]

EndPackage[(* "PseudoAnosov`" *)]
