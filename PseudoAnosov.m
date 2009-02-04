(* PseudoAnosov Mathematica Package *)

BeginPackage["PseudoAnosov`", {"Global`"}]


PseudoAnosov::usage = "Functions for manipulating characteristic polynomials of pseudo-Anosov maps."

TracesPower::usage = "TracesPower[p,x,m], where p(x) is the characteristic polynomial of a matrix M, lists the traces Tr[M^k] for 1 <= k <= m."

LefschetzNumbers::usage = "LefschetzNumbers[p,x,m], where p(x) is the characteristic polynomial of a matrix M, lists the Lefschetz numbers 2-Tr[M^k] for 1 <= k <= m."


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


End[(* "`Private`" *)]

EndPackage[(* "PseudoAnosov`" *)]
