ReciprocalPolynomialCoefficientBounds::usage = "ReciprocalPolynomialCoefficientBounds[n,t] lists the bound on the magnitude of the coefficients {a[1],...,a[Floor[n/2]]} of a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1, given that its largest eigenvalue is h with t = h + 1/h.";

ReciprocalPolynomialTracesBounds::usage = "ReciprocalPolynomialTracesBounds[n,r] lists the bound on the magnitude of the traces {T[1],...,T[Floor[n/2]]}, where T[k]=Trace[M^k], of the matrix M of a reciprocal polynomial x^n + a[1] x^(n-1) + a[2] x^(n-2) + ... + a[2] x^2 + a[1] x + 1, given that its largest eigenvalue is r."


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
ReciprocalPolynomialBoundedList2[x_,n_,a_List] := Module[
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


(* StrataList::usage = "StrataList[g,n] gives the list of strata for a ...
 hyperbolic surface of genus g with n boundary components (default n=0).  If n>0, we need at least one singularity per boundary component.  Note that a p-pronged singularity on the boundary is really a (p-2)-prong when the boundary components are shrunk to punctures.  Each stratum in the list is of the form {k_1,...,k_m}, where k_i is the degree of each singularity, and the sum over the k_i gives -2(Euler Characteristic).  Use Tally/@OrientableStrataList[g] to group singularities by multiplicity."; *)

StrataList[g_Integer,n_Integer:0] := Module[{},
    (* If n>0, we need at least one singularity per boundary
       component.  Note that a p-pronged singularity on the boundary
       is really a (p-2)-prong when the boundary components are shrunk
       to punctures. *)
    Message[StrataList::nottested];
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
StrataList::nottested = "Warning: This function is not well tested."


(*
   Lefschetz tests for negative Perron root
*)

(****** Buggy ******)
LefschetzSingularityPermutationsNegativeQ[s_List,L_List] :=
    (* Call the normal test but only apply to even iterates *)
    LefschetzSingularityPermutationsAQ[s,L, IterateTest->EvenQ]


PseudoAnosov::notminimal = "Warning: Not a minimal polynomial."
