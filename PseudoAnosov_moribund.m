ReciprocalPolynomialFromTraces::usage = "ReciprocalPolynomialFromTraces[x,n,T] creates a reciprocal polynomial of degree n from a list of traces of powers of its associated matrix.  ReciprocalPolynomialFromTraces[x,T] creates a polynomial of degree 2 Length[T]."

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


(*
   Lefschetz tests for negative Perron root
*)

(****** Buggy ******)
LefschetzSingularityPermutationsNegativeQ[s_List,L_List] :=
    (* Call the normal test but only apply to even iterates *)
    LefschetzSingularityPermutationsAQ[s,L, IterateTest->EvenQ]
