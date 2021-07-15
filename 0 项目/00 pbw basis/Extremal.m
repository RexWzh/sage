(* ::Package:: *)

(* ::Section:: *)
(*Functions about extremal graphs*)


BeginPackage["Extremal`"];
Print["Loading functions about extremal graphs..."];


(* ::Subsubsection:: *)
(*Main function*)


SimpleGraphs;
(*IGSubisomorphicQ;*)
SubgraphQ;
GraphF;
GraphK3;
ExtremalGraphs;
path;


Begin["Private`"];


(* ::Subsubsection:: *)
(*functions*)


(*Import functions*)
Needs["IGraphM`"];
SubgraphQ=IGSubisomorphicQ;


path="~/desktop/work_space/1 MMA/0 pkg/simplegraphs";
(*simple graphs data*)
SimpleGraphs[n_]:=Import[path<>"/graph"<>ToString@n<>".g6"];
(*get part of the simple graphs -- too slow*)
SimpleGraphs[n_,a_,b_]:=Import[path<>"/graph"<>ToString@n<>".g6",{"GraphList",Range[a,b]}];
(*get part of the simple graphs*)
SimpleGraphs[n_,k_]:=Import[ToString@StringForm["`1`/graph`2`/graph`2`-`3`.g6",path,n,k]];

(*Graph Subscript[F, k] of 2k+1 vertices (k\[GreaterEqual]1) *)
GraphF[k_?Positive]:=Module[{vertices1,vertices2,edges},
vertices1=Range@k;
vertices2=Range[-1,-k,-1];
edges=(0<->#&)/@Join[vertices1,vertices2];
edges=Join[edges,Thread[TwoWayRule[Range@k,Range[-1,-k,-1]]]];
Graph@edges];

(*Disjoint union of triangules*)
GraphK3[k_?Positive]:=GraphDisjointUnion@@ConstantArray[GraphF@1,k];


(*search for graphs in "graphs"*)
ExtremalGraphs[graphs_,forbidden_,edge_:0]:=Module[{selected,maxedge},
selected=Select[graphs,EdgeCount@#>=edge&&!IGSubisomorphicQ[forbidden,#]&];
If[Length@selected==0,Return@{{},edge}];
maxedge=EdgeCount/@selected//Max;
{Select[selected,EdgeCount@#==maxedge&],maxedge}]


End[];
EndPackage[];
