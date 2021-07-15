(* ::Package:: *)

BeginPackage["Kostant`"];Print["Loading kostant functions..."];


(* ::Subsubsection::Closed:: *)
(*\:4e3b\:8981\:51fd\:6570*)


(*\:6570\:503c\:8ba1\:7b97*)
Kc;KcMove;KcNormal;
(*\:8bb0\:53f7*)
RR;S;
(*\:5206\:89e3\:51fd\:6570*)
subtract;ConsecutiveQ;RootCQ;KosDec;
(*NextLevel; KosDec\:8ba1\:7b97\:8fc7\:7a0b\:4f7f\:7528\:7684\:51fd\:6570*)
(*\:4e0b\:6807\:8f6c\:5316*)
RootConC;

(*\:57fa\:4e8e\:53cc\:53c2\:9012\:63a8\:7684\:6b63\:6839\:5206\:89e3*)
R;ReverseRoots;
KDa;KDc;
(*\:7b26\:53f7\:6839 \:8f6c \:5b57\:7b26\:4e32*)
RootString;RootsSeperate;
(*\:57fa\:4e8e\:5355\:4f4d\:56fe\:7684\:6b63\:6839\:5206\:89e3*)
DecK;
DecD;

(*\:7b26\:53f7\:6839 \[TwoWayRule] \:7d22\:5f15\:6839*)
Roots2Index;Index2Roots;Index2String;
(*\:5bf9\:79f0\:5207\:53e3*)
SymmetricCuts;
(*\:6b63\:5219\:5207\:53e3*)
RegularCuts;
(*\:5224\:65ad\:662f\:5426\:57fa\:5143*)
StandardIndexQ;
(*\:8ba1\:7b97\:7b49\:4ef7\:56fe-\:548c\:4e0b\:8fb9\:51fd\:6570\:91cd\:590d\:4e86*)
(*EquivalentIndexs;*)
(*SequenceFlips;(*\:5e8f\:5217\:5c40\:90e8\:7ffb\:8f6c\:ff0c\:8fc7\:7a0b\:51fd\:6570*)*)
(*\:7f3a\:4e00\:6839\:7684\:6b63\:6839\:5206\:89e3*)
StandardIndexBasis;StandardRootBasis;
(*DeleteLastCase;(*\:4ece\:540e\:5f80\:524d\:5220\:9664\:ff0c\:8fc7\:7a0b\:51fd\:6570*)*)

(*index to standard form*)
(*LocalFlip;*)
IndexStandardize;
(*indexs \[TwoWayRule] fold indexs \[Rule] string*)
Index2Fold;Fold2Index;
Fold2String;(*Block2String;*)

(*Map inverse*)
ActionOnRoots;
Map2Matrix;MapInverse;MapInverseQ;Matrix2Map;
MapInverseViaMatrix;
InverseActionMap;

(*solution of basis relations*)
BasisRelations;
RelationFilter;AvailableRelations;
SolveRelations;SolveRelationsWithSteps;GraphRelations;

(*rule of games*)
IndexRelations;
FoldOrdering;
BasisOrdering;ReverseOrdering;


Begin["Private`"]


(* ::Subsubsection::Closed:: *)
(*Part 1- Values of Kostant partition*)


(*C\:578b\:8fde\:7eed\:6839\:7684kostant\:51fd\:6570*)
Kc::usage="\:8fde\:7eed\:6839\:7684\:5206\:89e3\:6570\:76ee\:ff0c\:5141\:8bb8\:65ad\:4e00\:70b9";
Kc[0,0]=1;
Kc[1,1]=3;
Kc[left_,right_]:=
Module[{l=left,r=right},
(*\:4e0d\:59a8\:8bbel\[GreaterEqual]r*)
If[l<r,Kc[r,l],
(*\:8ba1\:7b97\:7ed3\:679c\:4fdd\:5b58*)
Kc[l,r]=
Which[
l>r+1,2^(l-r-1) Kc[r+1,r], 
l==r+1,Kc[l,l]-Kc[r,r],
l==r>1,5 Kc[l-1,r-2]
]]];


(* C\:578b\:51cf\:4e00\:6839\:7684kostant\:51fd\:6570 *)
Kc[1,1,0]=1;
Kc[0,0,0]=0;
Kc[left_,right_,minus_]:=
Module[{l=left,r=right,k=minus},
If[k>Max[l,r],Print["k\:6570\:503c\:5927\:4e8e\:4e24\:7aef\:6570\:503c"];Return[None]];
If[k<0,Print["k\:4e3a\:8d1f\:6570"];Return[None]];
(*\:4e0d\:59a8\:8bbel\[GreaterEqual]r*)
If[l<r,Kc[r,l,k],
(*\:8ba1\:7b97\:7ed3\:679c\:4fdd\:5b58*)
Kc[l,r,k]=
Which[
(**** \:5e73\:51e1\:60c5\:5f62 ****)
(*k\[GreaterEqual]r*)
k==r==0,2^(l-1), 
0<k==l,Kc[k-1,r], 
0<k==r,Kc[l,k-1],
r<k<l,2^(l-k-1) Kc[k-1,r],
(*k<r\[Equal]l*)
1<k+1==r==l,2 Kc[k,k-1],
(**** \:975e\:5e73\:51e1\:60c5\:5f62 ****)
k+1<l==r,Kc[r-1,r,k]+Kc[r-1,r-1,k],
(*k<r<l-1*)
k<r<l-1,2^(l-r-1) Kc[r+1,r,k],
(*k<r\[Equal]l-1*)
k<r==l-1,Sum[Kc[r,i,k],{i,r}]+2^(r-1)
(*,True,Print["\:6f0f\:7b97"];Abort[]*)
]]];

(*C\:578b\:51cf\:4e00\:6839\:7684\:51fd\:6570\:ff0c\:7d22\:5f15\:5e73\:79fb*)
KcMove::usage="C\:578b\:51cf\:4e00\:6839\:7684\:51fd\:6570\:ff0c\:7d22\:5f15\:5e73\:79fb";
KcMove[l_,r_,k_]:=Kc[k+l,k+r,k];


(*C\:578bkostant\:51fd\:6570\:500d\:6570*)
KcNormal::usage="C\:578bkostant\:51fd\:6570\:500d\:6570";
KcNormal[1,1]=2;
KcNormal[2,2]=7;
KcNormal[left_,right_]:=
Module[{l=left,r=right},
(*\:4e0d\:59a8\:8bbel\[GreaterEqual]r*)
If[l<r,KcNormal[r,l],
(*\:8ba1\:7b97\:7ed3\:679c\:4fdd\:5b58*)
KcNormal[l,r]=
Which[
r==0,2^l,
l>r+1>=2,2^(l-r-1) KcNormal[r+1,r], 
l==r+1>=2,KcNormal[l,l]-KcNormal[r,r],
l==r>2,5 KcNormal[l-1,r-2]
]]];


(* ::Subsubsection::Closed:: *)
(*Part 2- Kostant partition |  \[Beta] - \[Alpha] \[Element] \[CapitalPhi] *)


SetAttributes[RR,{Protected,Orderless}]; (*\:6839\:96c6*)
SetAttributes[S,{Protected,Orderless}];(*\:89e3\:96c6*)

(*\:96c6\:5408\:4f5c\:5dee*)
subtract[l1_,l2_]:=Fold[DeleteCases[#1,#2,1,1]&,l1,List@@l2];

(*\:5224\:65ad\:662f\:5426\:4e3a\:7b49\:5dee1\:7684\:6570\:5217*)
ConsecutiveQ[list_]:=list==Range@@MinMax@list;

(*\:5224\:65ad\:662f\:5426\:4e3aC\:578b\:5355\:6839*)
RootCQ[root_RR]:=With[{r=List@@root,max=Max@@root},
ConsecutiveQ[r]||(SubsetQ[r,Range[0,max]]&&(r~subtract~Range@max//ConsecutiveQ))]


(*\:4e0b\:5c42\:89e3\:51fd\:6570
1.\:7b5b\:9009\:548c\:4e3a\:6839\:7684\:4e8c\:5143\:5b50\:96c6
2.\:51cf\:4e24\:6839\:ff0c\:52a0\:4e00\:6839*)
NextLevel::usage="\:7531\:4e00\:89e3\:6c42\:4e0b\:5c42\:89e3\:ff0c\:7b2c\:4e09\:53c2\:53ef\:8bbe\:7f6e\:54c8\:5e0c\:503c"
SetAttributes[NextLevel,HoldRest];(*\:5141\:8bb8\:4fee\:6539\:53c2\:6570*)
NextLevel[root_S]:=root~subtract~#~Append~(Join@@#)&/@ Select[Subsets[root,{2}],RootCQ@(Join@@#)&];
NextLevel[root_S,hash_]:=Select[NextLevel[root],
With[{h=Hash@#},If[MemberQ[hash,h],False,AppendTo[hash,h];True]]&]

(*\:6b63\:6839\:5206\:89e3\:ff0c\:9010\:7ea7\:751f\:6210*)
KosDec::usage="\:8f93\:5165\:4e00\:6839\:ff0c\:6c42\:6b63\:6839\:5206\:89e3"
KosDec[RR[]]=0;(*\:7a7a\:96c6\:60c5\:5f62*)
KosDec[root_RR]:=Module[{r=root,hash={},sol},
sol = RR/@S@@r;(*\:521d\:59cb\:89e3*)
AppendTo[hash,Hash@sol];
NestWhileList[Function[x,NextLevel[#,hash]&/@x//Flatten],{sol},Length@#!=0&][[;;-2]]]

(*\:4e0b\:6807\:8f6c\:5bf9\:5e94\:8fde\:7eed\:6839*)
RootConC[left_,right_]:=(RR@@Range[0,left])~Join~(RR@@Range[right]);
RootConC[left_,right_,minus_]:=RootConC[left,right]//DeleteCases[#,minus,1,1]&;


(* ::Subsubsection::Closed:: *)
(*Part 3- Kostant partition | 2-para-recursion*)


SetAttributes[R,{Protected}]

(*\:6839\:53cd\:5411*)
ReverseRoots=Reverse/@Reverse@#&;


(*\:9012\:5f52\:7b97\:6cd5\:6784\:9020\:7ebf\:6bb5\:89e3*)


(*A\:578b Kostant \:8fde\:7eed\:6839*)
KDa[roots_]:=Switch[Length@roots,
1,{roots},
2,{roots,{Join@@roots}},(*Join\:5408\:5e76\:6839*)
_,(Prepend[#,First@roots]&/@KDa@roots [[2;;]])
~Join~
(KDa[{Join@@roots[[;;2]]}~Join~roots[[3;;]]])];
(*KDa[R@Range[3]]*)


(*C\:578b Kostant \:8fde\:7eed\:6839*)
KDc[0,0,roots_]:={roots}; (*\:5355\:89e3*)

(*\:521d\:59cb\:5316\:7b2c\:4e09\:53c2*)
KDc[left_,right_]:=With[{roots=R/@(Range[left,0,-1]~Join~Range@right)},
KDc[left,right,roots]];

(*\:4e3b\:51fd\:6570*)
KDc[left_,right_,roots_,begin_:{},end_:{}]:=
Module[{sol1,sol2},(*\:51fd\:6570\:4e0d\:68c0\:67e5roots\:ff0c\:89c6\:9700\:8981\:53ef\:6dfb\:52a0*)
If[left==0||right==0,Return[Join[begin,#,end]&/@KDa[roots]]];
(*\:4e0b\:6807\:4e92\:6362\:ff0c\:5e76\:53cd\:5411\:4e24\:6b21*)
If[left<right,Return[ReverseRoots/@KDc[right,left,ReverseRoots@roots,ReverseRoots@end,ReverseRoots@begin]]];
If[
(*left\[Equal]right \:9010\:9636\:9012\:5f52\:ff0c\:5e76\:6dfb\:52a0\:6700\:957f\:6839\:89e3*)
left==right,
sol1 =Table[KDc[left-i, right-i+1, Append[roots[[i+1;;-i-1]],Join@@roots[[-i;;]]], {Join@@roots[[;;i]]}]
,{i,left}]//Flatten[#,1]&;
sol2={{Join@@roots}},
(*left>right \:9010\:9636\:9012\:5f52\:ff0c\:5e76\:6dfb\:52a0\:ff1a\:7ed1\:5b9a\:524dleft+1\:4e2a\:540eA\:578b\:9012\:5f52\:7684\:7ed3\:679c*)
sol1=Table[KDc[left-i,right,roots[[i+1;;]],{Join@@roots[[;;i]]}],{i,left}]//Flatten[#,1]&;
sol2=KDa[{Join@@roots[[;;left+1]]}~Join~roots[[left+2;;]]]];
sol1 = Join[begin,#,end]&/@sol1;
sol2 = Join[begin,#,end]&/@sol2;
sol1~Join~sol2];


(*\:7528\:5b57\:7b26\:4e32\:683c\:5f0f\:67e5\:770b\:6b63\:6839*)
SetAttributes[RootString,Listable]
RootString[R_] := ToString/@List@@R//Riffle[#,"-"]&//StringJoin;

(*\:5b57\:7b26\:4e32\:5f62\:5f0f\:67e5\:770b\:6b63\:6839\:5206\:89e3,\:6b63\:6839\:4e4b\:95f4\:7528|\:9694\:5f00*)
RootsSeperate[roots_]:= roots//RootString//Riffle[#,"|"]&//StringJoin;


(* ::Subsubsection::Closed:: *)
(*Part 4- Kostant partition | two 1-para-recursion*)


(*\:521d\:59cb\:5316\:7b2c\:4e8c\:53c2*)
DecK[n_]:=With[{roots=R/@Range[n,0,-1]~Join~Range@n},
DecK[n,roots]];
DecD[n_]:=With[{roots=R/@Range[n,0,-1]~Join~Range@n},
DecD[n,roots]];

(*\:8bbe\:7f6e\:521d\:503c*)
DecK[0,roots_,begin_:{},end_:{}]:={Join[begin,roots,end]};

DecK[1,roots_,begin_:{},end_:{}]:={
Join[begin,roots,end],
Join[begin,roots[[;;-3]],{Join@@roots[[-2;;]]},end],
Join[begin,{Join@@roots},end]};

DecD[1,roots_,begin_:{},end_:{}]:={Join[begin,roots,end],
Join[begin,roots[[;;-3]],{Join@@roots[[-2;;]]},end],
Join[begin,{Join@@roots[[;;2]]},roots[[3;;]],end],
Join[begin,{Join@@roots},end]};

(*K\:578b\:5206\:89e3*)
DecK[n_,roots_,begin_:{},end_:{}]:=Module[{r1,r2,r3,sol1,sol2,sol3},
r1= roots[[2;;-2]];(*\:4e24\:4fa7\:65ad\:5f00*)
r2=roots[[2;;-3]]~Join~{Join@@roots[[-2;;]]};(*\:65ad\:5f00\:5de6\:4fa7,\:5408\:5e76\:53f3\:4fa7*)
r3=Join[{Join@@roots[[;;2]]},roots[[3;;-3]],{Join@@roots[[-2;;]]}];(*\:4e24\:4fa7\:5408\:5e76*)
sol1 =DecK[n-1,r1,Append[begin,First@roots],Prepend[end,Last@roots]];
sol2=DecD[n-1,r2,Append[begin,First@roots],end];
sol3 = DecK[n-1,r3,begin,end];
Join[sol1,sol2,sol3]];

(*D\:578b\:5206\:89e3*)
DecD[n_,roots_,begin_:{},end_:{}]:=Module[{r1,r2,r3,r4,sol1,sol2,sol3,sol4},
r1= roots[[2;;-2]];(*\:4e24\:4fa7\:65ad\:5f00*)
r2=roots[[2;;-3]]~Join~{Join@@roots[[-2;;]]};(*\:65ad\:5f00\:5de6\:4fa7*)
r3 = {Join@@roots[[;;2]]}~Join~roots[[3;;-2]];(*\:65ad\:5f00\:53f3\:4fa7*)
r4=Join[{Join@@roots[[;;2]]},roots[[3;;-3]],{Join@@roots[[-2;;]]}];(*\:4e24\:4fa7\:5408\:5e76*)
sol1=DecK[n-1,r1,Append[begin,First@roots],Prepend[end,Last@roots]];
sol2=DecD[n-1,r2,Append[begin,First@roots],end];
sol3=DecD[n-1,r3,begin,Prepend[end,Last@roots]];
sol4=DecD[n-1,r4,begin,end];
Join[sol1,sol2,sol3,sol4]];


(* ::Subsubsection::Closed:: *)
(*Part 5- Regular cuts*)


(*\:6570\:636e\:683c\:5f0f\:76f8\:5173*)

(*\:7b26\:53f7\:6839 \[Rule] \:7d22\:5f15\:6839*)
Roots2Index[roots_]:=Module[{ori,res},
(*\:5b9a\:4f4d\:5305\:542b\:5355\:68390\:7684\:6b63\:6839*)
ori=Position[roots,_?(MemberQ[#,0]&)][[1,1]];
(*\:5185\:4fa7\:6839*)
res=(-Last/@roots[[;;ori-1]])~Join~(Last/@roots[[ori;;-2]]+1);
(*\:52a0\:4e0a\:4e24\:4fa7*)
Join[{-roots[[1,1]]-1},res,{roots[[-1,-1]]+1}]];

(*\:5bf9\:79f0\:5207\:53e3
\:8fd4\:56de\:683c\:5f0f\:ff1aAssociation \:5de6\:53f3\:6b63\:5219\:5207\:53e3\[Rule]\:7d22\:5f15
*)
SymmetricCuts[ind_]:=SymmetricCuts[Sort@-Select[ind,Negative],Select[ind,Positive]];
SymmetricCuts[left_,right_]:=Module[{cuts,neglen,linds,rinds},
cuts=left\[Intersection]right;
neglen=Length@left;
linds = Position[left,#][[1,1]]&/@cuts;
rinds = Position[right,#][[1,1]]&/@cuts;
Association@MapThread[#1->{neglen+1-#2,neglen+#3}&,{cuts,linds,rinds}]];

(*\:6b63\:5219\:5207\:53e3*)
RegularCuts[ind_]:=RegularCuts[Sort@-Select[ind,Negative],Select[ind,Positive]];
RegularCuts[left_,right_]:=Module[{cuts,linds,rinds,rcuts=<||>,lcuts=<||>,f,neglen},
(*\:5207\:53e3\:53ca\:7d22\:5f15*)
cuts=(left\[Intersection]right);
neglen=Length@left;
linds = Position[left,#][[1,1]]&/@cuts;
rinds = Position[right,#][[1,1]]&/@cuts;
(*\:5206\:914d\:51fd\:6570*)
f[l_,r_,x_]:=Which[
l==r==1,Nothing,(*\:5185\:90e8\:6ca1\:6709\:5207\:53e3*)
l==1,AppendTo[rcuts,x->{neglen-l+1,neglen+r}],(*\:5185\:90e8\:5207\:53e3\:9760\:53f3*)
r==1,AppendTo[lcuts,x->{neglen-l+1,neglen+r}],(*\:5185\:90e8\:5207\:53e3\:9760\:5de6*)
left[[l-1]]<right[[r-1]],AppendTo[rcuts,x->{neglen-l+1,neglen+r}],(*\:5185\:90e8\:5207\:53e3\:9760\:53f3*)
left[[l-1]]>right[[r-1]],AppendTo[lcuts,x->{neglen-l+1,neglen+r}](*\:5185\:90e8\:5207\:53e3\:9760\:5de6*)];
MapThread[f,{linds,rinds,cuts}];
{lcuts,rcuts}];

(*\:5224\:65ad\:662f\:5426\:6807\:51c6\:57fa\:5143\:ff08\:53f3\:6b63\:5219\:5207\:6570\:4e3a0\:ff09*)
StandardIndexQ[ind_]:=(RegularCuts@ind//Last//Length)==0;

(*\:7d22\:5f15\:5f62\:5f0f\:8f6c\:7b26\:53f7\:6839\:5f62\:5f0f\:ff0c\:7d22\:5f15\:6570\:9700\[GreaterEqual]2*)
Index2Roots[ind_]:=Module[{f},
f[l_,r_]:=Which[
r<0,R@@Range[-l-1,-r,-1],
l<0<r,R@@(Range[-l-1,0,-1]~Join~Range[r-1]),
0<l,R@@Range[l,r-1]];
MapThread[f,{ind[[;;-2]],ind[[2;;]]}]];
(*\:7d22\:5f15\:8f6c\:5b57\:7b26\:4e32*)
Index2String=RootsSeperate@Index2Roots[#]&;


(*\:7d22\:5f15\:65b9\:6cd5\:6c42\:6b63\:6839\:5206\:89e3*)

(*\:8fde\:7eed\:7684\:6b63\:6839\:5206\:89e3-\:7d22\:5f15\:65b9\:6cd5*)
SetAttributes[J,{Orderless,Protected}];
StandardIndexBasis[l_,r_]:=StandardIndexBasis[l,r,J[-l-1,r+1]];
StandardIndexBasis[l_,r_,ind_]:=Module[{basis,f,NextLevel,hash},
(*\:5224\:65ad\:5e76\:8bb0\:5f55\:5f97\:5230\:7684\:5143\:7d20*)
hash={};
f[x_]:=With[{hy=Hash@x},If[StandardIndexQ[List@@x]&&!MemberQ[hash,hy],AppendTo[hash,hy];True,False]];
(*\:9012\:5f52\:51fd\:6570*)
NextLevel[x_]:=Table[
If[MemberQ[x,i]||i==0,Nothing,Append[x,i]],
{i,-l,r}]//Select[f];
(*Nest \:9012\:5f52*)
basis=NestWhileList[Function[x,NextLevel/@x//Flatten],{ind},Length@#!=0&][[;;-2]]//Flatten;
List@@@basis];

(*\:7f3a\:4e00\:6839\:7684\:6b63\:6839\:5206\:89e3*)
StandardIndexBasis[l_,r_,k_Integer]:=Module[{chkQ},
chkQ[ind_,0]:=(MemberQ[ind,1]&&MemberQ[ind,-1]);
chkQ[ind_,x_]:=(MemberQ[ind,x]&&MemberQ[ind,x+1])||(MemberQ[ind,-x]&&MemberQ[ind,-x-1]);
StandardIndexBasis[l,r]//Select[chkQ[#,k]&]];

(*\:7d22\:5f15\:65b9\:6cd5-\:8fd4\:56de\:7b26\:53f7\:6839\:683c\:5f0f*)
StandardRootBasis[l_,r_]:=Index2Roots/@StandardIndexBasis[l,r];
DeleteLastCase[list_,case_]:=Delete[list,Position[list,case]//Flatten//Last];
StandardRootBasis[l_,r_,k_]:=DeleteLastCase[#,R[k]]&/@(Index2Roots/@StandardIndexBasis[l,r,k]);


(*\:6821\:9a8c\:6b63\:5219\:5207\:53e3\:6570\:76ee\:4e0e\:540c\:89e3\:7b49\:4ef7\:7c7b*)

(*\:5c06\:5e8f\:5217\:5c40\:90e8\:7ffb\:8f6c\:ff0c\:540c\:65f6\:8fdb\:884c\:53cd\:53f7*)
SequenceFlips[seq_,flips_]:=Module[{flipseq},
(*\:4e0d\:9700\:8981\:7ffb\:8f6c*)
If[Length@flips==0,Return@{seq}];
(*\:4e2d\:95f4\:7ffb\:8f6c\:ff0c\:5916\:8fb9\:4e0d\:52a8*)
flipseq = Join[seq[[;;flips[[1,1]]]],
-Reverse@seq[[flips[[1,1]]+1;;flips[[1,2]]-1]],
seq[[flips[[1,2]];;]]];
SequenceFlips[seq,flips[[2;;]]]~Join~SequenceFlips[flipseq,flips[[2;;]]]];

(*\:8ba1\:7b97\:7b49\:4ef7\:56fe*)
EquivalentIndexs[ind_]:=Module[{flips},
(*\:6b63\:5219\:5207\:53e3\:7684\:4f4d\:7f6e*)
flips = List@@@RegularCuts@ind//Flatten[#,1]&;
SequenceFlips[ind,flips]
];


(* ::Subsubsection::Closed:: *)
(*Part 6- Folded diagram*)


(*\:5c40\:90e8\:7ffb\:8f6c\:ff0c\:7528\:4e8e\:7d22\:5f15\:6807\:51c6\:5316*)
LocalFlip[ind_,cut_]:=LocalFlip[ind,cut,SymmetricCuts[ind]];
LocalFlip[ind_,cut_,cuts_]:=Module[{l,r,u,v},
{l,r}=cuts[cut];
(*\:5185\:90e8\:6ca1\:6709\:5bf9\:79f0\:5207\:53e3*)
If[cut==Min@Keys[cuts],
Return@Join[ind[[;;l]],-Reverse@ind[[l+1;;r-1]],ind[[r;;]]] ];
(*\:5185\:90e8\:5b58\:5728\:5bf9\:79f0\:5207\:53e3*)
{u,v}=With[{keys=Keys@cuts},cuts[keys[[Position[keys,cut][[1,1]]-1]]]];
(*Join[ind\[LeftDoubleBracket];;l\[RightDoubleBracket],ind\[LeftDoubleBracket]l+1;;u-1\[RightDoubleBracket],ind\[LeftDoubleBracket]u;;v\[RightDoubleBracket],ind\[LeftDoubleBracket]v+1;;r-1\[RightDoubleBracket],ind\[LeftDoubleBracket]r;;\[RightDoubleBracket]]*)
Join[ind[[;;l]],-Reverse@ind[[v+1;;r-1]],ind[[u;;v]],-Reverse@ind[[l+1;;u-1]],ind[[r;;]]]]

(*\:7d22\:5f15\:6807\:51c6\:5316*)
IndexStandardize[ind_]:=IndexStandardize[ind,SymmetricCuts[ind]];
IndexStandardize[ind_,cuts_]:=Module[{rightcuts},
rightcuts=Keys@RegularCuts[ind][[2]];
If[Length@rightcuts==0,Return @ind];
Fold[LocalFlip[#1,#2,cuts]&,ind,rightcuts]]


(*index to fold index*)
Index2Fold[ind_]:=Module[{first,last,f,cuts=List@@SymmetricCuts@ind},
(*if ind has no symmetric cut*)
If[Length@cuts==0,Return[
{S[ReverseSort@-Select[ind,Negative],ReverseSort@Select[ind,Positive]]}]];
(*function that generates a block via two cut-indexs*)
f[l2_,l1_,r1_,r2_]:=S[ReverseSort@-ind[[l2;;l1]],ReverseSort@ind[[r1;;r2]]];
(*the outside layer*)
first=f[1,cuts[[-1,1]],cuts[[-1,2]],-1];
(*the inner layer*)
last=With[{list=ind[[cuts[[1,1]];;cuts[[1,2]]]]},S[ReverseSort@-Select[list,Negative],ReverseSort@Select[list,Positive]]];
Join[{first},Reverse@MapThread[f,{cuts[[2;;,1]],cuts[[;;-2,1]],cuts[[;;-2,2]],cuts[[2;;,2]]}],{last}]]

(*fold index to standardize index*)
Fold2Index[fold_]:=Module[{left,right},
left=Join@@fold[[;;,1]]//DeleteDuplicates;
right=Join@@fold[[;;,2]]//DeleteDuplicates;
(*left side \[GreaterEqual] right side*)
If[Max@left>Max@right,
right~Join~-left,left~Join~-right]//Sort//IndexStandardize]


(*fold blocks to string form*)
Block2String[block_]:={Index2String[-First@block],Index2String[-Last@block]};
(*fold indexs to string form*)
Fold2String[fold_]:=Block2String/@fold//Transpose//MatrixForm;


(* ::Subsubsection::Closed:: *)
(*Part 7- Map inverse*)


(*Replace one element*)
ReplaceOne[list_,old_,new_]:=Append[DeleteCases[list,old,1,1],new];
(*cut out one Subscript[\[Alpha], k] from roots*)
ActionOnRoots[k_Integer][roots_]:=ActionOnRoots[roots,k];
ActionOnRoots[roots_,k_]:=Module[{chkQ,f,newroots={},noncomms,rootsRR},
rootsRR=S@@roots/.R->RR;
chkQ[root_RR]:=root==RR[k]||(MemberQ[root,k]&&RootCQ[DeleteCases[root,k,1,1]]);
noncomms =Select[rootsRR,chkQ[#]&]//Union;
f[root_]:=AppendTo[newroots,
If[Length@root==1,DeleteCases[rootsRR,root,1,1],
ReplaceOne[rootsRR,root,DeleteCases[root,k,1,1]]
]];
f/@noncomms;
newroots];


(* Map matrix of f, where f:domain \[Rule] image*)
Map2Matrix[f_,domain_,image_]:=Module[{line},
Assert[SubsetQ[image,Union@Flatten[f/@domain,1]],{0,"f[domain] is not contained in codomain\n"}];
line[x_]:=With[{img=f[x]},If[MemberQ[img,#],1,0]&/@image];
line/@domain]

(*left inverse of the map f*)
MapInverse[f_,domain_]:=Module[{func},
func[x_]:=Select[domain,MemberQ[f[#],x]&];
func]

(*check that x\[Element]invf(f(x)) for all x\[Element]domain 
i.e. invf is a left inverse of f*)
MapInverseQ[invf_,f_,domain_]:=Module[{chkQ},
chkQ[x_]:=AllTrue[f[x],MemberQ[invf[#],x]&];
AllTrue[domain,chkQ]];

(*construct a map via matrix*)
Matrix2Map[mat_,domain_,image_]:=Module[{func},
func[x_]:=MapThread[If[#1==1,#2,Nothing]&,{mat[[Position[domain,x][[1,1]]]],image}];
func];
(*map inverse via matrix*)
MapInverseViaMatrix[f_,domain_,image_]:=Matrix2Map[Transpose@Map2Matrix[f,domain,image],image,domain];


(*Inverse action map*)
InverseActionMap[l_,r_,k_]:=Module[{basis,basisk,images,mat},
basis =StandardRootBasis[l,r];
basisk=StandardRootBasis[l,r,k];
images = S@@@(basisk/.R->RR);
mat=Map2Matrix[ActionOnRoots[k],basis,images];
Matrix2Map[Transpose@mat,basisk,basis]]


(* ::Subsubsection:: *)
(*Part 8- Solution of basis relations*)


(*relations induced by simple roots
If fromone is True, then the simple root Subscript[\[Alpha], 0] will be omitted*)
BasisRelations[l_,r_]:=BasisRelations[l,r,True];
BasisRelations[l_,r_,fromone_?BooleanQ,sort_:True]:=Module[{basis,basisk,pos,relmap,mat,posbasisk,k,relations={}},
(*Return the (first) position of # in basis*)
pos=Position[basis,#][[1,1]]&;
basis=Index2Roots/@If[sort,SortBy[StandardIndexBasis[l,r],FoldOrdering],StandardIndexBasis[l,r]]//Reverse;
StandardRootBasis[l,r]; 
For[k=Boole@fromone,k<=l,k++,
basisk=StandardRootBasis[l,r,k];
mat=Map2Matrix[ActionOnRoots[k],basis,S@@@basisk/.R->RR];
posbasisk=pos/@(Index2Roots/@StandardIndexBasis[l,r,k]);
relmap=Matrix2Map[Transpose@mat,posbasisk,Range[Length@basis]];
AppendTo[relations,#->DeleteCases[relmap@#,#]&/@posbasisk]];
Join@@relations];


GraphRelations[l_,r_]:=Module[{basis,graphs,ordering,relations,rule,len,graphrelations},
basis=StandardIndexBasis[l,r];
graphs=RootsSeperate/@StandardRootBasis[l,r];
(*\:57fa\:5143\:91cd\:65b0\:6392\:5217*)
ordering=OrderingBy[basis,FoldOrdering]//Reverse;
basis=Part[basis,#]&/@ordering;
graphs=Part[graphs,#]&/@ordering;
(*\:5173\:7cfb\:8f6c\:5316\:56fe*)
relations=BasisRelations[l,r];
rule=#->{#,graphs[[#]]}&/@Range[Length@basis];
graphrelations=relations/.rule;
len=MapThread[ConstantArray,{Range@l,Table[Kc[l,r,i],{i,l}]}]//Flatten;
ReplacePart[graphrelations,{i_,2}:>Prepend[graphrelations[[i,2]],len[[i]]]]];


IndexRelations[l_,r_]:=Module[{basis,relations,rule,len,indexrelations},
(*\:57fa\:5143\:91cd\:65b0\:6392\:5217*)
basis=SortBy[StandardIndexBasis[l,r],FoldOrdering]//Reverse;
(*\:5173\:7cfb\:8f6c\:5316\:56fe*)
relations=BasisRelations[l,r];
rule=#->{#,basis[[#]]}&/@Range[Length@basis];
indexrelations=relations/.rule;
len=MapThread[ConstantArray,{Range@l,Table[Kc[l,r,i],{i,l}]}]//Flatten;
ReplacePart[indexrelations,{i_,2}:>Prepend[indexrelations[[i,2]],len[[i]]]]];


(*tools for solving the relations*)
(*delete useless relations of the current situation*)
RelationFilter[relations_,keys_]:=Module[{exclusions,chkQ},
chkQ=MemberQ[keys,First@#]&&SubsetQ[keys,Last@#]&;
exclusions=Select[relations,chkQ];
relations~subtract~exclusions]

(*available relations of the current situation*)
AvailableRelations[relations_,keys_]:=Module[{available,chkQ},
chkQ=MemberQ[keys,First@#]&&Length@(Last@#~subtract~keys)==1&;
Select[relations,chkQ]];


(*solve the relations with start value keys*)
SolveRelations[relations_,keys_,sort_:0]:=Module[{filtered,available,sols=keys},
filtered=RelationFilter[relations,keys];
While[Length@(available=AvailableRelations[filtered,sols])!=0,
available=First@(Last@#~subtract~sols)&/@available;
AppendTo[sols,Switch[sort,
1,Max@available,-1,Min@available,_,First@available]];
filtered=RelationFilter[relations,sols]];
sols]

SolveRelationsWithSteps[relations_,keys_,sort_:0]:=Module[{filtered,available,sols=keys,steps={},f},
f=First@(Last@#~subtract~sols)&;
filtered=RelationFilter[relations,keys];
While[Length@(available=AvailableRelations[filtered,sols])!=0,
AppendTo[sols,Switch[sort,
(*maximal one*)
1,With[{avail=Last@SortBy[available,{f,First}]},AppendTo[steps,avail];f@avail],
(*minimal one*)
-1,With[{avail=First@SortBy[available,{f,-First[#]&}]},AppendTo[steps,avail];f@avail],
(*first one*)
_,With[{avail=First@available},AppendTo[steps,avail];f@avail]]];
filtered=RelationFilter[relations,sols]];
{sols,steps}];
(*(*solution with steps*)
SolveRelationsWithSteps[relations_,keys_,sort_:0]:=Module[{filtered,available,sols=keys,steps={},f},
f=First@(Last@#~subtract~sols)&;
filtered=RelationFilter[relations,keys];
While[Length@(available=AvailableRelations[filtered,sols])!=0,
AppendTo[sols,Switch[sort,
(*maximal one*)
1,With[{avail=Last@SortBy[available,f]},AppendTo[steps,avail];f@avail],
(*minimal one*)
-1,With[{avail=First@SortBy[available,f]},AppendTo[steps,avail];f@avail],
(*first one*)
_,With[{avail=First@available},AppendTo[steps,avail];f@avail]]];
filtered=RelationFilter[relations,sols]];
{sols,steps}]*)


(* ::Subsubsection::Closed:: *)
(*Part 9- Game rules*)


FoldOrdering[ind_]:=FoldOrdering[ind,Max[-ind]-1,Max@ind-1];
FoldOrdering[ind_,l_,r_]:=Module[{list},
list=Riffle[Range[r],Range[-1,-r,-1]]~Join~Range[-r-1,-l,-1];
Count[ind,#]&/@list//FromDigits[#,2]&];
BasisOrdering[l_,r_]:=OrderingBy[StandardIndexBasis[l,r],(FoldOrdering[#,l,r])&];
ReverseOrdering[ordering_]:=Position[ordering,#][[1,1]]&/@Range[Length@ordering];


(* ::Subsubsection::Closed:: *)
(*EndPackage*)


End[];
EndPackage[]
