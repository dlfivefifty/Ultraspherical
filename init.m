(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



BeginPackage["Ultraspherical`","RiemannHilbert`Common`"];
BandedOperator;
ReplaceEntry;
IncreaseLength;
ApplyToRows;
LowerIndex;
LeftIndex;
RightIndex;
IncreaseSize;
GetFill;
SetFill;
GetRowGenerator;
GetFillValue;
Filler;
Givens;
GivensReduce;
DerivativeOperator;
ConversionOperator;
DirichletOperator;
IdentityOperator;
ZeroOperator;
LaplaceOperator;
ListRowReduce;
RightHandSide;
Begin["Private`"];


GetFill[BandedOperator[A_List,fill_List,rowgen_,___],i_]:=If[i>Length[fill],0 First[fill],fill[[i]]];
GetFillValue[BandedOperator[A_List,fill_List,rowgen_],{i_,_}]:=fill[[i]];
GetFillValue[BandedOperator[A_List,fill_List,rowgen_,Filler->fls_],{i_,j_}]:=fill[[i]].fls[j];
SetFill[BandedOperator[A_List,fill_List,rowgen_,opts___],i_,val_]:=BandedOperator[A,ReplacePart[fill,i->val],rowgen,opts];
GetRowGenerator[BandedOperator[A_List,fill_List,rowgen_,___]]:=rowgen;

GetRow[BandedOperator[A_List,fill_List,rowgen_,___],i_]:=If[i>Length[A],rowgen[i],A[[i]]];

BandedOperator/:bnd_BandedOperator[[i_Integer,j_Integer]]:=\[Piecewise]{
 {First[bnd][[i]][[j]], i<=Length[bnd]&&LeftIndex[bnd,i]<=j<=RightIndex[bnd,i]},
 {GetFillValue[bnd,{i,j}], i<=Length[bnd]&&j>RightIndex[bnd,i]},
 {GetRowGenerator[bnd][i][[j]], LeftIndex[bnd,i]<=j<=RightIndex[bnd,i]},
 {0, True}
};
BandedOperator/:bnd_BandedOperator[[i_List,j_Integer]]:=bnd[[#,j]]&/@i;
BandedOperator/:sp_BandedOperator[[i_Integer,jm_Integer;;jM_Integer]]:=Module[{j},
	Table[sp[[i,j]],{j,jm,jM}]];
BandedOperator/:sp_BandedOperator[[im_Integer;;iM_Integer,j_Integer]]:=Module[{i},
	Table[sp[[i,j]],{i,im,iM}]];
BandedOperator/:sp_BandedOperator[[im_Integer;;iM_Integer,jm_Integer;;jM_Integer]]:=Module[{i},
	Table[sp[[i,j]],{i,im,iM},{j,jm,jM}]];

LeftIndex[bnd_BandedOperator,row_]:=Max[GetRow[bnd,row]//FirstIndex,1];
RightIndex[bnd_BandedOperator,row_]:=GetRow[bnd,row]//LastIndex;

LowerIndex[bnd_BandedMatrix,col_]:=Throw["not implemented"];


BandedOperator/:Length[BandedOperator[A_List,___]]:=Length[A];

ToArray[bnd_BandedOperator]:=bnd[[;;Length[bnd],;;RightIndex[bnd,Length[bnd]]]];
BandedOperator/:MatrixForm[bnd_BandedOperator]:=
If[GetRowGenerator[bnd]===Null,
MatrixMap[MatrixForm,bnd[[;;Length[bnd],;;RightIndex[bnd,Length[bnd]]+3]]]//MatrixForm,
MatrixMap[MatrixForm,bnd[[;;Length[bnd]+3,;;RightIndex[bnd,Length[bnd]]+3]]]//MatrixForm];


IncreaseLength[BandedOperator[A_List,fill_List,rowgen_,opts___]]:=BandedOperator[Join[A,{rowgen[Length[A]+1]}],Append[fill,0 Last[fill]],rowgen,opts];


(** Set Length to be at least n **)
IncreaseLength[bnd_BandedOperator,n_]:=If[Length[bnd]>=n,bnd,IncreaseLength[bnd//IncreaseLength,n]];



(** Set Row Length at i to be at least j **)

IncreaseRightIndex[bnd:BandedOperator[A_List,fil_List,rowgen_,opts___],i_,j_]:=Module[{B},
If[j<=RightIndex[bnd,i],
bnd
,
B=A;
B[[i]]=Append[B[[i]],bnd[[i,RightIndex[bnd,i]+1]]];
IncreaseRightIndex[BandedOperator[B,fil,rowgen,opts],i,j]
]
];

IncreaseRowIndexRange[bnd:BandedOperator[A_List,fil_List,rowgen_,opts___],i_,{li_,lr_}]:=Module[{B},
If[lr<=Length[A[[i]]],
B=A;
B[[i]]=IncreaseIndexRange[B[[i]],{li,lr}];
BandedOperator[B,fil,rowgen,opts]
,
B=A;
B[[i]]=Append[B[[i]],bnd[[i,RightIndex[bnd,i]+1]]];
IncreaseRowIndexRange[BandedOperator[B,fil,rowgen,opts],i,{li,lr}]
]
];


RowLength[BandedOperator[A_List,fil_List,rowgen_,opts___],i_]:=A[[i]]//Length;

IncreaseDimension[bnd_BandedOperator,bnd2_BandedOperator]:=Module[{i,Bn},
Bn=IncreaseLength[bnd,Length[bnd2]];
Do[
Bn=IncreaseRowIndexRange[Bn,i,GetRow[bnd2,i]//IndexRange];,
{i,Length[bnd2]}];
Bn];



ReplaceEntry[bnd:BandedOperator[A_List,fil_List,rowgen_,fls:OptionsPattern[]],{i_,j_},p_,opts:OptionsPattern[IncreaseSize->False]]:=Module[{B,nfil},

If[j<LeftIndex[bnd,i],
Throw["left of left entry"]];

If[i>Length[A],
If[!OptionValue[IncreaseSize],
Throw["Replacing entry past size"]];
ReplaceEntry[bnd//IncreaseLength,{i,j},p,opts]
,
If[j<=RightIndex[bnd,i],
B=A;
B[[i]]=ReplaceEntry[B[[i]],j,p];
BandedOperator[B,fil,rowgen,fls]
,
If[!OptionValue[IncreaseSize],
Throw["Replacing entry past size"]];

ReplaceEntry[IncreaseRightIndex[bnd,i,j],{i,j},p,opts]
]
]
];



ApplyToRows[G_,Bn_BandedOperator,{row1_,row2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,rhs},
Bn1=Bn;

Do[
vals=Bn1[[{row1,row2},i]];
If[!NZeroQ[vals[[1]]]&&!NZeroQ[vals[[2]]],
vals=G.vals;
Bn1=ReplaceEntry[Bn1,{row1,i},vals[[1]],IncreaseSize->True];
Bn1=ReplaceEntry[Bn1,{row2,i},vals[[2]],IncreaseSize->True];
];

,{i,LeftIndex[Bn1,row2],RightIndex[Bn1,row2]}];


vals=G.{GetFill[Bn1,row1],GetFill[Bn1,row2]};

Bn1=SetFill[Bn1,row1,vals[[1]]];
Bn1=SetFill[Bn1,row2,vals[[2]]];

If[OptionValue[RightHandSide]===Null,
Bn1,
{Bn1,ApplyToRows[G,{row1,row2},OptionValue[RightHandSide]]}
]
];



BandedOperator/:c_?NumberQ BandedOperator[A_List,fill_List,rowgen_,opts:OptionsPattern[]]:=BandedOperator[c A,c fill, c rowgen[#]&,opts];
BandedOperator/:c_?NumberQ BandedOperator[A_List,fill_List,Null,opts:OptionsPattern[]]:=BandedOperator[c A,c fill, Null,opts];
BandedOperator/: (bndA:BandedOperator[A_List,fill1_List,Null,opts:OptionsPattern[]])+(bndB:BandedOperator[B_List,fill2_List,Null,opts:OptionsPattern[]]):=Module[{bnA,bnB},

bnA=IncreaseDimension[bndA,bndB];
bnB=IncreaseDimension[bndB,bnA];

BandedOperator[First[bnA]+First[bnB],bnA[[3]]+bnB[[3]],  Null,opts]

];

BandedOperator/:( bndA:BandedOperator[A_List,fill1_List,rowgen1_?(!(#===Null)&),opts:OptionsPattern[]])+(bndB:BandedOperator[B_List,fill2_List,rowgen2_,opts:OptionsPattern[]]):=Module[{bnA,bnB,rowlength},
bnA=IncreaseDimension[bndA,bndB];
bnB=IncreaseDimension[bndB,bnA];

rowlength=First[bnA][[-1]]//Length;

BandedOperator[First[bnA]+First[bnB],bnA[[3]]+bnB[[3]],  PadRight[rowgen1[#],rowlength]+PadRight[rowgen2[#],rowlength]&,opts]
];


ApplyToRows[G_,{row1_Integer,row2_Integer},rhsin_List]:=Module[{rhs},
rhs=PadRight[rhsin,Max[row1,row2,Length[rhsin]]];
rhs[[{row1,row2}]]=G.rhs[[{row1,row2}]];
rhs
];


ApplyToRows[G_,{row1_Integer,row2_Integer},{srow1_Integer,srow2_Integer},rhsin_]:=Module[{rhs,rhs1,rhs2},
rhs=PadRight[rhsin,Max[row1,row2,Length[rhsin]],{{}}];


If[row1==row2,
rhs[[row1]]=ApplyToRows[G,{srow1,srow2},rhs[[row1]]];
,
{rhs[[row1]],rhs[[row2]]}=ApplyToRows[G,{srow1,srow2},rhs[[row1]],rhs[[row2]]];
];

rhs
];


ApplyToRows[G_,{row1_Integer,row2_Integer},rhsin1_List,rhsin2_List]:=Module[{rhs1,rhs2},
rhs1=PadRight[rhsin1,Max[row1,Length[rhsin1]]];
rhs2=PadRight[rhsin2,Max[row2,Length[rhsin2]]];

{rhs1[[row1]],rhs2[[row2]]}=G.{rhs1[[row1]],rhs2[[row2]]};


{rhs1,rhs2}
];


ApplyToRows[G_,{row1_Integer,row2_Integer},{srow1_Integer,srow2_Integer},rhsin1_,rhsin2_]:=Module[{rhs1,rhs2,srhs1,srhs2},
rhs1=PadRight[rhsin1,Max[row1,Length[rhsin1]],{{}}];
rhs2=PadRight[rhsin2,Max[row2,Length[rhsin2]],{{}}];

{rhs1[[row1]],rhs2[[row2]]}=ApplyToRows[G,{srow1,srow2},rhs1[[row1]],rhs2[[row2]]];


{rhs1,rhs2}
];

ApplyToRows[G_,{row1_Integer,row2_Integer},{srow1_Integer,srow2_Integer},{ssrow1_Integer,ssrow2_Integer},rhsin_]:=Module[{rhs,rhs1,rhs2},
rhs=rhsin;


If[row1==row2,
rhs[[row1]]=ApplyToRows[G,{srow1,srow2},{ssrow1,ssrow2},rhs[[row1]]];
,
{rhs[[row1]],rhs[[row2]]}=ApplyToRows[G,{srow1,srow2},{ssrow1,ssrow2},rhs[[row1]],rhs[[row2]]];
];

rhs
];



RowZeroQ[0,_]:=True;
RowZeroQ[bnd_BandedOperator,row_]:=If[row<=Length[bnd],
	{First[bnd][[row]],GetFill[bnd,row]},
{GetRowGenerator[bnd][row]}]//Flatten//Abs//Total//NZeroQ;
RowZeroQ[_,_]:=False;
NZeroQ[bnd_BandedOperator]:=(GetRowGenerator[bnd]===Null)&&NZeroQ[{First[bnd],GetFill[bnd,1]}//Flatten//Abs//Total];


ApplyToRows[G_,Bn_BandedOperator,Bnn_BandedOperator,{row1_,row2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,Bn2,i,rhs1,rhs2},
Bn1=Bn;
Bn2=Bnn;

Do[
If[!(NZeroQ[Bn1[[row1,i]]]&&NZeroQ[Bn2[[row2,i]]]),
vals=G.{Bn1[[row1,i]],Bn2[[row2,i]]};



If[vals[[1]]!=Bn1[[row1,i]],
Bn1=ReplaceEntry[Bn1,{row1,i},vals[[1]],IncreaseSize->True];
];
If[vals[[2]]!=Bn2[[row2,i]],
Bn2=ReplaceEntry[Bn2,{row2,i},vals[[2]],IncreaseSize->True];
];
];

,{i,Min[LeftIndex[Bn1,row1],LeftIndex[Bn2,row2]],Max[RightIndex[Bn1,row1],RightIndex[Bn2,row2]]}];


vals=G.{GetFill[Bn1,row1],GetFill[Bn2,row2]};
Bn1=SetFill[Bn1,row1,vals[[1]]];
Bn2=SetFill[Bn2,row2,vals[[2]]];


If[OptionValue[RightHandSide]===Null,
{Bn1,Bn2},
{{Bn1,Bn2},ApplyToRows[G,{row1,row2},OptionValue[RightHandSide]]}
]
];

(*This is for operator of operators *)

NullOperatorQ[G_,B1_,B2_,{srow1_,srow2_}]:=(((RowZeroQ[B1,srow1]||NZeroQ[B1])&&(RowZeroQ[B2,srow2]||NZeroQ[B2]))||(NZeroQ[G[[1,2]]]&&(G[[2,2]]~NEqual~1)&&RowZeroQ[B1,srow1]));

ApplyToRows[G_,BDx_BandedOperator,{row1_,row2_},{srow1_,srow2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,B1,B2,i,rhs,rhs1,rhs2},
Bn1=BDx;

Do[

If[row1==row2,
B1=Bn1[[row1,i]];
If[!NullOperatorQ[G,B1,B1,{srow1,srow2}],
B1=ApplyToRows[G,B1,{srow1,srow2}];
Bn1=ReplaceEntry[Bn1,{row1,i},B1,IncreaseSize->True];
];
,
{B1,B2}=Bn1[[{row1,row2},i]];
If[!NullOperatorQ[G,B1,B2,{srow1,srow2}],
{B1,B2}=ApplyToRows[G,B1,B2,{srow1,srow2}];
Bn1=ReplaceEntry[Bn1,{row1,i},B1,IncreaseSize->True];
Bn1=ReplaceEntry[Bn1,{row2,i},B2,IncreaseSize->True];
];
];
,{i,Min[LeftIndex[Bn1,row1],LeftIndex[Bn1,row2]],Max[RightIndex[Bn1,row1],RightIndex[Bn1,row2]]}];


{B1,B2}={GetFill[Bn1,row1],GetFill[Bn1,row2]};
{B1,B2}=ApplyToRows[G,#[[1]],#[[2]],{srow1,srow2}]&/@Thread[{B1,B2}]//Thread;

Bn1=SetFill[Bn1,row1,B1];
Bn1=SetFill[Bn1,row2,B2];

If[OptionValue[RightHandSide]===Null,
Bn1,
{Bn1,ApplyToRows[G,{row1,row2},{srow1,srow2},OptionValue[RightHandSide]]}
]
];




ApplyToRows[G_,BDx_BandedOperator,Bd2_BandedOperator,{row1_,row2_},{srow1_,srow2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,B1,B2,Bn2,i,rhs1,rhs2,srhs1,srhs2},
Bn1=BDx;
Bn2=Bd2;

Do[
{B1,B2}={Bn1[[row1,i]],Bn2[[row2,i]]};
If[!NullOperatorQ[G,B1,B2,{srow1,srow2}],
{B1,B2}=ApplyToRows[G,B1,B2,{srow1,srow2}];
Bn1=ReplaceEntry[Bn1,{row1,i},B1,IncreaseSize->True];
Bn2=ReplaceEntry[Bn2,{row2,i},B2,IncreaseSize->True];
];

,{i,Min[LeftIndex[Bn1,row1],LeftIndex[Bn2,row2]],Max[RightIndex[Bn1,row1],RightIndex[Bn2,row2]]}];

{B1,B2}={GetFill[Bn1,row1],GetFill[Bn2,row2]};
{B1,B2}=ApplyToRows[G,#[[1]],#[[2]],{srow1,srow2}]&/@Thread[{B1,B2}]//Thread;

Bn1=SetFill[Bn1,row1,B1];
Bn2=SetFill[Bn2,row2,B2];

If[OptionValue[RightHandSide]===Null,
{Bn1,Bn2},
{{Bn1,Bn2},ApplyToRows[G,{row1,row2},{srow1,srow2},OptionValue[RightHandSide]]}
]
];

(*This is for list operators of operators *)

ApplyToRows[G_,BL:{__BandedOperator},{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,B1,B2,rhs},
Bn1=BL;

If[row1===row2,
Bn1[[row1]]=ApplyToRows[G,Bn1[[row1]],{srow1,srow2},{ssrow1,ssrow2}];
,
{B1,B2}=Bn1[[{row1,row2}]];
{B1,B2}=ApplyToRows[G,B1,B2,{srow1,srow2},{ssrow1,ssrow2}];
Bn1[[row1]]=B1;
Bn1[[row2]]=B2;
];

If[OptionValue[RightHandSide]===Null,
Bn1
,
{Bn1,ApplyToRows[G,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},OptionValue[RightHandSide]]}
]
];



Givens[Bin_,Binn_,{i_,j_},k_]:=Module[{A,a,bB},
a=Bin[[i,k]];
bB=Binn[[j,k]];
({
 {a, bB},
 {-bB, a}
})/Norm[{a,bB}]
];
Givens[Bn_,{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},col_,scol_]:=
Givens[Bn[[row1]][[srow1,col]],Bn[[row2]][[srow2,col]],{ssrow1,ssrow2},scol]//N;


GivensReduce[BDx_,{row1_,row2_},opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx,BDx,{row1,row2},row1]//N;
ApplyToRows[G,BDx,{row1,row2},opts]
];
GivensReduce[BDx_,{row1_,row2_},{srow1_,srow2_},opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx[[row1,row1]],BDx[[row2,row1]],{srow1,srow2},srow1]//N;
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},opts]
];
GivensReduce[BDx_,{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx[[row1]][[srow1,srow1]],BDx[[row2]][[srow2,srow1]],{ssrow1,ssrow2},ssrow1]//N;
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},opts]
];


GivensReduce[BDx_,{row1_,row2_},col_,opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx,BDx,{row1,row2},col]//N;
ApplyToRows[G,BDx,{row1,row2},opts]
];
GivensReduce[BDx_,{row1_,row2_},{srow1_,srow2_},col_,scol_,opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx[[row1,col]],BDx[[row2,col]],{srow1,srow2},scol]//N;
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},opts]
];
GivensReduce[BDx_,{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},col_,scol_,opts:OptionsPattern[RightHandSide->Null]]:=ApplyToRows[Givens[BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},col,scol],BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},opts];


RowReduceMatrix[Bin_,Binn_,{i_,j_},k_]:=Module[{A,a,bB},
a=Bin[[i,k]];
bB=Binn[[j,k]];
({
 {1, 0},
 {-bB/a, 1}
})
];
RowReduceMatrix[Bn_,{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},col_,scol_]:=
RowReduceMatrix[Bn[[row1]][[srow1,col]],Bn[[row2]][[srow2,col]],{ssrow1,ssrow2},scol]//N;


BandedOperator/:RowReduce[BDx_BandedOperator,{row1_,row2_},opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=RowReduceMatrix[BDx,BDx,{row1,row2},row1]//N;
ApplyToRows[G,BDx,{row1,row2},opts]
];
BandedOperator/:RowReduce[BDx_BandedOperator,{row1_,row2_},{srow1_,srow2_},opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=RowReduceMatrix[BDx[[row1,row1]],BDx[[row2,row1]],{srow1,srow2},srow1]//N;
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},opts]
];
ListRowReduce[BDx_BandedOperator,{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=RowReduceMatrix[BDx[[row1]][[srow1,srow1]],BDx[[row2]][[srow2,srow1]],{ssrow1,ssrow2},ssrow1]//N;
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},opts]
];


BandedOperator/:RowReduce[BDx_BandedOperator,{row1_,row2_},col_Integer,opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=RowReduceMatrix[BDx,BDx,{row1,row2},col]//N;
ApplyToRows[G,BDx,{row1,row2},opts]
];
BandedOperator/:RowReduce[BDx_BandedOperator,{row1_,row2_},{srow1_,srow2_},col_Integer,scol_Integer,opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=RowReduceMatrix[BDx[[row1,col]],BDx[[row2,col]],{srow1,srow2},scol]//N;
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},opts]
];
ListRowReduce[BDx_,{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},col_Integer,scol_Integer,opts:OptionsPattern[RightHandSide->Null]]:=ApplyToRows[RowReduceMatrix[BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},col,scol],BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},opts];


DerivativeOperator[2,Filler->fls_]:=BandedOperator[{ShiftList[{0,0,4},0]},{0 fls[1]},ShiftList[{0,0,2 (#+1)},1-#]&,Filler->fls];
DerivativeOperator[1,Filler->fls_]:=BandedOperator[{ShiftList[{1},-1]},{0 fls[1]},ShiftList[{#},-#]&,Filler->fls];
DerivativeOperator[1]:=DerivativeOperator[1,Filler->({(-1)^(#-1),1}&)];
DerivativeOperator[2]:=DerivativeOperator[2,Filler->({(-1)^(#-1),1}&)];

ConversionOperator[1,Filler->fls_]:=BandedOperator[{ShiftList[{1,0,-1/2},0]},{0 fls[1]},ShiftList[{1/2,0,-1/2},1-#]&,Filler->fls];
ConversionOperator[1]:=ConversionOperator[1,Filler->({(-1)^(#-1),1}&)];
ConversionOperator[2,Filler->fls_]:=BandedOperator[{ShiftList[{1,0,-2/3,0,1/6},0]},{0 fls[1]},ShiftList[{1/(2 #),0,-(1/(2 (#+2)))-1/(2 #),0,1/(2(#+2))},1-#]&,Filler->fls];
ConversionOperator[2]:=ConversionOperator[2,Filler->({(-1)^(#-1),1}&)];

DirichletOperator[-1]:=BandedOperator[{ShiftList[{1},0]},{{1,0}},Null,Filler->({(-1)^(#-1),1}&)];
DirichletOperator[1]:=BandedOperator[{ShiftList[{1},0]},{{0,1}},Null,Filler->({(-1)^(#-1),1}&)];


IdentityOperator[Filler->fls_]:=BandedOperator[{ShiftList[{1,0,0},0]},{0 fls[1]},ShiftList[{1,0,0},1-#]&,Filler->fls];
IdentityOperator[]:=IdentityOperator[Filler->({(-1)^(#-1),1}&)];

ZeroOperator[Filler->fls_]:=BandedOperator[{ShiftList[{0,0,0},0]},1,{0 fls[1]},ShiftList[{0,0,0},1-#]&,Filler->fls];
ZeroOperator[]:=ZeroOperator[Filler->({(-1)^(#-1),1}&)];

ZeroOperator[1,\[Infinity],Filler->fls_]:=BandedOperator[{ShiftList[{0},0]},1,{0 fls[1]},Null,Filler->fls];
ZeroOperator[1,\[Infinity]]:=ZeroOperator[1,\[Infinity],Filler->({(-1)^(#-1),1}&)];






DirichletOperator[-1,All]:=BandedOperator[{{IdentityOperator[]}},1,{{IdentityOperator[],ZeroOperator[]}},Null,Filler->({(-1)^(#-1),1}&)];
DirichletOperator[1,All]:=BandedOperator[{{IdentityOperator[]}},1,{{ZeroOperator[],IdentityOperator[]}},Null,Filler->({(-1)^(#-1),1}&)];
DirichletOperator[All,-1]:=BandedOperator[{{DirichletOperator[-1]}},1,{{ZeroOperator[1,\[Infinity]],ZeroOperator[1,\[Infinity]]}},{DirichletOperator[-1]}&,Filler->({(-1)^(#-1),1}&)];
DirichletOperator[All,1]:=BandedOperator[{{DirichletOperator[1]}},1,{{ZeroOperator[1,\[Infinity]],ZeroOperator[1,\[Infinity]]}},{DirichletOperator[1]}&,Filler->({(-1)^(#-1),1}&)];


DerivativeOperator[0,1]:=BandedOperator[{{ConversionOperator[1]}},0,{ZeroOperator[]},{(#)ConversionOperator[1]}&];
DerivativeOperator[1,0]:=BandedOperator[{{DerivativeOperator[1],ZeroOperator[],-DerivativeOperator[1]/2}},1,{ZeroOperator[]},{DerivativeOperator[1]/2,ZeroOperator[],-DerivativeOperator[1]/2}&];


DerivativeOperator[0,2]:=BandedOperator[{{4 ConversionOperator[2]}},-1,{ZeroOperator[]},{2(#+1)ConversionOperator[2]}&];
DerivativeOperator[2,0]:=BandedOperator[{{DerivativeOperator[2],ZeroOperator[],-2/3 DerivativeOperator[2],ZeroOperator[],DerivativeOperator[2]/6}},1,{ZeroOperator[]},{DerivativeOperator[2]/(2 #),ZeroOperator[],(-(1/(2 (#+2)))-1/(2 #))DerivativeOperator[2],ZeroOperator[],DerivativeOperator[2]/(2(#+2))}&];


LaplaceOperator:=BandedOperator[{{DerivativeOperator[2],ZeroOperator[],BandedOperator[{(-2/3) {0,0,4,0,0}+4 {1,0,-2/3,0,1/6}},1,{{0,0}},(-2/3){0,0,2 (#+1),0,0}+4{1/(2 #),0,-(1/(2 (#+2)))-1/(2 #),0,1/(2(#+2))}&,Filler->({(-1)^(#-1),1}&)],ZeroOperator[],DerivativeOperator[2]/6}},1,{{ZeroOperator[],ZeroOperator[]}},{DerivativeOperator[2]/(2 #),ZeroOperator[],BandedOperator[{(-(1/(2 (#+2)))-1/(2 #)) {0,0,4,0,0}+2(#+1){1,0,-2/3,0,1/6}},1,{{0,0}},Function[rw,(-(1/(2 (#+2)))-1/(2 #)){0,0,2 (rw+1),0,0}+2(#+1){1/(2rw),0,-(1/(2 (rw+2)))-1/(2 rw),0,1/(2(rw+2))}],Filler->({(-1)^(#-1),1}&)],ZeroOperator[],DerivativeOperator[2]/(2(#+2))}&,Filler->({(-1)^(#-1),1}&)];


End[];
EndPackage[];
