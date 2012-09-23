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
NeumannOperator;
ZeroOperator;
LaplaceOperator;
ListRowReduce;
RightHandSide;
GetRow;
NullOperatorQ;
RowZeroQ;
MultiplicationMatrix;
FillGenerator;
DropColumn;
DropRow;
AlternatingFiller;
BasisOperator;
Begin["Private`"];


CDerivativeFiller[2]:=Function[k,{\[Piecewise]{
 {0, EvenQ[k]},
 {((k-1)/2)^3, True}
},\[Piecewise]{
 {0, EvenQ[k]},
 {((k-1)/2), True}
},\[Piecewise]{
 {1/4 k (2-3 k+k^2), EvenQ[k]},
 {0, True}
},\[Piecewise]{
 {(k-1), EvenQ[k]},
 {0, True}
}}];
AlternatingFiller:=({\[Piecewise]{
 {-1, EvenQ[#]},
 {1, OddQ[#]}
},1}&);
NeumannFiller:=({\[Piecewise]{
 {1, EvenQ[#]},
 {-1, OddQ[#]}
},1}(#-1)^2&);


GetFill[BandedOperator[A_List,fill_List,rowgen_,Filler->fls_,FillGenerator->fgen_],i_]:=If[i>Length[fill],fgen[i],fill[[i]]];
GetFill[BandedOperator[A_List,fill_List,rowgen_,Filler->fls_],i_]:=If[i>Length[fill],0 First[fill],fill[[i]]];
GetFillValue[bnd:BandedOperator[_List,_List,_],{i_,_}]:=GetFill[bnd,i];
GetFillValue[bnd:BandedOperator[_List,_List,_,Filler->fls_,___],{i_,j_}]:=GetFill[bnd,i].fls[j];
SetFill[BandedOperator[A_List,fill_List,rowgen_,opts___],i_,val_]:=BandedOperator[A,ReplacePart[fill,i->val],rowgen,opts];
GetRowGenerator[BandedOperator[A_List,fill_List,rowgen_,___]]:=rowgen;

GetRow[BandedOperator[A_List,fill_List,rowgen_,___],i_]:=If[i>Length[A],rowgen[i],A[[i]]];

BandedOperator/:bnd_BandedOperator[[i_Integer,j_Integer]]:=\[Piecewise]{
 {First[bnd][[i]][[j]], i<=Length[bnd]&&LeftIndex[bnd,i]<=j<=RightIndex[bnd,i]},
 {GetRowGenerator[bnd][i][[j]], LeftIndex[bnd,i]<=j<=RightIndex[bnd,i]},
 {GetFillValue[bnd,{i,j}], j>RightIndex[bnd,i]},
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
MatrixMap[MatrixForm,bnd[[;;Length[bnd],;;RightIndex[bnd,Length[bnd]]+3]]//PadRight[#,{Length[#],Length[#[[1]]]+1},\[Ellipsis]]&]//MatrixForm,
MatrixMap[MatrixForm,bnd[[;;Length[bnd]+3,;;RightIndex[bnd,Length[bnd]]+3]]//PadRight[#,{Length[#],Length[#[[1]]]+1},\[Ellipsis]]&//PadRight[#,{Length[#]+1,Length[#[[1]]]},\[AliasIndicator]]&]//MatrixForm];


IncreaseLength[bnd:BandedOperator[A_List,fill_List,rowgen_,opts___]]:=BandedOperator[Join[A,{rowgen[Length[A]+1]}],Append[fill,GetFill[bnd,Length[A]+1]],rowgen,opts];


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
If[lr<=LastIndex[A[[i]]],
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


RowZeroQ[0,_]:=True;
RowZeroQ[bnd_BandedOperator,row_]:=If[row<=Length[bnd],
	{ToList[First[bnd][[row]]],GetFill[bnd,row]},
{GetRowGenerator[bnd][row]//ToList}]//Flatten//Abs//Total//NZeroQ;
RowZeroQ[_,_]:=False;
NZeroQ[bnd_BandedOperator]:=(GetRowGenerator[bnd]===Null)&&NZeroQ[{First[bnd],GetFill[bnd,1]}//Flatten//Abs//Total];


ReplaceEntry[bnd:BandedOperator[A_List,fil_List,rowgen_,fls:OptionsPattern[]],{i_,j_},p_,opts:OptionsPattern[IncreaseSize->False]]:=Module[{B,nfil},

Which[
((p~NEqual~bnd[[i,j]])===True) ||(NZeroQ[p]&&NZeroQ[bnd[[i,j]]]),
bnd
,
i>Length[A]&&!OptionValue[IncreaseSize],
Throw["Replacing entry past size"]
,
i>Length[A],
ReplaceEntry[bnd//IncreaseLength,{i,j},p,opts]
,
(True===((NZeroQ[p]&&j==LeftIndex[bnd,i])&&(Length[A[[i]]]>1))&&OptionValue[IncreaseSize]),
(** Shorten list **)
B=A;
B[[i]]=Drop[B[[i]],1];
BandedOperator[B,fil,rowgen,fls]
,
j<=RightIndex[bnd,i],
B=A;
B[[i]]=ReplaceEntry[B[[i]],j,p,opts];
BandedOperator[B,fil,rowgen,fls]
,
!OptionValue[IncreaseSize],
Throw["Replacing entry past size"]
,
True 
,
ReplaceEntry[IncreaseRightIndex[bnd,i,j],{i,j},p,opts]
]
];



ApplyToRows[G_,Bn_BandedOperator,{row1_,row2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,rhs},
Bn1=Bn;

Do[
vals=Bn1[[{row1,row2},i]];

vals=G.vals;
Bn1=ReplaceEntry[Bn1,{row1,i},vals[[1]],IncreaseSize->True];
Bn1=ReplaceEntry[Bn1,{row2,i},vals[[2]],IncreaseSize->True];


,{i,Min[LeftIndex[Bn1,row2],LeftIndex[Bn1,row1]],Max[RightIndex[Bn1,row2],RightIndex[Bn1,row1]]}];


vals=G.{GetFill[Bn1,row1],GetFill[Bn1,row2]};

Bn1=SetFill[Bn1,row1,vals[[1]]];
Bn1=SetFill[Bn1,row2,vals[[2]]];

If[OptionValue[RightHandSide]===Null,
Bn1,
{Bn1,ApplyToRows[G,{row1,row2},OptionValue[RightHandSide]]}
]
];



BandedOperator/:c_?NumberQ BandedOperator[A_List,fill_List,rowgen_,Filler->flr_,FillGenerator->flg_]:=BandedOperator[c A,c fill, c rowgen[#]&,Filler->flr,FillGenerator->(c flg[#]&)];


BandedOperator/:( bndA:BandedOperator[A_List,fill1_List,rowgen1_?(!(#===Null)&),Filler->flr1_,FillGenerator->flg1_])+(bndB:BandedOperator[B_List,fill2_List,rowgen2_,Filler->flr2_,FillGenerator->flg2_]):=Module[{bnA,bnB,rowlength},
bnA=IncreaseDimension[bndA,bndB];
bnB=IncreaseDimension[bndB,bnA];



BandedOperator[First[bnA]+First[bnB],bnA[[2]]+bnB[[2]], IncreaseIndexRange[ rowgen1[#],rowgen2[#]//IndexRange]+IncreaseIndexRange[rowgen2[#],rowgen1[#]//IndexRange]&,Filler->flr1,FillGenerator->(flg1[#]+flg2[#]&)]
];

BandedOperator/:( bndA:BandedOperator[A_List,fill1_List,rowgen1_?(!(#===Null)&),Filler->flr1_,FillGenerator->flg1_])+(bndB:BandedOperator[B_List,fill2_List,rowgen2_,Filler->flr2_]):=Module[{bnA,bnB,rowlength},
bnA=IncreaseDimension[bndA,bndB];
bnB=IncreaseDimension[bndB,bnA];



BandedOperator[First[bnA]+First[bnB],bnA[[2]]+bnB[[2]], Function[row,
Module[{lowi,highi},
{lowi,highi}=Thread[{rowgen1[row]//IndexRange,
rowgen2[row]//IndexRange}]//{Min[#[[1]]],Max[#[[2]]]}&;
ShiftList[bndA[[row,lowi;;highi]]+bndB[[row,lowi;;highi]],1-lowi]
]],Filler->flr1,FillGenerator->flg1]
];



BandedOperator/:c_?NumberQ BandedOperator[A_List,fill_List,rowgen_,opts:OptionsPattern[]]:=BandedOperator[c A,c fill, c rowgen[#]&,opts];
BandedOperator/:c_?NumberQ BandedOperator[A_List,fill_List,Null,opts:OptionsPattern[]]:=BandedOperator[c A,c fill, Null,opts];
BandedOperator/: (bndA:BandedOperator[A_List,fill1_List,Null,opts:OptionsPattern[]])+(bndB:BandedOperator[B_List,fill2_List,Null,opts:OptionsPattern[]]):=Module[{bnA,bnB},

bnA=IncreaseDimension[bndA,bndB];
bnB=IncreaseDimension[bndB,bnA];

BandedOperator[First[bnA]+First[bnB],bnA[[2]]+bnB[[2]],  Null,opts]

];

BandedOperator/:( bndA:BandedOperator[A_List,fill1_List,rowgen1_?(!(#===Null)&),opts:OptionsPattern[]])+(bndB:BandedOperator[B_List,fill2_List,rowgen2_,opts:OptionsPattern[]]):=Module[{bnA,bnB,rowlength},
bnA=IncreaseDimension[bndA,bndB];
bnB=IncreaseDimension[bndB,bnA];



BandedOperator[First[bnA]+First[bnB],bnA[[2]]+bnB[[2]], IncreaseIndexRange[ rowgen1[#],rowgen2[#]//IndexRange]+IncreaseIndexRange[rowgen2[#],rowgen1[#]//IndexRange]&,opts]
];


SparseRule[sl_ShiftList,rw_]:={rw,#}->sl[[#]]&/@(Range@@IndexRange[sl]);
BandedOperator/:SparseArray[bnd_BandedOperator,{n_,m_}]:=SparseArray[MapIndexed[SparseRule[#1,#2[[1]]]&,IncreaseLength[bnd,n][[1]]]//Flatten][[;;n,;;m]]


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




ApplyToRows[G_,Bn_BandedOperator,Bnn_BandedOperator,{row1_,row2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,Bn2,i,rhs1,rhs2},
Bn1=Bn;
Bn2=Bnn;

Do[
vals=G.{Bn1[[row1,i]],Bn2[[row2,i]]};


Bn1=ReplaceEntry[Bn1,{row1,i},vals[[1]],IncreaseSize->True];
Bn2=ReplaceEntry[Bn2,{row2,i},vals[[2]],IncreaseSize->True];
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

NullOperatorQ[G_,B1_,B2_,{srow1_,srow2_}]:=((((RowZeroQ[B1,srow1]||(NZeroQ[B1]))&&(RowZeroQ[B2,srow2]||(NZeroQ[B2])))||((NZeroQ[G[[1,2]]])&&(G[[2,2]]~NEqual~1)&&RowZeroQ[B1,srow1])))===True;

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

{B1,B2}=ApplyToRows[G,B1,B2,{srow1,srow2}];
Bn1=ReplaceEntry[Bn1,{row1,i},B1,IncreaseSize->True];
Bn2=ReplaceEntry[Bn2,{row2,i},B2,IncreaseSize->True];

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


ApplyToRows[G_,BL:{__BandedOperator},{row1_,row2_},{srow1_,srow2_},OptionsPattern[RightHandSide->Null]]:=Module[{vals,Bn1,B1,B2,rhs},
Bn1=BL;

If[row1===row2,
Bn1[[row1]]=ApplyToRows[G,Bn1[[row1]],{srow1,srow2}];
,
{B1,B2}=Bn1[[{row1,row2}]];
{B1,B2}=ApplyToRows[G,B1,B2,{srow1,srow2}];
Bn1[[row1]]=B1;
Bn1[[row2]]=B2;
];

If[OptionValue[RightHandSide]===Null,
Bn1
,
{Bn1,ApplyToRows[G,{row1,row2},{srow1,srow2},OptionValue[RightHandSide]]}
]
];


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
If[NZeroQ[a]&&NZeroQ[bB],
IdentityMatrix[2],
({
 {a, bB},
 {-bB, a}
})/Norm[{a,bB}]
]
];
Givens[Bn:{__BandedOperator},{row1_,row2_},{srow1_,srow2_},col_]:=Givens[Bn[[row1]],Bn[[row2]],{srow1,srow2},col];

Givens[Bn:{__BandedOperator},{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},col_,scol_]:=
Givens[Bn[[row1]][[srow1,col]],Bn[[row2]][[srow2,col]],{ssrow1,ssrow2},scol];



GivensReduce[BDx_,{row1_,row2_},col_,opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx,BDx,{row1,row2},col];
ApplyToRows[G,BDx,{row1,row2},opts]
];
GivensReduce[BDx_BandedOperator,{row1_,row2_},{srow1_,srow2_},col_,scol_,opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx[[row1,col]],BDx[[row2,col]],{srow1,srow2},scol];
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},opts]
];

GivensReduce[BDx:{__BandedOperator},{row1_,row2_},{srow1_,srow2_},scol_,opts:OptionsPattern[RightHandSide->Null]]:=Module[{G},
G=Givens[BDx,{row1,row2},{srow1,srow2},scol];
ApplyToRows[G,BDx,{row1,row2},{srow1,srow2},opts]
];

GivensReduce[BDx:{__BandedOperator},{row1_,row2_},{srow1_,srow2_},{ssrow1_,ssrow2_},col_,scol_,opts:OptionsPattern[RightHandSide->Null]]:=ApplyToRows[Givens[BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},col,scol],BDx,{row1,row2},{srow1,srow2},{ssrow1,ssrow2},opts];


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
DerivativeOperator[1]:=DerivativeOperator[1,Filler->AlternatingFiller];
DerivativeOperator[2]:=DerivativeOperator[2,Filler->AlternatingFiller];

ConversionOperator[1,Filler->fls_]:=BandedOperator[{ShiftList[{1,0,-1/2},0]},{0 fls[1]},ShiftList[{1/2,0,-1/2},1-#]&,Filler->fls];
ConversionOperator[1]:=ConversionOperator[1,Filler->AlternatingFiller];
ConversionOperator[2,Filler->fls_]:=BandedOperator[{ShiftList[{1,0,-2/3,0,1/6},0]},{0 fls[1]},ShiftList[{1/(2 #),0,-(1/(2 (#+2)))-1/(2 #),0,1/(2(#+2))},1-#]&,Filler->fls];
ConversionOperator[2]:=ConversionOperator[2,Filler->AlternatingFiller];

DirichletOperator[-1]:=BandedOperator[{ShiftList[{1},0]},{{1,0}},Null,Filler->AlternatingFiller];
DirichletOperator[1]:=BandedOperator[{ShiftList[{1},0]},{{0,1}},Null,Filler->AlternatingFiller];
NeumannOperator[-1]:=BandedOperator[{ShiftList[{0},0]},{{1,0}},Null,Filler->NeumannFiller];
NeumannOperator[1]:=BandedOperator[{ShiftList[{0},0]},{{0,1}},Null,Filler->NeumannFiller];


DirichletOperator[-1,Filler->fls_]:=BandedOperator[{ShiftList[{1},0]},{{1,0}~Join~(0 fls[1])},Null,Filler->(({(-1)^(#-1),1}~Join~fls[#])&)];
DirichletOperator[1,Filler->fls_]:=BandedOperator[{ShiftList[{1},0]},{{0,1 }~Join~(0 fls[1])},Null,Filler->(({(-1)^(#-1),1}~Join~fls[#])&)];

BasisOperator[1]:=BandedOperator[{ShiftList[{1},0]},{{0,0}},Null,Filler->AlternatingFiller];
BasisOperator[2]:=BandedOperator[{ShiftList[{0,1},0]},{{0,0}},Null,Filler->AlternatingFiller];


IdentityOperator[Filler->fls_]:=BandedOperator[{ShiftList[{1,0,0},0]},{0 fls[1]},ShiftList[{1,0,0},1-#]&,Filler->fls];
IdentityOperator[]:=IdentityOperator[Filler->AlternatingFiller];

ZeroOperator[Filler->fls_]:=BandedOperator[{ShiftList[{0,0,0},0]},{0 fls[1]},ShiftList[{0,0,0},1-#]&,Filler->fls];
ZeroOperator[]:=ZeroOperator[Filler->AlternatingFiller];

ZeroOperator[1,\[Infinity],Filler->fls_]:=BandedOperator[{ShiftList[{0},0]},{0 fls[1]},Null,Filler->fls];
ZeroOperator[1,\[Infinity]]:=ZeroOperator[1,\[Infinity],Filler->AlternatingFiller];






DirichletOperator[-1,All]:=BandedOperator[{ShiftList[{IdentityOperator[]},0]},{{IdentityOperator[],ZeroOperator[]}},Null,Filler->AlternatingFiller];
DirichletOperator[1,All]:=BandedOperator[{ShiftList[{IdentityOperator[]},0]},{{ZeroOperator[],IdentityOperator[]}},Null,Filler->AlternatingFiller];
DirichletOperator[All,-1]:=BandedOperator[{ShiftList[{DirichletOperator[-1]},0]},{{ZeroOperator[1,\[Infinity]],ZeroOperator[1,\[Infinity]]}},ShiftList[{DirichletOperator[-1]},1-#]&,Filler->AlternatingFiller];
DirichletOperator[All,1]:=BandedOperator[{ShiftList[{DirichletOperator[1]},0]},{{ZeroOperator[1,\[Infinity]],ZeroOperator[1,\[Infinity]]}},ShiftList[{DirichletOperator[1]},1-#]&,Filler->AlternatingFiller];


DerivativeOperator[0,1]:=BandedOperator[{ShiftList[{ConversionOperator[1]},-1]},{ZeroOperator[]},ShiftList[(#){ConversionOperator[1]},-#]&];
DerivativeOperator[1,0]:=BandedOperator[{ShiftList[{DerivativeOperator[1],ZeroOperator[],-DerivativeOperator[1]/2},0]},{ZeroOperator[]},ShiftList[{DerivativeOperator[1]/2,ZeroOperator[],-DerivativeOperator[1]/2},1-#]&];


DerivativeOperator[0,2]:=BandedOperator[{ShiftList[{4 ConversionOperator[2]},-2]},{ZeroOperator[]},ShiftList[{2(#+1)ConversionOperator[2]},-1-#]&];
DerivativeOperator[2,0]:=BandedOperator[{ShiftList[{DerivativeOperator[2],ZeroOperator[],-2/3 DerivativeOperator[2],ZeroOperator[],DerivativeOperator[2]/6},0]},{ZeroOperator[]},ShiftList[{DerivativeOperator[2]/(2 #),ZeroOperator[],(-(1/(2 (#+2)))-1/(2 #))DerivativeOperator[2],ZeroOperator[],DerivativeOperator[2]/(2(#+2))},1-#]&];


LaplaceOperator:=BandedOperator[{ShiftList[{DerivativeOperator[2],ZeroOperator[],BandedOperator[{ShiftList[(-2/3) {0,0,4,0,0}+4 {1,0,-2/3,0,1/6},0]},{{0,0}},ShiftList[(-2/3){0,0,2 (#+1),0,0}+4{1/(2 #),0,-(1/(2 (#+2)))-1/(2 #),0,1/(2(#+2))},1-#]&,Filler->AlternatingFiller],ZeroOperator[],DerivativeOperator[2]/6},0]},{{ZeroOperator[],ZeroOperator[]}},ShiftList[{DerivativeOperator[2]/(2 #),ZeroOperator[],BandedOperator[{ShiftList[(-(1/(2 (#+2)))-1/(2 #)) {0,0,4,0,0}+2(#+1){1,0,-2/3,0,1/6},0]},{{0,0}},Function[rw,ShiftList[(-(1/(2 (#+2)))-1/(2 #)){0,0,2 (rw+1),0,0}+2(#+1){1/(2rw),0,-(1/(2 (rw+2)))-1/(2 rw),0,1/(2(rw+2))},1-rw]],Filler->AlternatingFiller],ZeroOperator[],DerivativeOperator[2]/(2(#+2))},1-#]&,Filler->AlternatingFiller];


SparseHankelMatrix[f_List,n_]:=PadRight[HankelMatrix[f]//SparseArray,{n,n}];SparseHankelMatrix[f_IFun,n_]:=SparseHankelMatrix[f//DCT,n];
HankelZeroFirst[M_,m_]:=Module[{A},
A=M;
A[[1,;;Min[m,Length[M]]]]=0.;
A];

SparseToeplitzMatrix[g_LFun?ScalarFunQ,n_]:=SparseToeplitzMatrix[g//FFT,n];
SparseToeplitzMatrix[f_IFun?ScalarFunQ,n_]:=SparseToeplitzMatrix[ShiftList[Reverse[f/2.//DCT//Rest],f/2.//DCT//DoubleFirst],n];
SparseToeplitzMatrix[sl_ShiftList,n_]:=SparseArray[ToList[MapIndexed[{i_,j_}/;i==j-{##}[[2,1]]->{##}[[1]]&,ChopDrop[sl//Reverse,$MachineTolerance]]],{n,n}];

SparseToeplitzMatrix[l1_List,l2_List,n_]:=Module[{k},
SparseArray[Table[Band[{k,1}]->l1[[k]],{k,Min[n,Length[l1]]}]~Join~Table[Band[{1,k}]->l1[[k]],{k,2,Min[n,Length[l2]]}],n]
];

MultiplicationMatrix[v_List,n_]:=(HankelZeroFirst[SparseHankelMatrix[v/2.,n],v//Length]+SparseToeplitzMatrix[(v//DoubleFirst)/2,(v//DoubleFirst)/2,n]);

MultiplicationMatrix[f_IFun,n_]:=MultiplicationMatrix[f//DCT,n];


CDerivativeOperator[2,Filler->fls_]:=BandedOperator[{ShiftList[{0},0]},{(0 fls[1])~Join~{4,0,0,0}},ShiftList[{0},1-#]&,Filler->((fls[#]~Join~CDerivativeFiller[2][#])&),FillGenerator->Function[row,(0 fls[1])~Join~\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{"{", 
RowBox[{"0", ",", "0", ",", "4", ",", 
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"row", "-", "2"}], ")"}]}], 
RowBox[{"(", "row", ")"}]}]}], "}"}], 
RowBox[{"EvenQ", "[", "row", "]"}]},
{
RowBox[{"{", 
RowBox[{"8", ",", 
RowBox[{
RowBox[{"-", "2"}], " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{"row", "-", "1"}], ")"}], "2"]}], ",", "0", ",", "0"}], "}"}], 
RowBox[{"OddQ", "[", "row", "]"}]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False]\)]];
CDerivativeOperator[2]:=BandedOperator[{ShiftList[{0},0]},{{4,0,0,0}},ShiftList[{0},1-#]&,Filler->CDerivativeFiller[2],FillGenerator->Function[row,\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{"{", 
RowBox[{"0", ",", "0", ",", "4", ",", 
RowBox[{
RowBox[{"-", 
RowBox[{"(", 
RowBox[{"row", "-", "2"}], ")"}]}], 
RowBox[{"(", "row", ")"}]}]}], "}"}], 
RowBox[{"EvenQ", "[", "row", "]"}]},
{
RowBox[{"{", 
RowBox[{"8", ",", 
RowBox[{
RowBox[{"-", "2"}], " ", 
SuperscriptBox[
RowBox[{"(", 
RowBox[{"row", "-", "1"}], ")"}], "2"]}], ",", "0", ",", "0"}], "}"}], 
RowBox[{"OddQ", "[", "row", "]"}]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False]\)]];


DropRow[bnd:BandedOperator[A_List,fill_List,rowgen_,Filler->fls_],row_]:=If[row>Length[A]||Length[A]==1,DropRow[bnd//IncreaseLength,row],BandedOperator[Drop[A,{row}],Drop[fill,{row}],rowgen[#+1]&,Filler->fls]];

DropColumn[bnd:BandedOperator[A_List,fill_List,rowgen_,Filler->fls_],col_]:=BandedOperator[Drop[#,col]&/@A,fill,rowgen,Filler->fls];


End[];
EndPackage[];
