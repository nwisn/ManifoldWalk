(* ::Package:: *)

LeviCivitaCoefficients[metric__,coords__]:=Module[{mu,nu,lambda},Table[1/2 Sum[Inverse[metric][[lambda,dummy]](\!\(
\*SubscriptBox[\(\[PartialD]\), \(coords[\([nu]\)]\)]\(metric[\([dummy, mu]\)]\)\)+\!\(
\*SubscriptBox[\(\[PartialD]\), \(coords[\([mu]\)]\)]\(metric[\([dummy, nu]\)]\)\)-\!\(
\*SubscriptBox[\(\[PartialD]\), \(coords[\([dummy]\)]\)]\(metric[\([mu, nu]\)]\)\)),{dummy,1,Dimensions[coords][[1]]}],{mu,1,Dimensions[coords][[1]]},{nu,1,Dimensions[coords][[1]]},{lambda,1,Dimensions[coords][[1]]}]]//FullSimplify


EulerLagrangeEquations[coefficients__,coords__,param_]:=Module[{mu,nu,lambda},Table[coords[[lambda]]''[param]+Sum[coefficients[[mu]][[nu]][[lambda]] coords[[mu]]'[param]coords[[nu]]'[param],{mu,1,Dimensions[coords][[1]]},{nu,1,Dimensions[coords][[1]]}],{lambda,1,Dimensions[coords][[1]]}]]//FullSimplify


FisherMetric[pdf_,coords__,assumptions__,limits__]:=Module[{mu,nu},Table[Assuming[assumptions,-\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(limits[\([1]\)]\), \(limits[\([2]\)]\)]\(pdf\ 
\*SubscriptBox[\(\[PartialD]\), \(coords[\([mu]\)]\)]\(
\*SubscriptBox[\(\[PartialD]\), \(coords[\([nu]\)]\)]Log[pdf]\) \[DifferentialD]x\)\)],{mu,1,Dimensions[coords][[1]]},{nu,1,Dimensions[coords][[1]]}]]//FullSimplify


AlphaCoefficients[pdf_,coords__,assumptions__,limits__]:=Module[{mu,nu,lambda},Table[Assuming[assumptions,\!\(
\*SubsuperscriptBox[\(\[Integral]\), \(limits[\([1]\)]\), \(limits[\([2]\)]\)]\(pdf \((\((
\*SubscriptBox[\(\[PartialD]\), \(coords[\([mu]\)]\)]\(
\*SubscriptBox[\(\[PartialD]\), \(coords[\([nu]\)]\)]Log[pdf]\) + FractionBox[\((1 - \[Alpha])\), 2] 
\*SubscriptBox[\(\[PartialD]\), \(coords[\([mu]\)]\)]Log[pdf] 
\*SubscriptBox[\(\[PartialD]\), \(coords[\([nu]\)]\)]Log[pdf])\) \((
\*SubscriptBox[\(\[PartialD]\), \(coords[\([lambda]\)]\)]Log[pdf])\))\) \[DifferentialD]x\)\)],{mu,1,Dimensions[coords][[1]]},{nu,1,Dimensions[coords][[1]]},{lambda,1,Dimensions[coords][[1]]}]]//FullSimplify
