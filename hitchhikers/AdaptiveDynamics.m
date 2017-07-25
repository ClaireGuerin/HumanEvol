(* :Title: Adaptive Dynamics Toolkit *)

(* :Context: "AdaptiveDynamics`" *)

(* :Authors:
        Ake Brannstrom, Nicolas Loeuille
	Bug reports to brnstrom@iiasa.ac.at
*)

(* :Summary:
	This package defines functions for drawing pairwise invasibility 
	plots (PIPs) and mutual invasibility plots (MIPs).
*)
(* :History:
	051221 Package first sees the light of day.
*)
(* :Mathematica Version: 4.0 *)
(* :Keywords:    *) 
(* :Suggestions: Add arrows and plus-signs using mouse. *) 
(* :Discussion:  *)

BeginPackage["AdaptiveDynamics`", {"Graphics`Arrow`", "Graphics`ImplicitPlot`", "NumericalMath`NLimit`"}]


TestInvExp::usage =
	"TestInvExp[S, r0, r1] runs a series of tests on the invasion exponent S."

NFindSP::usage =
	"NFindSP[S, r0, r1] finds critical points of S."

DrawPIP::usage =
	"DrawPIP[S, r0, r1] draws a pairwise invasibility plot."

DrawRC::usage =
	"DrawRC[S, r0, r1] draws the region of coexistence."

DrawMIP::usage =
	"DrawMIP[S, r0, r1] draws a mutual invasibility plot."

AddPlus::usage = 
        "AddPlus[gfx, x, y, font_size] adds the plus sign at the specified coordinates."

AddArrow::usage = 
	"AddArrow[gfx_, S_, x1_, x2_, h_] adds arrows showing the direction
	 of evolution from the point (x1, x2) and hopefully also the symmetric case (x2, x1)."

Begin["`Private`"]

AddPlus[gfx_, x_, y_, fs_] :=
  Show[{gfx, Graphics[{Text[StyleForm["+", FontSize ->fs], {x, y}]}]}]
  
(* Default argument for AddPlus *)
AddPlus[gfx_, x_, y_] := AddPlus[gfx, x, y, 40]


AddArrow[gfx_, S_, x1_, x2_, h_] := 
  Module[{m, h1, h2},
    h1 := h*Sign[ND[S[x1, x2, m], m, x1]];
    h2 := h*Sign[ND[S[x1, x2, m], m, x2]];
     Show[gfx, 
      Graphics[{Arrow[{x1, x2}, {x1 + h1, x2}], 
          Arrow[{x1, x2}, {x1, x2 + h2}]}]]
    ]

DrawPIP[S_, r0_, r1_, numPoints_] := 
  Module[{r, m}, 
    ContourPlot[S[r, m], {r, r0, r1}, {m, r0, r1},
      Contours->{0},
      ColorFunction->(GrayLevel[1 - 0.1 (#)]&),
      PlotPoints->numPoints,
      FrameLabel -> {"Resident trait", "Mutant trait"}
    ]
  ]

(* Default argument for numPoints *)
DrawPIP[S_, r0_, r1_] := DrawPIP[S, r0, r1, 50]


(* Draw region of coexistence *)
DrawRC[S_, r0_, r1_, numPoints_, extraPlotOptions_] := 
    Module[{x, y, Q}, 
      Q[x_, y_] = If[Positive[S[x, y]] && Positive[S[y, x]], 1.0, -1.0];
      ContourPlot[Q[x, y], {x, r0, r1}, {y, r0, r1}, Contours->{0},
		  PlotPoints->numPoints, ColorFunction->(GrayLevel[1-0.1(#)] &),
                  extraPlotOptions]
    ]

(* TODO: DrawRC: Allow specification of color *)

DrawRC[S_, r0_, r1_, numPoints_] := 
    DrawRC[S, r0, r1, numPoints, DisplayFunction->$DisplayFunction]

(* Default argument for numPoints *)
DrawRC[S_, r0_, r1_] := DrawRC[S, r0, r1, 200]


DrawMIP[S_, r0_, r1_, numPointsCoex_, numPointsIso_] := 
  Module[{m, S1Grad, S2Grad, regCoex, iso1, iso2, pickMaxima1, pickMaxima2, isoData1, 
          isoData2, allIsos1, allIsos2, maxIsos1, maxIsos2, allSel1, allSel2, epsilon, maxSel1, maxSel2}, 

    S1Grad[x1_, x2_] := 
      If[S[x1, x2] <= 0. || S[x2, x1] <= 0.,
        100., 
        ND[S[x1,x2, m], m, x1]
      ];

    S2Grad[x1_, x2_] := 
      If[S[x1, x2] <= 0. || S[x2, x1] <= 0.,
        100.,
        ND[S[x1, x2, m], m, x2]
      ];

    regCoex := DrawRC[S, r0, r1, numPointsCoex, DisplayFunction->Identity]; 

    Off[FindRoot::"lstol"];

    (* Should be easy to exploit symmetry. Use just one and then 
       flip x,y to y,x. *)
    iso1 := ImplicitPlot[S1Grad[x, y] == 0, {x, r0, r1}, {y, r0, r1}, PlotPoints->numPointsIso, DisplayFunction->Identity]; 
    iso2 := ImplicitPlot[S2Grad[x, y] == 0, {x, r0, r1}, {y, r0, r1}, PlotPoints->numPointsIso, DisplayFunction->Identity];

    isoData1 = Map[#[[3, 1]] &, First[Graphics[iso1]]];
    isoData2 = Map[#[[3, 1]] &, First[Graphics[iso2]]];

    epsilon := 10^(-4);
    PickZeros[pointData_] := Select[pointData, Abs[ND[S[#[[1]], #[[2]], m], m, #[[1]]]] < epsilon || Abs[ND[S[#[[1]], #[[2]], m], m, #[[2]]]] < epsilon &];
    PickMaxima1[pointData_] := Select[pointData, ND[S[#[[1]], #[[2]], m], {m, 2}, #[[1]]] < 0 &];
    PickMaxima2[pointData_] := Select[pointData, ND[S[#[[1]], #[[2]], m], {m, 2}, #[[2]]] < 0 &];
    
    allSel1 = DeleteCases[Map[PickZeros, isoData1], {}_];
    allSel2 = DeleteCases[Map[PickZeros, isoData2], {}_];

    maxSel1 = DeleteCases[Map[PickMaxima1, allSel1], {}_];
    maxSel2 = DeleteCases[Map[PickMaxima2, allSel2], {}_];

    allIsos1 = Map[ListPlot[#, PlotStyle -> PointSize[0.01], DisplayFunction->Identity] &, allSel1];
    allIsos2 = Map[ListPlot[#, PlotStyle -> PointSize[0.01], DisplayFunction->Identity] &, allSel2];

    maxIsos1 = Map[ListPlot[#, PlotStyle -> PointSize[0.02], DisplayFunction->Identity] &, maxSel1];
    maxIsos2 = Map[ListPlot[#, PlotStyle -> PointSize[0.02], DisplayFunction->Identity] &, maxSel2];

  (* Find the points that are thick and plot them as large points *)
(*     contourLines = Join[First[Graphics[iso1]], First[Graphics[iso2]]]; *)
    On[FindRoot::"lstol"];

    Show[regCoex, allIsos1, allIsos2, maxIsos1, maxIsos2, DisplayFunction->$DisplayFunction]
(*     Show[regCoex, iso1, iso2, maxIsos, DisplayFunction->$DisplayFunction] *)
]
  
    
(* Default argument for numPoints *)
DrawMIP[S_, r0_, r1_] := DrawMIP[S, r0, r1, 200, 200];
DrawMIP[S_, r0_, r1_, numPointsCoex_] := DrawMIP[S, r0, r1, numPointsCoex, 200];
End[]


EndPackage[]
