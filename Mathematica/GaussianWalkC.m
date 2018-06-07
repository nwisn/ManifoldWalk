(* Mathematica package *)

gaussianwalk[dimension_, datasamplesize_, outputsize_, oderunsperpoint_, stepsperoderun_] :=
    Block[ {outputmu,outputsigma, h,s,stretch,zeroarray,basemu,baseproduct,
    		outputnumber,currentmu,bigproduct,
    		odenumber,mudirection,sigmadirection,
    		deltamu,L},
    	
    	outputmu = ConstantArray[0.,{outputsize,dimension}];
    	outputsigma = ConstantArray[0.,{outputsize,dimension,dimension}];
    	
        h = N[1/stepsperoderun];
        s = N[Sqrt[1/(datasamplesize*oderunsperpoint)]];
        stretch = N[Sqrt[2]*s];
        zeroarray = ConstantArray[0., {dimension, dimension}];
        basemu = ConstantArray[0., dimension];
        baseproduct = N[IdentityMatrix[dimension]];
        
        For[outputnumber = 1, outputnumber <= outputsize, outputnumber++,
           	currentmu = basemu;
           	bigproduct = baseproduct;
           	For[odenumber = 1, odenumber <= oderunsperpoint, odenumber++,
            	mudirection = getmudirection[s,dimension];
            	sigmadirection = getsigmadirection[s, dimension];
            	{deltamu, L} = simulatethendecompose[{basemu,baseproduct},
									            	 {mudirection, sigmadirection},
									            	 eulerlagrangefield,
									            	 h,
									            	 stepsperoderun];
            	currentmu = L.(currentmu + LTsolve[L, deltamu]);
            	bigproduct = L.bigproduct;
            ];
            outputmu[[outputnumber]] = currentmu;
            outputsigma[[outputnumber]] = bigproduct.Transpose[bigproduct];
        ];
   		Return[{outputmu,outputsigma}];
	];
	
getsigmadirection[std_,dim_] := 
 	Block[{mat, upper, diag}, 
  		mat = RandomVariate[NormalDistribution[0, std], {dim, dim}];
  		upper = UpperTriangularize[mat, 1];
  		diag = Sqrt[2.]*DiagonalMatrix[Diagonal[mat]];
  		Return[diag + upper + Transpose[upper]];
  	];

getmudirection[std_, dim_] := 
	Block[{},
   		Return[RandomVariate[NormalDistribution[0, std], dim]];
   	];
	
simulatethendecompose[initialvalue_, initialtangent_, accelfield_, stepsize_, numsteps_] := 	
	Block[ {currentvalue, currenttangent, doublederiv, t},
  		currentvalue = initialvalue;
  		currenttangent = initialtangent;
  		doublederiv = accelfield[currentvalue, currenttangent];
     	currentvalue = currentvalue + stepsize*currenttangent;
  		For[t = 1, t <= numsteps - 1, t++,
   			currenttangent = currenttangent + stepsize*doublederiv;
   			doublederiv = accelfield[currentvalue, currenttangent];
   			currentvalue = currentvalue + stepsize*currenttangent;
   		];
  		Return[{currentvalue[[1]], Transpose[CholeskyDecomposition[currentvalue[[2]]]]}];
  	];
	
	
eulerlagrangefield[point_, tangent_] := 
	Block[ {mu, sigma, dmu, dsigma, logarithmicderivative},
  		{mu, sigma} = point;
  		{dmu, dsigma} = tangent;
  		logarithmicderivative = dsigma.Inverse[sigma];
  		Return[{logarithmicderivative.dmu, 
  				logarithmicderivative.dsigma - dmu.dmu
  				}];
  	];
  	
LTsolve[a_, b_] := 
	Block[ {x,d},
		d = Dimensions[a][[1]];
  		x = ConstantArray[0, d];
  		For[i = 1, i <= d, i++,
   			x[[i]] = (b[[i]] - ((a[[i, 1 ;; i]].x[[1 ;; i]])))/a[[i, i]]
   		];
  	Return[x];
  	];
  	
LTsolve2[a_, b_] := LinearSolve[LowerTriangularize[a], b]  	