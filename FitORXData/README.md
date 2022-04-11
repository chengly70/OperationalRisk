Scripts to fit our model to ORX data

fitFreq_Parms.m -- main script to fit, loads ORX data datORX_meanVar.mat, saves results in fittedModelParams.mat; calls objFunc.m

objFunc.m -- defines the objective function, passing in ORX data (defines L1-norm)
