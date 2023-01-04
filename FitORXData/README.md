Scripts to fit our model to ORX data

MORE_fitFreq_Parms.m -- main script to fit WITH 1000 (or however many you want) random initial paramertizations, loads ORX data datORX_meanVar.mat, saves results in fittedModelParams.mat; calls objFunc.m

objFunc.m -- defines the objective function, passing in ORX data (defines L1-norm)

Don't use this, old routine. 
fitFreq_Parms.m -- main script to fit, loads ORX data datORX_meanVar.mat, saves results in fittedModelParams.mat; calls objFunc.m

