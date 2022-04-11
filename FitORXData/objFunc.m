function err=objFunc(fPrm,orxDat)
% INPUT: fPrm is [a,tau,gam,c_inp], 42x1, orxDat is a structure
%

numCat=7; %!!!assuming this wont change!!!
%dt=1/365; %not used?
Tw=1; %1 year b/c data is year-to-year freq


%passed in, input from fPrm
a_prm=fPrm(1:numCat);
tau_prm=fPrm(numCat+1:2*numCat);
gam_prm=fPrm(2*numCat+1:3*numCat);
cij_prm=fPrm(3*numCat+1:end); %21 x 1 vect

% get ORX data
avgFreq=orxDat.avgF;
varFreq=orxDat.varF;
covFrks=orxDat.covF;

% -- for cov matrix; putting entries into 1 col vec -- 
%indicies for Cov (non-var)
linInd=[];
for rwI=1:(numCat-1)
    for clI=rwI+1:numCat
        linInd=[linInd; sub2ind([numCat numCat],rwI,clI)];
    end
end
[rwInd,clInd]=ind2sub([numCat numCat],linInd); clear rwI clI;

%35 x 1
rhsData=[avgFreq; covFrks(linInd); varFreq]; %values from ORXdata

lhs_Model=zeros(size(rhsData));
% theoret mean, 1st 7 entries are mean
lhs_Model(1:numCat) = a_prm.*tau_prm.*gam_prm; %set that here, use a_prm, tau_prm, gam_prm
%theor cov
for j=1:length(linInd) 
    %set all 21 theor cov here, use all parms (a,tau,gam,cij_prm)
    tau1=tau_prm(rwInd(j)); tau2=tau_prm(clInd(j));  %simpler to write tau1, tau2
    a1=a_prm(rwInd(j)); a2=a_prm(clInd(j)); 
    gam1=gam_prm(rwInd(j)); gam2=gam_prm(clInd(j));
    
    lhs_Model(numCat+j) = cij_prm(j)*min(gam1,gam2)*a1*a2*((tau1*tau2)/(tau1+tau2))*(tau1*(Tw+tau1*(exp(-Tw/tau1)-1))+tau2*(Tw+tau2*(exp(-Tw/tau2)-1)));
end
% theoret var; last 7 entries are var, set here using a_prm,tau_prm,gam_prm
lhs_Model(end-numCat+1:end) = a_prm.^2.*gam_prm.*tau_prm.^2.*(Tw+tau_prm.*(exp(-Tw./tau_prm)-1)); 

err = sum( abs(lhs_Model - rhsData) ); %L-1 norm


