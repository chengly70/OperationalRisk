%script to setup fitting of params, starting with 
% many more initial conditions


nmParms=1000; %fitting this many times

numCat=7;
Tw=1; %1 year b/c data is year-to-year freq
dt=1/365;

% outputs to save
CovMod=zeros(numCat*(numCat-1)/2 , nmParms );
cij_prm=zeros(numCat*(numCat-1)/2 , nmParms );
VarMod=zeros(numCat, nmParms);
a_prm=zeros(numCat, nmParms);
tau_prm=zeros(numCat, nmParms);
gam_prm=zeros(numCat, nmParms);

totParms=size(cij_prm,1)+3*numCat; %42=21+3*7

% -- for cov matrix; putting entries into 1 col vec -- 
%indicies for Cov (non-var)
linInd=[];
for rwI=1:(numCat-1)
    for clI=rwI+1:numCat
        linInd=[linInd; sub2ind([numCat numCat],rwI,clI)];
    end
end
[rwInd,clInd]=ind2sub([numCat numCat],linInd); clear rwI clI;

load datORX_meanVar.mat %need avgFreq,varFreq,covM_freq

covFreak=triu(covM_freq,1); %21 cov of freqs used to fit

%put all ORX data into 1 structure
orxDat=struct('avgF',avgFreq,'varF',varFreq,'covF',covFreak);
% put all model params into 1 struct
%fParms=struct('a',a_prm,'gam',gam_prm,'tau',tau_prm);

fnToMin=@(x)objFunc(x,orxDat);

optins=optimoptions('fmincon','Display','final-detailed','MaxIteration',1e6,'MaxFunctionEvaluations',1e6);
lowB=eps*ones(totParms,1); %lower bound, a,tau,gam all > 0
lowB( end-numCat*(numCat-1)/2+1 : end )=-1+eps; %for cin
upB=Inf(totParms,1); 
upB( end-numCat*(numCat-1)/2+1 : end )= 1-eps; %for cin
nonlcon=[]; %no nonlinear Constr
A=[]; b=[]; Aeq=[]; beq=[];

%seed rand # gen so reproducible
rng(1880);

for rnic=1:nmParms
    %setup initial guess
    a_0=15*rand(numCat,1); tau_0=2*rand(numCat,1); gam_0=40*rand(numCat,1);
    cin_0=2*rand( numCat*(numCat-1)/2 , 1) - 1; 
    %var_0=(a_0.*tau_0).^2.*gam_0.*(1+tau_0.*(exp(-1./tau_0)-1)); %not used
    fPrm0=[a_0; tau_0; gam_0; cin_0]; %is 42x1 vector, start the optimz here.
    % perform optimization
    [fPrm,fval,exitflag,output]=fmincon(fnToMin,fPrm0,A,b,Aeq,beq,lowB,upB,nonlcon,optins);
    
    % get pieces of model fits
    a_prm(:,rnic)=fPrm(1:numCat);
    tau_prm(:,rnic)=fPrm(numCat+1:2*numCat);
    gam_prm(:,rnic)=fPrm(2*numCat+1:3*numCat);
    cij_prm(:,rnic)=fPrm(3*numCat+1:end); %21 x 1 vect
    %theor cov
    for j=1:length(linInd)
        %set all 21 theor cov here, use all parms (a,tau,gam,cij_prm)
        tau1=tau_prm(rwInd(j),rnic); tau2=tau_prm(clInd(j),rnic);  %simpler to write tau1, tau2
        a1=a_prm(rwInd(j),rnic); a2=a_prm(clInd(j),rnic);
        gam1=gam_prm(rwInd(j),rnic); gam2=gam_prm(clInd(j),rnic);
        
        CovMod(j,rnic) = cij_prm(j,rnic)*min(gam1,gam2)*a1*a2*((tau1*tau2)/(tau1+tau2))*(tau1*(Tw+tau1*(exp(-Tw/tau1)-1))+tau2*(Tw+tau2*(exp(-Tw/tau2)-1)));
    end
    VarMod(:,rnic)=a_prm(:,rnic).^2.*tau_prm(:,rnic).^2.*gam_prm(:,rnic).*(Tw+tau_prm(:,rnic).*(exp(-Tw./tau_prm(:,rnic))-1));
    
end

save('RandPrm0_fits.mat','CovMod','VarMod','a_prm','tau_prm','gam_prm','cij_prm')




