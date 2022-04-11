%script to setup fitting of params

numCat=7;
Tw=1; %1 year b/c data is year-to-year freq
dt=1/365;

% -- for cov matrix; putting entries into 1 col vec -- 
%indicies for Cov (non-var)
linInd=[];
for rwI=1:(numCat-1)
    for clI=rwI+1:numCat
        linInd=[linInd; sub2ind([numCat numCat],rwI,clI)];
    end
end
[rwInd,clInd]=ind2sub([numCat numCat],linInd); clear rwI clI;

%freq distr params
a_prm=zeros(numCat,1);
gam_prm=zeros(numCat,1);
tau_prm=zeros(numCat,1);
cInp_prm=zeros(numCat); %only upper trian filled in

load datORX_meanVar.mat %need avgFreq,varFreq,covM_freq

covFreak=triu(covM_freq,1); %21 cov of freqs used to fit
%avgFreq
%varFreq

%put all ORX data into 1 structure
orxDat=struct('avgF',avgFreq,'varF',varFreq,'covF',covFreak);
% put all model params into 1 struct
%fParms=struct('a',a_prm,'gam',gam_prm,'tau',tau_prm);

%setup initial guess
a_0=4*ones(numCat,1); tau_0=0.5*ones(numCat,1); gam_0=avgFreq./tau_0; 
cin_0=zeros(numCat*(numCat-1)/2,1); %21 x 1
denom_C0=min(gam_0(rwInd),gam_0(clInd)).*a_0(rwInd).*a_0(clInd).*tau_0(rwInd).*tau_0(clInd)./...
    (tau_0(rwInd)+tau_0(clInd)).*(tau_0(rwInd).*(1+tau_0(rwInd).*(exp(-1./tau_0(rwInd))-1))+...
    tau_0(clInd).*(1+tau_0(clInd).*(exp(-1./tau_0(clInd))-1)));
cin_0=covFreak(linInd)./denom_C0;
%var_0=(a_0.*tau_0).^2.*gam_0.*(1+tau_0.*(exp(-1./tau_0)-1)); %not used

fPrm0=[a_0; tau_0; gam_0; cin_0]; %is 42x1 vector, start the optimz here.

fnToMin=@(x)objFunc(x,orxDat);

optins=optimoptions('fmincon','Display','final-detailed','MaxIteration',1e6,'MaxFunctionEvaluations',1e6);
lowB=eps*ones(size(fPrm0)); %lower bound, a,tau,gam all > 0
lowB( end-numCat*(numCat-1)/2+1 : end )=-1+eps; %for cin
upB=Inf(size(fPrm0)); 
upB( end-numCat*(numCat-1)/2+1 : end )= 1-eps; %for cin
nonlcon=[]; %no nonlinear Constr
A=[]; b=[]; Aeq=[]; beq=[];

% perform optimization
[fPrm,fval,exitflag,output]=fmincon(fnToMin,fPrm0,A,b,Aeq,beq,lowB,upB,nonlcon,optins);

%% get pieces of model fits
CovMod=zeros( numCat*(numCat-1)/2 , 1 );
a_prm=fPrm(1:numCat);
tau_prm=fPrm(numCat+1:2*numCat);
gam_prm=fPrm(2*numCat+1:3*numCat);
cij_prm=fPrm(3*numCat+1:end); %21 x 1 vect
%theor cov
for j=1:length(linInd) 
    %set all 21 theor cov here, use all parms (a,tau,gam,cij_prm)
    tau1=tau_prm(rwInd(j)); tau2=tau_prm(clInd(j));  %simpler to write tau1, tau2
    a1=a_prm(rwInd(j)); a2=a_prm(clInd(j)); 
    gam1=gam_prm(rwInd(j)); gam2=gam_prm(clInd(j));
    
    CovMod(j) = cij_prm(j)*min(gam1,gam2)*a1*a2*((tau1*tau2)/(tau1+tau2))*(tau1*(Tw+tau1*(exp(-Tw/tau1)-1))+tau2*(Tw+tau2*(exp(-Tw/tau2)-1)));
end
VarMod=a_prm.^2.*tau_prm.^2.*gam_prm.*(Tw+tau_prm.*(exp(-Tw./tau_prm)-1));

save fittedModelParams CovMod VarMod a_prm tau_prm gam_prm cij_prm

figure(1);
hold on
plot(a_prm.*tau_prm.*gam_prm,'r*-');
plot(orxDat.avgF,'bo-');
set(gca,'FontSize',18)
set(gca,'XLim',[1-eps numCat+eps])
ylabel('Mean Freq')

figure(2);
hold on
plot(VarMod,'r*-');
plot(orxDat.varF,'bo-');
set(gca,'FontSize',18)
set(gca,'XLim',[1-eps numCat+eps])
ylabel('Var Freq')

figure(3);
hold on
plot(CovMod,'r*-')
plot(orxDat.covF(linInd),'bo-');
set(gca,'FontSize',18)
set(gca,'XLim',[1-eps 21+eps])



