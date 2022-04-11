%driver script that calls InhomTwoPP_mr.m -- freq=PP & Severity is 1 of 5 dists
% here varying c_inp for fixed Sevs

flName='dVary_cinA';

tEnd=10000; %time in years
numrel=2000;
dt=0.005; %less than a day
tmvc=(0:dt:tEnd)';
Lt=round(tEnd/dt)+1;

%Tw=2; %larger time window, %for dVary_cin.mat
Tw=1; %usual time window   %for dVary_cinA.mat

c_inp=[-.8; -.6; -.4; -.2; .2; .4; .6; .8]; %correl of 2 frequencies; has to be in [-1,1]
lenCin=length(c_inp);

%severDistr=0;  %0,gamm=gamma; 1,logn=lognormal; 2,gpd=GPD; 3,wbl=Weibull; 4,burr=Burr
%severDistr='gamm'; %can also use char
%params=[ 2; 1]; %gamm[alph;bet], logn[mu;sig], gpd[kpar,sig], wbl[a,b], burr[alph,cpar,kpar]
% severDistr=1; params=[1; .5];
% severDistr=4;  params=[1; 2; 2]; %
paramsL1=struct('a',1.5,'nu_i',30,'tau_d',1.3,'SeverDist',2,'Sparams',[.15; 50]);
paramsL2=struct('a',2,'nu_i',40,'tau_d',.75,'SeverDist',3,'Sparams',[5; .4]);

%outputs 
mnR=zeros(lenCin,2);
vrR=zeros(lenCin,2);
mnTw=zeros(lenCin,2);
vrTw=zeros(lenCin,2);
mnR_an=zeros(lenCin,2);
vrR_an=zeros(lenCin,2);
mnTw_an=zeros(lenCin,2);
vrTw_an=zeros(lenCin,2);
covR=zeros(lenCin,1);
covR_an=zeros(lenCin,1);
covTw=zeros(lenCin,1);
covTw_an=zeros(lenCin,1);

for jin=1:lenCin
tic
    [mean_R,var_R,mean_Tw,var_Tw,an_mean,an_var,an_meanTw,an_varTw,...
        cov_R,an_cov,cov_Tw,an_covTw]=InhomTwoPP_mr(dt,tEnd,Tw,c_inp(jin),paramsL1,paramsL2,numrel);
toc
    %save data
        mnR(jin,:)=mean_R';
        vrR(jin,:)=var_R';
        mnTw(jin,:)=mean_Tw';
        vrTw(jin,:)=var_Tw';
        mnR_an(jin,:)=an_mean';
        vrR_an(jin,:)=an_var';
        mnTw_an(jin,:)=an_meanTw';
        vrTw_an(jin,:)=an_varTw';
        covR(jin,1)=cov_R;
        covR_an(jin,1)=an_cov;
        covTw(jin,1)=cov_Tw;
        covTw_an(jin,1)=an_covTw;

    %save mat file
    save([pwd,'/',flName],'mn*','vr*','covR*','covTw*','Tw','c_inp','paramsL1','paramsL2');
end
    

