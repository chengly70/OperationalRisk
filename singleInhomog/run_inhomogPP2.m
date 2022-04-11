% script to get data for Fig 2 (homog PP)
% calls homogPoissonP.m -- freq=PP & Severity is 1 of 5 dists
% here using GPD (see run_homogPP.m and drv_homogPoissP.m)

tEnd=100; %time in years
numrel=50000;
dt=0.005; %less than a month
tmvc=(0:dt:tEnd)';
ap=1;
tau_d=1.2;
lgamm=75/2;  
%avg rate (# per year), keep ap*tau_d*nu_i=75
Tw=1; 

%severDistr=0;  %0,gamm=gamma; 1,logn=lognormal; 2,gpd=GPD; 3,wbl=Weibull; 4,burr=Burr
%severDistr='gamm'; %can also use char
%params=[ 2; 1]; %gamm[alph;bet], logn[mu;sig], gpd[kpar,sig], wbl[a,b], burr[alph,cpar,kpar]
%params=[1; 2; 2];

flName='dInhomPP_gpd';
severDistr=2; 
kpar=[0.1; 0.15; 0.2];
spar=.85*(60:5:110)';

%outputs 
mnR=zeros(length(spar),length(kpar));
vrR=zeros(length(spar),length(kpar));
mnTw=zeros(length(spar),length(kpar));
vrTw=zeros(length(spar),length(kpar));
mnR_an=zeros(length(spar),length(kpar));
vrR_an=zeros(length(spar),length(kpar));
mnTw_an=zeros(length(spar),length(kpar));
vrTw_an=zeros(length(spar),length(kpar));

for kin=1:length(kpar)
tic
    for jin=1:length(spar)
    
        params=[ kpar(kin) ; spar(jin)];
        
    [mean_R,var_R,mean_Tw,var_Tw,an_mean,an_var,an_meanTw,an_varTw]=InhomogPoissonP_mr(dt,tEnd,ap,lgamm,tau_d,Tw,severDistr,params,numrel);
        %save data
        mnR(jin,kin)=mean_R;
        vrR(jin,kin)=var_R;
        mnTw(jin,kin)=mean_Tw;
        vrTw(jin,kin)=var_Tw;
        mnR_an(jin,kin)=an_mean;
        vrR_an(jin,kin)=an_var;
        mnTw_an(jin,kin)=an_meanTw;
        vrTw_an(jin,kin)=an_varTw;
    end
toc
    %save mat file
    save([pwd,'/',flName],'mn*','vr*','spar','kpar','lgamm','ap','tau_d','Tw');
end



