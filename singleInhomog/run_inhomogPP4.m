% script to get data for Fig 2 (homog PP)
% calls homogPoissonP.m -- freq=PP & Severity is 1 of 5 dists
% here using Burr (see run_homogPP.m and drv_homogPoissP.m)

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

flName='dInhomPP_burr';
severDistr=4; 
alph=(20:10:100)';
cpar=[2; 5; 2; 5];
kpar=[1.5; 1.5; 2; 2]; %making c & k same size

% a=100; c=3; k=1.5; 
% mn_sev=k*a*beta(k-1/c,1+1/c); sec_sev=k*a^2*beta(k-2/c,1+2/c);
% var_sev=sec_sev-mn_sev^2;

%outputs 
mnR=zeros(length(alph),length(cpar));
vrR=zeros(length(alph),length(cpar));
mnTw=zeros(length(alph),length(cpar));
vrTw=zeros(length(alph),length(cpar));
mnR_an=zeros(length(alph),length(cpar));
vrR_an=zeros(length(alph),length(cpar));
mnTw_an=zeros(length(alph),length(cpar));
vrTw_an=zeros(length(alph),length(cpar));

for kin=1:length(cpar)
tic
    for jin=1:length(alph)
    
        params=[ alph(jin) ; cpar(kin); kpar(kin)]; 
        
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
    save([pwd,'/',flName],'mn*','vr*','alph','cpar','kpar','lgamm','ap','tau_d','Tw');
end



