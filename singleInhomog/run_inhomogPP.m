% script to get data for Fig 4 (inhomog PP)
% calls InhomogPoissonP.m -- freq=PP & Severity is 1 of 5 dists
% here using gamma. VERY similar to run_imhomogPP_1realz.m BUT 
% adding more realiz & less overall tEnd

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

flName='dInhomPP_gamma';
severDistr=0; 
bet=[2;3;4];
alph=(3:2:25)';

%outputs 
mnR=zeros(length(alph),length(bet));
vrR=zeros(length(alph),length(bet));
mnTw=zeros(length(alph),length(bet));
vrTw=zeros(length(alph),length(bet));
mnR_an=zeros(length(alph),length(bet));
vrR_an=zeros(length(alph),length(bet));
mnTw_an=zeros(length(alph),length(bet));
vrTw_an=zeros(length(alph),length(bet));

for kin=1:length(bet)
tic
    for jin=1:length(alph)
    
        params=[ alph(jin); bet(kin)];
        
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
    save([pwd,'/',flName],'mn*','vr*','alph','bet','lgamm','ap','tau_d','Tw');
end



