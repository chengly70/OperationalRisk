function [mean_R,var_R,mean_Tw,var_Tw,an_mean,an_var,an_meanTw,an_varTw]=...
    InhomogPoissonP_mr(dt,tEnd,ap,lgamm,tau_d,Tw,severDistr,params,numrel)
% Aggregates all of the pieces from newTwin_logn.m, allowing general
% several distr, but assuming freq distr ~ Poisson Process
% INPUTS: tEnd=total time (in years), ap=jump,tau_d=timescale, lgamm=rate of PP(freq) [#/years]
% Tw=larger time window, .. numrel=# of realizatins
% VERY MUCH like InhomogPoissonP.m EXCEPT allow many realiz at once to avg.
% severDistr=indicates 1 of 5 heavy-tailed severity, with:
%   severDistr=0 or as string='gamm', params(1)=alph, params(2)=bet
%   severDistr=1 or as string='logn', params(1)=mu, params(2)=sig
%   severDistr=2 or as string='gpd', params(1)=mu, params(2)=sig
%   severDistr=3 or as string='wbl', params(1)=kpar, params(2)=sig (thet=0)
%   severDistr=4 or as string='burr', params(1)=alph, params(2)=cpar,params(3)=kpar
% OUTPUTS: tmvc=time vector (dt) to tEnd, mean_R/var_R for time series (dt)
% mean_Tw/var_Tw for cumul losses in Tw;
% an_mean/var an_meanTw/varTw are analytic calc of mean/vars in dt & Tw
% similar to homogPoissonP.m BUT the freq distriubtion follows ODE with delta kicks

rng('shuffle') %seed random number generator

% -- OUTPUTS (initialize) ---
mean_R=0;   %time avg in dt window
var_R=0;    %var in time in dt window
mean_Tw=0;  %mean of cumul loss in Tw>dt
var_Tw=0;

%tmvc=(0:dt:tEnd)';
Lt=round(tEnd/dt)+1; %total number of time steps

% -- CHECKS ---
if( Tw<=dt ) disp(['Make Tw > dt=',num2str(dt)]); return; end
if( rem(Tw,dt)~= 0 ) disp(['dt=',num2str(dt),' doesnt divide Tw!']); return; end
if( lgamm*dt > 0.3 ) disp('Either mean nu is too small or dt too large, change!'); return; end

%instance of severity; draw a head of time (depends on severDistr)
numSevs=numrel*poissinv(.999999,tEnd*ap*tau_d*lgamm); %produce numSevs severity amounts; 99.99% of time will be enough
severInstance=zeros(numSevs,1); 

switch severDistr
    case {0,'gamm'} %gamma (NOT heavy-tailed)
        if(length(params)~=2) disp('Gamma has 2 pars, alph=params(1), bet=params(2)'); return; end
        alph=params(1); bet=params(2);
        if(alph<=0 || bet<=0) disp('Both parms has to be >0'); return; end
        mn_sev=alph*bet; var_sev=alph*bet^2; sec_sev=var_sev+mn_sev^2;
        
        severInstance=gamrnd(alph,bet,numSevs,1); %get severities upfront
        
    case {1,'logn'} %lognormal
        if(length(params)~=2) disp('Logn has 2 pars, mu=params(1), sig=params(2)'); return; end
        mu=params(1); sig=params(2); %set logn parameters
        if(sig<=0) disp('sig needs to be >0'); return; end
        mn_sev = exp(mu+(sig^2)/2);
        var_sev= (exp(sig^2)-1)*(exp(2*mu+sig^2));
        sec_sev= var_sev+mn_sev^2;
        
        severInstance=lognrnd(mu,sig,numSevs,1);
        
    case {2,'gpd'} %GPD, generalized pareto distribution
        if(length(params)~=2) disp('GPD has 2 pars, kpar=params(1), sig=params(2)'); return; end
        kpar=params(1); sig=params(2); %set GPD params
        if(kpar>1/2) disp('VarS is infin, change kpar=params(1)<1/2'); return; end
        if(kpar>1) disp('MeanS is infin, change kpar=params(1)<1'); return; end
        if(kpar==0) %exp(1)
            mn_sev=sig; var_sev=sig^2; sec_sev=var_sev+mn_sev^2;
        else
            mn_sev=sig/(1-kpar); var_sev=sig^2/((1-kpar)^2*(1-2*kpar)); 
            sec_sev=var_sev+mn_sev^2;
        end
        
        severInstance=gprnd(kpar,sig,0,numSevs,1);
        
    case {3,'wbl'} %Weibull
        if(length(params)~=2) disp('Weibull has 2 pars, apar=params(1), bpar=params(2)'); return; end
        apar=params(1); bpar=params(2); %scale & shape, resp.
        if(apar<=0 || bpar<=0) disp('Both parms needs to be >0'); return; end
        mn_sev=apar*gamma(1+1/bpar); sec_sev=apar^2*gamma(1+2/bpar); var_sev=sec_sev-mn_sev^2;
        
        severInstance=wblrnd(apar,bpar,numSevs,1);
        
    case {4,'burr'} %Burr
        if(length(params)~=3) disp('Burr has 3 pars, alph=params(1), cpar=params(2), kpar=params(3)'); return; end
        alph=params(1); cpar=params(2); kpar=params(3);
        if(alph<=0 || cpar<=0 || kpar<=0) disp('All parms needs to be >0'); return; end
        if(kpar<1/cpar) disp('Need k larger or c smaller, cant have k<1/c (mean inf)'); return; end
        if(kpar<2/cpar) disp('Need k larger or c smaller, cant have k<2/c (var inf)'); return; end
        mn_sev=kpar*alph*beta(kpar-1/cpar,1+1/cpar); sec_sev=kpar*alph^2*beta(kpar-2/cpar,1+2/cpar);
        var_sev=sec_sev-mn_sev^2;
        
        severInstance=random('burr',alph,cpar,kpar,numSevs,1);
end

nu_z=zeros(numrel,Lt);

mean_R=0;
var_R=0;

mean_Tw=0;
var_Tw=0;

t_row = round(Tw/dt); %# dt steps in each Tw (time window)
num_Tw = floor(Lt/t_row);    %tot. # bins of size Tw; throw out last one
        
LossTmSeries=zeros(numrel,Lt); %stores severity of losses (2nd col) (time is by index location)

cnt=1;
for j=2:Lt
        nu_z(:,j)=nu_z(:,j-1)+dt/tau_d*(-nu_z(:,j-1));
        xiv=rand(numrel,1);
        nu_z(xiv < dt*lgamm,j)=nu_z(xiv < dt*lgamm,j)+ap; %Poisson Process; jumps
        
        xiv=rand(numrel,1);
        siz=sum(xiv < dt*nu_z(:,j));
        if(siz >0)
        LossTmSeries(xiv < dt*nu_z(:,j),j) = severInstance(cnt:cnt+siz-1);
        cnt=cnt+siz; %keeps track of pre-drawn Severities
        end
end

lssReshap=sum(reshape(LossTmSeries(:,1:t_row*num_Tw)',t_row,num_Tw*numrel));
Loss_Twin=reshape(lssReshap,num_Tw,numrel)'; %vector of size numrel x num_Tw

%save the simulation results
mean_R=mean(mean(LossTmSeries,2)); %average over time, mean(numrel x 1)
var_R=mean(var(LossTmSeries'));    %pt-wise var in time, mean(numrel x 1)

mean_Tw=mean(mean(Loss_Twin,2));
var_Tw=mean(var(Loss_Twin'));


% -- Analytic (mathematical) calculations --
an_mean= mn_sev*(ap*lgamm*tau_d)*dt;
an_var=sec_sev*(ap*lgamm*tau_d*dt) - (mn_sev*lgamm*tau_d*ap*dt)^2;

an_meanTw= Tw*an_mean/dt;
an_varTw=2*tau_d*an_var*(Tw+tau_d*(exp(-Tw/tau_d)-1))/dt;


end