function [mean_R,var_R,mean_Tw,var_Tw,an_mean,an_var,an_meanTw,an_varTw,...
    cov_R,an_cov,cov_Tw,an_covTw]=InhomTwoPP_mr(dt,tEnd,Tw,c_inp,paramsL1,paramsL2,numrel)
% Aggregates all of the pieces from newTwin_logn.m, allowing general
% several distr, but assuming freq distr ~ Poisson Process
% INPUTS: tEnd=total time (in years), Tw=larger time window, 
% paramsL1, paramsL2
% severDistr=indicates 1 of 5 heavy-tailed severity, with:
%   severDistr=0 or as string='gamm', params(1)=alph, params(2)=bet
%   severDistr=1 or as string='logn', params(1)=mu, params(2)=sig
%   severDistr=2 or as string='gpd', params(1)=mu, params(2)=sig
%   severDistr=3 or as string='wbl', params(1)=kpar, params(2)=sig (thet=0)
%   severDistr=4 or as string='burr', params(1)=alph, params(2)=cpar,params(3)=kpar
% OUTPUTS: tmvc=time vector (dt) to tEnd, mean_R/var_R for time series (dt)
% mean_Tw/var_Tw for cumul losses in Tw;
% an_mean/var an_meanTw/varTw are analytic calc of mean/vars in dt & Tw

rng('shuffle') %seed random number generator

%getting parameters from the structures (inputs)
a1=paramsL1.a;
gamma1=paramsL1.nu_i;
tau_d1=paramsL1.tau_d;
params1=paramsL1.Sparams;
severDistr1=paramsL1.SeverDist;

a2=paramsL2.a;
gamma2=paramsL2.nu_i;
tau_d2=paramsL2.tau_d;
params2=paramsL2.Sparams;
severDistr2=paramsL2.SeverDist;

gamBar=min(gamma1,gamma2)*c_inp;
nu1aug=gamma1-gamBar;
nu2aug=gamma2-gamBar;
prb1=dt*nu1aug;
prb2=dt*nu2aug; %allow all 3 events at a time step
prb3=dt*gamBar;

% -- OUTPUTS (initialize) ---
mean_R=zeros(2,1);   %time avg in dt window
var_R=zeros(2,1);    %var in time in dt window
mean_Tw=zeros(2,1);  %mean of cumul loss in Tw>dt
var_Tw=zeros(2,1);
an_mean=zeros(2,1); an_meanTw=zeros(2,1); an_var=zeros(2,1); an_varTw=zeros(2,1); 
cov_R=0; an_cov=0; cov_Tw=0; an_covTw=0;

%tmvc=(0:dt:tEnd)';
Lt=round(tEnd/dt)+1; %total number of time steps
nu_1=zeros(Lt,1);
nu_2=zeros(Lt,1);

% -- CHECKS ---
if( Tw<=dt ) disp(['Make Tw > dt=',num2str(dt)]); return; end
if( rem(Tw,dt)~= 0 ) disp(['dt=',num2str(dt),' doesnt divide Tw!']); return; end

%instance of severity; draw a head of time (depends on severDistr)
numSevs=numrel*poissinv(.999,tEnd*gamma1*a1*tau_d1); %produce numSevs severity amounts; 99.9% of time will be enough
severInstance1=zeros(numSevs,1); 
switch severDistr1
    case {0,'gamm'} %gamma (NOT heavy-tailed)
        if(length(params1)~=2) disp('Gamma has 2 pars, alph=params(1), bet=params(2)'); return; end
        alph=params1(1); bet=params1(2);
        if(alph<=0 || bet<=0) disp('Both parms has to be >0'); return; end
        mn_sev=alph*bet; var_sev=alph*bet^2; sec_sev=var_sev+mn_sev^2;
        
        severInstance1=gamrnd(alph,bet,numSevs,1); %get severities upfront
        
    case {1,'logn'} %lognormal
        if(length(params1)~=2) disp('Logn has 2 pars, mu=params(1), sig=params(2)'); return; end
        mu=params1(1); sig=params1(2); %set logn parameters
        if(sig<=0) disp('sig needs to be >0'); return; end
        mn_sev = exp(mu+(sig^2)/2);
        var_sev= (exp(sig^2)-1)*(exp(2*mu+sig^2));
        sec_sev= var_sev+mn_sev^2;
        
        severInstance1=lognrnd(mu,sig,numSevs,1);
        
    case {2,'gpd'} %GPD, generalized pareto distribution
        if(length(params1)~=2) disp('GPD has 2 pars, kpar=params(1), sig=params(2)'); return; end
        kpar=params1(1); sig=params1(2); %set GPD params
        if(kpar>1/2) disp('VarS is infin, change kpar=params(1)<1/2'); return; end
        if(kpar>1) disp('MeanS is infin, change kpar=params(1)<1'); return; end
        if(kpar==0) %exp(1)
            mn_sev=sig; var_sev=sig^2; sec_sev=var_sev+mn_sev^2;
        else
            mn_sev=sig/(1-kpar); var_sev=sig^2/((1-kpar)^2*(1-2*kpar)); 
            sec_sev=var_sev+mn_sev^2;
        end
        
        severInstance1=gprnd(kpar,sig,0,numSevs,1);
        
    case {3,'wbl'} %Weibull
        if(length(params1)~=2) disp('Weibull has 2 pars, apar=params(1), bpar=params(2)'); return; end
        apar=params1(1); bpar=params1(2); %scale & shape, resp.
        if(apar<=0 || bpar<=0) disp('Both parms needs to be >0'); return; end
        mn_sev=apar*gamma(1+1/bpar); sec_sev=apar^2*gamma(1+2/bpar); var_sev=sec_sev-mn_sev^2;
        
        severInstance1=wblrnd(apar,bpar,numSevs,1);
        
    case {4,'burr'} %Burr
        if(length(params1)~=3) disp('Burr has 3 pars, alph=params(1), cpar=params(2), kpar=params(3)'); return; end
        alph=params1(1); cpar=params1(2); kpar=params1(3);
        if(alph<=0 || cpar<=0 || kpar<=0) disp('All parms needs to be >0'); return; end
        if(kpar<1/cpar) disp('Need k larger or c smaller, cant have k<1/c (mean inf)'); return; end
        if(kpar<2/cpar) disp('Need k larger or c smaller, cant have k<2/c (var inf)'); return; end
        mn_sev=kpar*alph*beta(kpar-1/cpar,1+1/cpar); sec_sev=kpar*alph^2*beta(kpar-2/cpar,1+2/cpar);
        var_sev=sec_sev-mn_sev^2;
        
        severInstance1=random('burr',alph,cpar,kpar,numSevs,1);
end
numSevs2=numrel*poissinv(.999,tEnd*gamma2*a2*tau_d2); %produce numSevs severity amounts; 99.9% of time will be enough
severInstance2=zeros(numSevs2,1); 
switch severDistr2
    case {0,'gamm'} %gamma (NOT heavy-tailed)
        if(length(params2)~=2) disp('Gamma has 2 pars, alph=params(1), bet=params(2)'); return; end
        alph=params2(1); bet=params2(2);
        if(alph<=0 || bet<=0) disp('Both parms has to be >0'); return; end
        mn_sev2=alph*bet; var_sev2=alph*bet^2; sec_sev2=var_sev2+mn_sev2^2;
        
        severInstance2=gamrnd(alph,bet,numSevs2,1); %get severities upfront
        
    case {1,'logn'} %lognormal
        if(length(params2)~=2) disp('Logn has 2 pars, mu=params(1), sig=params(2)'); return; end
        mu=params2(1); sig=params2(2); %set logn parameters
        if(sig<=0) disp('sig needs to be >0'); return; end
        mn_sev2 = exp(mu+(sig^2)/2);
        var_sev2= (exp(sig^2)-1)*(exp(2*mu+sig^2));
        sec_sev2= var_sev2+mn_sev2^2;
        
        severInstance2=lognrnd(mu,sig,numSevs2,1);
        
    case {2,'gpd'} %GPD, generalized pareto distribution
        if(length(params2)~=2) disp('GPD has 2 pars, kpar=params(1), sig=params(2)'); return; end
        kpar=params2(1); sig=params2(2); %set GPD params
        if(kpar>1/2) disp('VarS is infin, change kpar=params(1)<1/2'); return; end
        if(kpar>1) disp('MeanS is infin, change kpar=params(1)<1'); return; end
        if(kpar==0) %exp(1)
            mn_sev2=sig; var_sev2=sig^2; sec_sev2=var_sev2+mn_sev2^2;
        else
            mn_sev2=sig/(1-kpar); var_sev2=sig^2/((1-kpar)^2*(1-2*kpar)); 
            sec_sev2=var_sev2+mn_sev2^2;
        end
        
        severInstance2=gprnd(kpar,sig,0,numSevs2,1);
        
    case {3,'wbl'} %Weibull
        if(length(params2)~=2) disp('Weibull has 2 pars, apar=params(1), bpar=params(2)'); return; end
        apar=params2(1); bpar=params2(2); %scale & shape, resp.
        if(apar<=0 || bpar<=0) disp('Both parms needs to be >0'); return; end
        mn_sev2=apar*gamma(1+1/bpar); sec_sev2=apar^2*gamma(1+2/bpar); var_sev2=sec_sev2-mn_sev2^2;
        
        severInstance2=wblrnd(apar,bpar,numSevs2,1);
        
    case {4,'burr'} %Burr
        if(length(params2)~=3) disp('Burr has 3 pars, alph=params(1), cpar=params(2), kpar=params(3)'); return; end
        alph=params2(1); cpar=params2(2); kpar=params2(3);
        if(alph<=0 || cpar<=0 || kpar<=0) disp('All parms needs to be >0'); return; end
        if(kpar<1/cpar) disp('Need k larger or c smaller, cant have k<1/c (mean inf)'); return; end
        if(kpar<2/cpar) disp('Need k larger or c smaller, cant have k<2/c (var inf)'); return; end
        mn_sev2=kpar*alph*beta(kpar-1/cpar,1+1/cpar); sec_sev2=kpar*alph^2*beta(kpar-2/cpar,1+2/cpar);
        var_sev2=sec_sev2-mn_sev2^2;
        
        severInstance2=random('burr',alph,cpar,kpar,numSevs2,1);
end

t_row = round(Tw/dt); %# dt steps in each Tw (time window)
num_Tw = floor(Lt/t_row);    %tot. # bins of size Tw; throw out last one

LossTmSeries=zeros(numrel,Lt); %stores severity of losses (2nd col) (time is by index location)
LossTmSeries2=zeros(numrel,Lt);
Avloss=zeros(Lt,1);
Avloss2=zeros(Lt,1);
Loss_Twin=zeros(numrel,num_Tw);
Loss_Twin2=zeros(numrel,num_Tw);

cnt1=1;
cnt2=1;
if(c_inp>=0)
    %loop for each time series
    for j=2:Lt
        %sim nu_1/nu_2 FIRST
        nu_1(j)=nu_1(j-1)+dt/tau_d1*(-nu_1(j-1));
        if(rand < prb1)
            nu_1(j)=nu_1(j)+a1; %Poisson Process; jumps
        end
        
        nu_2(j)=nu_2(j-1)+dt/tau_d2*(-nu_2(j-1));
        if(rand < prb2)
            nu_2(j)=nu_2(j)+a2; %Poisson Process; jumps
        end
        
        if(rand<prb3) %common JUMP prob
            nu_1(j)=nu_1(j)+a1;
            nu_2(j)=nu_2(j)+a2;
        end
        
        % now sim LossTmSeries
        xiv=rand(numrel,1);
        siz=sum(xiv < dt*nu_1(j));
        if(siz >0)
            if(cnt1+siz-1 > numSevs) % get more Severities
                severInstance1=moreSevs(paramsL1,numSevs);
                cnt1=1;
                LossTmSeries(xiv < dt*nu_1(j),j) = severInstance1(cnt1:cnt1+siz-1);
            else
                LossTmSeries(xiv < dt*nu_1(j),j) = severInstance1(cnt1:cnt1+siz-1);
            end       
            cnt1=cnt1+siz; %keeps track of pre-drawn Severities
        end
        xiv=rand(numrel,1);
        siz=sum(xiv < dt*nu_2(j));
        if(siz >0)
            if(cnt2+siz-1 > numSevs2)
                severInstance2=moreSevs(paramsL2,numSevs2);
                cnt2=1;
                LossTmSeries2(xiv < dt*nu_2(j),j) = severInstance2(cnt2:cnt2+siz-1);
            else
                LossTmSeries2(xiv < dt*nu_2(j),j) = severInstance2(cnt2:cnt2+siz-1);
            end
            cnt2=cnt2+siz; %keeps track of pre-drawn Severities
        end
    end
else
    prb3=dt*abs(gamBar);
    prb1=dt*nu1aug;
    prb2=prb1+dt*nu2aug; %allow all 2 events at a time step
    %loop for each time series
    for j=2:Lt
        %sim nu_1/nu_2 FIRST
        nu_1(j)=nu_1(j-1) + dt/tau_d1*(-nu_1(j-1));
        nu_2(j)=nu_2(j-1) + dt/tau_d2*(-nu_2(j-1));
        xiv1=rand;   xivc=rand;
        if(xivc < prb3) %allow only 1 event at a time
            if(rand<=.5) %half go in 1, other in 2
                nu_1(j)=nu_1(j)-min(a1,nu_1(j));
                nu_2(j)=nu_2(j)+a2;
            else
                nu_1(j)=nu_1(j)+a1;
                nu_2(j)=nu_2(j)-min(a2,nu_2(j));
            end
        else
            if(xiv1 < prb1)
                nu_1(j)=nu_1(j)+a1;
            elseif(xiv1 < prb2)
                nu_2(j)=nu_2(j)+a2;
            end
        end
        
        % now sim LossTmSeries
        xiv=rand(numrel,1);
        siz=sum(xiv < dt*nu_1(j));
        if(siz >0)
            if(cnt1+siz-1 > numSevs) % get more Severities
                severInstance1=moreSevs(paramsL1,numSevs);
                cnt1=1;
                LossTmSeries(xiv < dt*nu_1(j),j) = severInstance1(cnt1:cnt1+siz-1);
            else
                LossTmSeries(xiv < dt*nu_1(j),j) = severInstance1(cnt1:cnt1+siz-1);
            end       
            cnt1=cnt1+siz; %keeps track of pre-drawn Severities
        end
        xiv=rand(numrel,1);
        siz=sum(xiv < dt*nu_2(j));
        if(siz >0)
            if(cnt2+siz-1 > numSevs2)
                severInstance2=moreSevs(paramsL2,numSevs2);
                cnt2=1;
                LossTmSeries2(xiv < dt*nu_2(j),j) = severInstance2(cnt2:cnt2+siz-1);
            else
                LossTmSeries2(xiv < dt*nu_2(j),j) = severInstance2(cnt2:cnt2+siz-1);
            end
            cnt2=cnt2+siz; %keeps track of pre-drawn Severities
        end
    end
end

%-get cumulative losses in window Tw
lssReshap=sum(reshape(LossTmSeries(:,1:t_row*num_Tw)',t_row,num_Tw*numrel));
Loss_Twin=reshape(lssReshap,num_Tw,numrel)'; %vector of size numrel x num_Tw
lssReshap2=sum(reshape(LossTmSeries2(:,1:t_row*num_Tw)',t_row,num_Tw*numrel));
Loss_Twin2=reshape(lssReshap2,num_Tw,numrel)'; %vector of size numrel x num_Tw

%save the simulation results
mean_R(1,1)=mean(mean(LossTmSeries,2)); %average over time, numrel x 1
mean_R(2,1)=mean(mean(LossTmSeries2,2)); %average over time, numrel x 1
var_R(1,1)=mean(var(LossTmSeries'));    %pt-wise var in time, numrel x 1
var_R(2,1)=mean(var(LossTmSeries2'));    %pt-wise var in time, numrel x 1

mean_Tw(1)=mean(mean(Loss_Twin,2));
mean_Tw(2)=mean(mean(Loss_Twin2,2));
var_Tw(1)=mean(var(Loss_Twin')');
var_Tw(2)=mean(var(Loss_Twin2')');

Avloss = mean(LossTmSeries)';%average over realizations(columns)
Avloss2 = mean(LossTmSeries2)';%average over realizations(columns)
Cc=xcorr(Avloss,Avloss2,round(tEnd*.05/dt),'unbiased');
Cc=Cc-mean(Avloss)*mean(Avloss2); %scale b/c makes peak=1; use avg of EACH realz' var
len_c=(length(Cc)-1)/2+1;

cov_R=Cc(len_c);

%cov in Tw from simulations, larger time windows

cov_Tw=zeros(numrel,1);
for j=1:numrel
    tmp=cov(Loss_Twin(j,:),Loss_Twin2(j,:));
    cov_Tw(j)=tmp(1,2)*dt; %only saving (1,2) entry
end
cov_Tw=mean(cov_Tw);

% -- Analytic (mathematical) calculations --
an_mean(1,1)=mn_sev*(a1*gamma1*tau_d1)*dt;
an_mean(2,1)=mn_sev2*(a2*gamma2*tau_d2)*dt;

an_var(1,1)=sec_sev*(a1*gamma1*tau_d1*dt) - (mn_sev*gamma1*tau_d1*a1*dt)^2;
an_var(2,1)=sec_sev2*(a2*gamma2*tau_d2*dt) - (mn_sev2*gamma2*tau_d2*a2*dt)^2;

an_meanTw=Tw*an_mean./dt;

an_varTw(1,1)=2*tau_d1*an_var(1,1)*(Tw+tau_d1*(exp(-Tw/tau_d1)-1))/dt;
an_varTw(2,1)=2*tau_d2*an_var(2,1)*(Tw+tau_d2*(exp(-Tw/tau_d2)-1))/dt;

an_cov=mn_sev*mn_sev2*gamBar*a1*a2*tau_d1*tau_d2/(tau_d1+tau_d2)*dt^2;%cvmt(1,2);
% Cc_theory=an_cov*exp(-abs(t_c2)./tau_d2);
% Cc_theory((length(t_c2)+1)/2:end)=...
%     an_cov*exp(-t_c2((length(t_c2)+1)/2:end)./tau_d1); %positive part, G_e1 is first

an_covTw = an_cov/dt*(tau_d2*(tau_d2*exp(-Tw/tau_d2)-tau_d2+Tw) + tau_d1*(tau_d1*exp(-Tw/tau_d1)-tau_d1+Tw));



% add sub-function in to draw more Severities
    function sevs=moreSevs(parmsL,numSevs)
        sevs=zeros(numSevs,1);
        switch parmsL.SeverDist
            case {0,'gamm'} %gamma (NOT heavy-tailed)
                sevs=gamrnd(parmsL.Sparams(1),parmsL.Sparams(2),numSevs,1); %get severities upfront
                
            case {1,'logn'} %lognormal
                sevs=lognrnd(parmsL.Sparams(1),parmsL.Sparams(2),numSevs,1);
                
            case {2,'gpd'} %GPD, generalized pareto distribution
                sevs=gprnd(parmsL.Sparams(1),parmsL.Sparams(2),0,numSevs,1);
                
            case {3,'wbl'} %Weibull
                sevs=wblrnd(parmsL.Sparams(1),parmsL.Sparams(2),numSevs,1);
                
            case {4,'burr'} %Burr
                sevs=random('burr',parmsL.Sparams(1),parmsL.Sparams(2),parmsL.Sparams(3),numSevs,1);
        end
    end

end