%driver script that calls InhomTwoPP_mr.m -- freq=PP & Severity is 1 of 5 dists
% here varying SEVS and cinp (!!NOT done!!, specify 40 sevs parms)

flName='dVary_CseverDist';

tEnd=10000; %time in years
numrel=2000;
dt=0.005; %less than a day
tmvc=(0:dt:tEnd)';
Lt=round(tEnd/dt)+1;

Tw=1; %usual time window 

c_inp=[-.8; -.6; -.4; -.2; .2; .4; .6; .8]; %correl of 2 frequencies; has to be in [-1,1]
lenCin=length(c_inp);

lenSpms=40;
hlfSmp=lenSpms/2; %needs to divide
%logn*gpd for first  , wbl*burr
Sdist_id = [ones(hlfSmp,1) 2*ones(hlfSmp,1);  3*ones(hlfSmp,1) 4*ones(hlfSmp,1)];
%Sparms = zeros(lenSpms,5);  %logn*gpd for first  , wbl*burr
Sparms = [rand(hlfSmp,1)*.487+1.177 rand(hlfSmp,1)*2.1+1.1 rand(hlfSmp,1)*.1+.1 51+rand(hlfSmp,1)*42.5 zeros(hlfSmp,1); ...
    3+rand(hlfSmp,1)*5 .25+rand(hlfSmp,1)*.45 20+rand(hlfSmp,1)*80 2+rand(hlfSmp,1)*3 1.5+rand(hlfSmp,1)*.5];

%severDistr=0;  %0,gamm=gamma; 1,logn=lognormal; 2,gpd=GPD; 3,wbl=Weibull; 4,burr=Burr
%severDistr='gamm'; %can also use char
%params=[ 2; 1]; %gamm[alph;bet], logn[mu;sig], gpd[kpar,sig], wbl[a,b], burr[alph,cpar,kpar]
% severDistr=1; params=[1; .5];
% severDistr=4;  params=[1; 2; 2]; %

%outputs 
mnR=zeros(lenSpms,lenCin,2);
vrR=zeros(lenSpms,lenCin,2);
mnTw=zeros(lenSpms,lenCin,2);
vrTw=zeros(lenSpms,lenCin,2);
mnR_an=zeros(lenSpms,lenCin,2);
vrR_an=zeros(lenSpms,lenCin,2);
mnTw_an=zeros(lenSpms,lenCin,2);
vrTw_an=zeros(lenSpms,lenCin,2);
covR=zeros(lenSpms,lenCin);
covR_an=zeros(lenSpms,lenCin);
covTw=zeros(lenSpms,lenCin);
covTw_an=zeros(lenSpms,lenCin);

for kin=1:lenCin
    tic
    for jin=1:lenSpms
        
        paramsL1=struct('a',1.5,'nu_i',30,'tau_d',1.3,'SeverDist',Sdist_id(jin,1),'Sparams',Sparms(jin,1:2));
        if(jin>hlfSmp) %2nd distr is Burr, need 3 parms
            paramsL2=struct('a',2,'nu_i',40,'tau_d',.75,'SeverDist',Sdist_id(jin,2),'Sparams',Sparms(jin,3:5));
        else %not Burr
            paramsL2=struct('a',2,'nu_i',40,'tau_d',.75,'SeverDist',Sdist_id(jin,2),'Sparams',Sparms(jin,3:4));
        end
        
        [mean_R,var_R,mean_Tw,var_Tw,an_mean,an_var,an_meanTw,an_varTw,...
            cov_R,an_cov,cov_Tw,an_covTw]=InhomTwoPP_mr(dt,tEnd,Tw,c_inp(kin),paramsL1,paramsL2,numrel);
        %save data
        mnR(jin,kin,:)=mean_R';
        vrR(jin,kin,:)=var_R';
        mnTw(jin,kin,:)=mean_Tw';
        vrTw(jin,kin,:)=var_Tw';
        mnR_an(jin,kin,:)=an_mean';
        vrR_an(jin,kin,:)=an_var';
        mnTw_an(jin,kin,:)=an_meanTw';
        vrTw_an(jin,kin,:)=an_varTw';
        covR(jin,kin)=cov_R;
        covR_an(jin,kin)=an_cov;
        covTw(jin,kin)=cov_Tw;
        covTw_an(jin,kin)=an_covTw;
        
        %save mat file
        save([pwd,'/',flName],'mn*','vr*','covR*','covTw*','Tw','c_inp','Sparms','Sdist_id');
    end
    toc
end


