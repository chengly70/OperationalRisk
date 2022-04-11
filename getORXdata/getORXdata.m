%script to import data from 'ORX dataA.xlsx', As of Nov 6, 2021

fileNameF='ORXdataFrq.xlsx';
fileNameS='ORXdataSev.xlsx';

sheetFrq=1;
sheetSev=1;

rngFrq='B2:K7';
rngSev='B2:K7';

%Internal Fraud, External Fraud, Employment Practices & Workplace Safety, 
% Clients, Products, & Business Practices (CPBP), Disasters & Public Safety, (DPS)
% Technology & Infrastructure Failure (TIF), 
% Executation, Delivery & Process Management (EDPM)
riskCats={'IF','Internal Fraud';'EF','External Fraud';'EPWS','Employment Practices & Workplace Safety'; ...
    'CPBP', 'Clients, Products, & Business Practices';
    'DPS','Disasters & Public Safety';'TIF','Technology & Infrastructure Failure'; 'EDPM','Executation, Delivery & Process Management'};

matDataF=xlsread(fileNameF,sheetFrq,rngFrq);
matDataS=xlsread(fileNameS,sheetSev,rngSev);

TotGrss=matDataS(:,1); %only getting bank losses
TotFreq=matDataF(:,1); 
numBnks=matDataF(:,2);  %same as matDataS(:,2)

%freq of all 7 risk categories, 1st col=IF, 2nd col=EF, etc
freqAll_rCats=matDataF(:,4:end); %get cols 8, 10, 12, etc.
%divide by numBnks to get per bank average (# events per year per bank):
freqAll_rCats=diag(1./numBnks)*freqAll_rCats;

%severity of all 7 risk categories, in euro-Billions
sevAll_rCats=matDataS(:,4:end); %get cols 7, 9, 11, etc.
%divide by numBnks to get per bank average, multl by 1000 to get Millions
%  so units are (millions-euro per year per bank):
sevAll_rCats=diag(1./numBnks)*sevAll_rCats*1000;

covM_freq=cov(freqAll_rCats); %7x7 cov matrix of frequencies
covM_sev=cov(sevAll_rCats); %7x7 cov matrix of severities

crrM_freq=corrcoef(freqAll_rCats); %7x7 correl matrix of frequencies
crrM_sev=corrcoef(sevAll_rCats); %7x7 correl matrix of frequencies

%avg freq for all 7 risk cat (avg over years)
avgFreq=mean(freqAll_rCats)';
%avg sev for all 7 risk cat (avg over years)
avgSev=mean(sevAll_rCats)';

