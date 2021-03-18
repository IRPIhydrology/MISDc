%% MISDc 2 layer daily
clc,clear,close all
inputdata=load('migi_0205_daily.txt');
area=137; % basin area
FIG=1;
NPAR=10;
X_ini=ones(NPAR,1)*.1;X_ini(1)=0.05;
[NS,ANSE,KGE,NSradQ,X,WW,Qsim]=...
    cal_MILc_WEB_2L_snow_IE(inputdata,...
    X_ini,area,cd,['Qsim_MISDc_2L_daily'],['_test'],FIG);
        
%% MISDc 2 layer hourly
clc,clear,close all
inputdata=load('migi_0205.txt');
area=137; % basin area
FIG=1;
NPAR=10;
X_ini=ones(NPAR,1)*.1;X_ini(1)=0.05;

[NS,ANSE,KGE,X,WW,Qsim]=...
    cal_MILc_WEB_2L_snow_IE_h(inputdata,...
    X_ini,area,cd,['Qsim_MISDc_2L_hourly'],['_test'],FIG);
