function [NS,ANSE,KGE,NS_radQ,X,WW,Qsim]=cal_MISDc_WEB_2L_snow_IE(data,X_ini,Ab,dir,name,name_suff,FIG)
NPAR=10;
if nargin<3,X_ini=ones(NPAR,1)*.5;end
[RES]=fmincon(@calibOK,X_ini,[],[],[],[],...
     zeros(NPAR,1),ones(NPAR,1),[],optimset('Display','iter','MaxIter',150,...
     'MaxFunEvals',500,'TolFun',1E-5,'TolCon',6,...
     'Largescale','off','Algorithm','active-set'),data,Ab,dir,name,name_suff);
X=convert_adim(RES);
[NS,ANSE,KGE,NS_radQ,Qsim,WW]= MISDc_WEB_2L_snow_IE(data,X,Ab,FIG,dir,name,name_suff);

%---------------------------------------------------------------------------------
function [err]=calibOK(X_0,data,Ab,dir,name,name_suff)

X=convert_adim(X_0);
[NS,ANSE,KGE,NS_radQ]=MISDc_WEB_2L_snow_IE(data,X,Ab,0,dir,name,name_suff);
err=1-KGE;
save X_PAR

%---------------------------------------------------------------------------------
function X=convert_adim(X_0)
%      W_p  W_max2    m2    Ks  gamma1   Kc   alpha      Cm   m22   Ks2     
LOW=[  0.1,    300,  2.0, 0.10,    0.5, 0.4,    1.0, 0.1/24,  5.0, 0.01]';
UP =[  0.9,   4000, 10.0, 20.0,    3.5, 2.0,   15.0,      3, 20.0, 65.0]';
X=LOW+(UP-LOW).*X_0;
