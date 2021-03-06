function [NS,ANSE,KGE,NS_radQ,X,WW,Qsim]=cal_MISDc_WEB_2L_snow_IE_h(data,X_ini,Ab,dirout,name,name_suff,FIG)
NPAR=10;
if nargin<3,X_ini=ones(NPAR,1)*.5;end
[RES]=fmincon(@calibOK,X_ini,[],[],[],[],...
     zeros(NPAR,1),ones(NPAR,1),[],optimset('Display','iter','MaxIter',200,...
     'MaxFunEvals',1200,'TolFun',1E-5,'TolCon',6,...
     'Largescale','off','Algorithm','active-set'),data,Ab,dirout,name,name_suff);
X=convert_adim(RES);
[NS,ANSE,KGE,NS_radQ,Qsim,WW]= MISDc_WEB_2L_snow_IE_h(data,X,Ab,FIG,dirout,name,name_suff);

%---------------------------------------------------------------------------------
function [err]=calibOK(X_0,data,Ab,dirout,name,name_suff)

X=convert_adim(X_0);
[NS,ANSE,KGE,NS_radQ]=MISDc_WEB_2L_snow_IE_h(data,X,Ab,0,dirout,name,name_suff);
err=1-KGE;

%---------------------------------------------------------------------------------
function X=convert_adim(X_0)
%      W_p  W_max2    m2    Ks  gamma1   Kc   alpha      Cm   m22   Ks2     
LOW=[  0.1,    100,  5.0, 0.10,    0.5, 0.4,    1.0, 0.1/24,  2.0, 0.01]';
UP =[  0.9,   1000, 60.0, 40.0,    3.5, 2.0,   35.0,      3, 40.0, 65.0]';
X=LOW+(UP-LOW).*X_0;
