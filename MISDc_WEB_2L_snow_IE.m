%-----------------------------------------------------------------------------------------------
%   MISDc rainfall-runoff model (lumped version, daily)
%------------------------------------------------------------------------------------------------

function [NS,ANSE,KGE,NS_radQ,Qsim,WW,WW2]=MISDc_WEB_2L_snow_IE(DPEQ,PAR,Ab,FIG,dir,name,name_suff)

% Loading input data
[M,~]=size(DPEQ);
D=DPEQ(:,1); PIO_=DPEQ(:,2); TEMPER=DPEQ(:,3); Qobs=DPEQ(:,4);
delta_T=round(nanmean(diff(D))*24*10000)/10000;
MESE=month(D);

% Model parameter
W_p       = PAR(1);  % initial conditions, fraction of W_max (0-1)
W_max2    = PAR(2);  % total water capacity of 2nd layer
m2        = PAR(3);  % exponent of drainage for 1st layer
Ks        = PAR(4);  % hydraulic conductivity for 1st layer
gamma1    = PAR(5);  % coefficient lag-time relationship
Kc        = PAR(6);  % parameter of potential evapotranspiration
alpha     = PAR(7);  % exponent runoff
Cm        = PAR(8);  % Snow module parameter degree-day
m22       = PAR(9);  % exponent of drainage for 2nd layer
Ks2       = PAR(10); % hydraulic conductivity for 2nd layer

W_max     = 150;     % FIXED WATER CAPACITY 1st LAYER

% Other data
dt      = 0.2;          % computation time step in hour
Ks      = Ks.*delta_T;  % mm/h --> mm/delta_T
Ks2     = Ks2.*delta_T; % mm/h --> mm/delta_T

% Snow Module
[PIO,SWE]=snow_model(PIO_, TEMPER, -0.5, 0.5, Cm);

% Potential Evapotranspiration parameter
% L=[0.2100;0.2200;0.2300;0.2800;0.3000;0.3100;
%     0.3000;0.2900;0.2700;0.2500;0.2200;0.2000];
L=[0.2500;0.2500;0.2500;0.2500;0.2500;0.2500;
    0.2500;0.2500;0.2500;0.2500;0.2500;0.2500];
Ka=1.26;
EPOT=(TEMPER>0).*(Kc*(Ka*L(MESE).*(0.46*TEMPER+8)-2))/(24/delta_T);
% EPOT=Kc.*EPOT;
clear DPEQ TEMPER

% Initialization
BF=zeros(M,1);
QS=zeros(M,1);
WW=zeros(M,1);
WW2=zeros(M,1);

% Main ROUTINE
W=W_p*W_max;
W2=W_p*W_max2;
S=NaN;
Pcum=0;
IE=0;
for t=1:M
    IE=PIO(t)*((W/W_max).^alpha);
        
    E=EPOT(t)*W/W_max;
    if  W2<W_max2
        PERC=Ks*(W/W_max).^(m2);
    else
        PERC=0;
    end
    PERC2=Ks2*(W2/W_max2).^(m22);

    PERC(PERC>0.6.*W_max)=0.6.*W_max;
%     PERC2(PERC2>0.6.*W_max2)=0.6.*W_max2;
    
    W=max(0,W+(PIO(t)-IE-PERC-E)+SWE(t));
    W2=max(0,W2+PERC-PERC2);
    
    if W>=W_max
        SE=W-W_max;
        W=W_max;
    else
        SE=0;
    end
    if W2>=W_max2
        SE2=W2-W_max2;
        W2=W_max2;
    else
        SE2=0;
    end
    if W<0, W=0; end
    if W2<0, W2=0; end

    WW(t)=W./W_max;
    WW2(t)=W2./W_max2;
    
    % Runoff contribution
    BF(t)=Ks2*((W+W2)/(W_max+W_max2)).^(m22);
    QS(t)=IE+SE+SE2;       % surface flow
    
end

% Convolution (GIUH)
IUH1=IUH_comp(gamma1,Ab,dt,delta_T)*dt;IUH1=IUH1./sum(IUH1);
IUH2=IUH_NASH(1,0.5*gamma1,Ab,dt,delta_T)*dt;IUH2=IUH2./sum(IUH2);
QSint=interp1(1:M,QS,1:dt:M)';
BFint=interp1(1:M,BF,1:dt:M)';
temp1=conv(IUH1,QSint);
temp2=conv(IUH2,BFint);
Qsim1=temp2(1:round(1/dt):M*round(1/dt)).*(Ab.*1000./delta_T./3600);
Qsim=(temp1(1:round(1/dt):M*round(1/dt))+temp2(1:round(1/dt):M*round(1/dt)))...
    .*(Ab.*1000./delta_T./3600);
Qsim=real(Qsim);

% Calculation of model performance
RMSE=nanmean((Qsim-Qobs).^2).^0.5;
NS=1-nansum((Qsim-Qobs).^2)./nansum((Qobs-nanmean(Qobs)).^2);
ANSE=1-nansum((Qobs+nanmean(Qobs)).*(Qsim-Qobs).^2)./...
    nansum((Qobs+nanmean(Qobs)).*(Qobs-nanmean(Qobs)).^2);
NS_radQ=1-nansum((sqrt(Qsim)-sqrt(Qobs)).^2)./nansum((sqrt(Qobs)-nanmean(sqrt(Qobs))).^2);
NS_lnQ=1-nansum((log(Qsim+0.00001)-log(Qobs+0.00001)).^2)...
    ./nansum((log(Qobs+0.00001)-nanmean(log(Qobs+0.00001))).^2);
X=[Qsim,Qobs]; X(any(isnan(X)'),:) = [];
RRQ=corrcoef(X).^2; RQ=RRQ(2);
KGE=klinggupta(Qsim,Qobs);

% Figure
if FIG==1
    clf
    % saturation degree
    set(gcf,'paperpositionmode','manual','paperposition',[1 1 24 13])
    set(gcf,'position',[50   50   900   500])
    h(1) = axes('Position',[0.1 0.77 0.8 0.10]);
    hold on
    plot(D,WW,'r');
    plot(D,WW2,'k');
    legend('SSM','RZSM','location','best');
    set(gca,'XAxisLocation','top')
    datetick('x',12,'keeplimits')
    ylabel({'saturation','[-]'});
    grid on, box on,
    s=(['NSE= ',num2str(NS,'%3.3f'),...
        ' ANSE= ',num2str(ANSE,'%3.3f'),...
        ' NSE(radQ)= ',num2str(NS_radQ,'%3.3f'),...        
        ' R^2= ',num2str(RQ,'%3.3f'),...
        ' KGE= ',num2str(KGE,'%3.3f')]);
    nname=[name,name_suff];nname(nname=='_')='-';
    title(['\bf',nname,': ',s]); hold on
    axis([D(1) D(end)+1 -0.05 1.05])
    %  observed rainfall
    h(2) = axes('Position',[0.1 0.65 0.8 0.10]);
    plot(D,PIO),hold on,
    datetick('x',12,'keeplimits')
    set(gca,'ydir','reverse')
    set(gca,'Xticklabel','')
    ylabel({'rainfall','[mm]'});
    grid on, box on
    axis([D(1) D(end)+1 0 nanmax(PIO).*1.05])
    % observed and simulated discharge
    h(3) = axes('Position',[0.1 0.1 0.8 0.54]);
    plot(D,Qobs,'-g','Linewidth',3);hold on
    plot(D,Qsim,'--r','Linewidth',1.0);
    legend('Q_o_b_s','Q_s_i_m','location','best');
    ylabel({'discharge',' [m^3/s]'})
    grid on, box on, axis tight
    axis([D(1) D(end)+1 0 nanmax(Qobs).*1.2]);
    datetick('x',11,'keeplimits')
    print(gcf,[dir,'\',name,name_suff],'-dpng','-r250')
end

% -------------------------------------------------------------------------------
% Calculation of Geomorphological Instantaneous Unit Hydrograph
% -------------------------------------------------------------------------------
function IUH=IUH_comp(gamma,Ab,dt,deltaT)

Lag=(gamma*1.19*Ab^0.33)/deltaT;
hp=0.8/Lag;
data=load('IUH.txt');
t=data(:,1)*Lag;IUH_0=data(:,2)*hp;
ti=0:dt:max(t);
IUH=interp1(t,IUH_0,ti)';

% -------------------------------------------------------------------------------
% Calculation of Nash Instantaneous Unit Hydrograph
% -------------------------------------------------------------------------------
function IUH=IUH_NASH(n,gamma,Ab,dt,deltaT)

K=(gamma*1.19*Ab.^.33)./deltaT;
time=0:dt:100;
IUH=((time/K).^(n-1).*exp(-time/K)/factorial(n-1)/K)';

%--------------------------------------------------------------------------
% Snow accumulation-melting MODEL
%--------------------------------------------------------------------------

function [rainfall,SWE_melting,SWE_snowpack]=snow_model(precipitation, temperature, temp_min, temp_max, Cm)

rainfall = zeros(length(precipitation),1);
snowfall = zeros(length(precipitation),1);
SWE_snowpack = zeros(length(precipitation),1);
SWE_melting = zeros(length(precipitation),1);

% The precipitation is divided into rainfall and snowfall
% REFERENCES:
% U.S. Army Corps of Engineers (1956)
% Boscarello, L., Ravazzani, G., Pellegrini, M., Dedieu, J. P., & Mancini, M. (2014). Calibration of hydrological model FEST from MODIS images in Alpine Catchments. Politecnico di Milano, Dipartimento di Ingegneria Idraulica, Ambientale, Infrastrutture viarie, Rilevamento.
% Degree Day Method (Mockus, 1964)

% INITIALIZATION

if precipitation(1,1) == NaN || temperature(1,1) == NaN
    rainfall(1,1) = NaN;
    snowfall(1,1) = NaN;
elseif temperature(1,1) <= temp_min
    snowfall(1,1) = precipitation(1,1);
    rainfall(1,1) = 0;
    SWE_snowpack(1,1) = snowfall(1,1); % [mm]
    SWE_melting(1,1) = 0; % [mm]
elseif temperature(1,1) >= temp_max
    snowfall(1,1) = 0;
    rainfall(1,1) = precipitation(1,1);
    SWE_snowpack(1,1) = 0; % [mm]
    SWE_melting(1,1) = 0; % [mm]
else
    rainfall(1,1) = precipitation(1,1) * ((temperature(1,1)-temp_min)/(temp_max-temp_min));
    snowfall(1,1) = precipitation(1,1) - rainfall(1,1);
    SWE_snowpack(1,1) = snowfall(1,1);
    SWE_melting(1,1) = 0;
end

% rho_fresh_snow(1,1) = 10^3 * (0.05 + ((temperature(1,1) * 1.8 + 32) / 100)^2); % Fresh snow density [kg/m3]

for i=2:length(precipitation)
    if precipitation(i,1) == NaN || temperature(i,1) == NaN
        % se manca il dato, metto i NaN
        rainfall(i,1) = NaN;
        snowfall(i,1) = NaN;
    elseif temperature(i,1) <= temp_min
        % if the temperature is less than the low threshold, 
        % the precipitation is entirely snowfall
        rainfall(i,1) = 0;
        snowfall(i,1) = precipitation(i,1);
        SWE_snowpack(i,1) = SWE_snowpack(i-1,1) + snowfall(i,1);
        SWE_melting(i,1) = 0;
    elseif temperature(i,1) > temp_max
        % if the temperature is more than the high threshold,
        % the precipitation is entirely rainfall
        rainfall(i,1) = precipitation(i,1);
        snowfall(i,1) = 0;
        SWE_melting(i,1) = Cm * (temperature(i,1) - temp_max);
        % h_melting(i,1) = rho_water * SWE_melting(i,1) / rho_snow;
        % Check the snowpack SWE
        if SWE_snowpack(i-1,1) >= SWE_melting(i,1)
            SWE_snowpack(i,1) = SWE_snowpack(i-1,1) - SWE_melting(i,1);
            % h_snowpack(i,1) = h_snowpack(i-1,1) - h_melting(i,1);
        else
            SWE_melting(i,1) = SWE_snowpack(i-1,1);
            % h_melting(i,1) = h_snowpack(i-1,1);
            % h_snowpack(i,1) = 0;
            SWE_snowpack(i,1) = 0;
        end
    else
        rainfall(i,1) = precipitation(i,1) * ((temperature(i,1)-temp_min)/(temp_max-temp_min));
        snowfall(i,1) = precipitation(i,1) - rainfall(i,1);
        SWE_snowpack(i,1) = SWE_snowpack(i-1,1) + snowfall(i,1);
        SWE_melting(i,1) = 0;
    end
end

function [ kge , r, relvar, bias ] = klinggupta(modelled,observed)
%Nash Sutcliffe Efficiency measure

modelled(isnan(observed))=NaN;

cflow=[modelled,observed];
cflow=rem_nan(cflow);

sdmodelled=nanstd(modelled);
sdobserved=nanstd(observed);
 
mmodelled=nanmean(modelled);
mobserved=nanmean(observed);

r=corrcoef(cflow,'rows','pairwise');
r=r(1,2);
relvar=sdmodelled/sdobserved;
bias=mmodelled/mobserved;

%KGE timeseries 
kge=1-  sqrt( ((r-1)^2) + ((relvar-1)^2)  + ((bias-1)^2) );
 
function [X,ID_NaN] = rem_nan(X)
ID_NaN=any(isnan(X)');
X(ID_NaN,:) = [];
% X(any(isnan(X),2),:) = [];
