%-----------------------------------------------------------------------------------------------
%   MILc rainfall-runoff model
%
%   Brocca, L., Liersch, S., Melone, F., Moramarco, T., Volk, M. (2013).
%   Application of a model-based rainfall-runoff database as efficient tool for flood risk management.
%   Hydrology and Earth System Sciences Discussion, 10, 2089-2115.
%------------------------------------------------------------------------------------------------

function [NS,NS_lnQ,NS_radQ,RQ,ANSE]=MILc_WEB(name,PAR,FIG)

% Loading input data
DPEQ=load([name,'.txt']);
[M,N]=size(DPEQ);
D=DPEQ(:,1); PIO=DPEQ(:,2); TEMPER=DPEQ(:,3); Qobs=DPEQ(:,4);
delta_T=round(nanmean(diff(D))*24*10000)/10000;
MESE=month(D);

% Model parameter
W_p       = PAR(1); % initial conditions, fraction of W_max (0-1)
W_max     = PAR(2); % Field capacity
m2        = PAR(3); % exponent of drainage
Ks        = PAR(4); % Ks parameter of infiltration and drainage
Nu        = PAR(5); % fraction of drainage verusu interflow
gamma1    = PAR(6); % coefficient lag-time relationship
Kc        = PAR(7); % parameter of potential evapotranspiration
lambda    = PAR(8); % initial abstraction coefficient
Sr_coeff  = PAR(9); % multiplicative coefficient for Sr

% Other data
FPAR=load('fixed_par.txt');
Ab      = FPAR(1);     % Basin area
dt      = FPAR(2);     % computation time step in hour
Ks      = Ks.*delta_T; % mm/h --> mm/delta_T

% Potential Evapotranspiration parameter
L=[0.2100;0.2200;0.2300;0.2800;0.3000;0.3100;
    0.3000;0.2900;0.2700;0.2500;0.2200;0.2000];
Ka=1.26;
EPOT=(TEMPER>0).*(Kc*(Ka*L(MESE).*(0.46*TEMPER+8)-2))/(24/delta_T);
% EPOT=Kc.*EPOT;
clear DPEQ TEMPER

% Initialization
BF=zeros(M,1);
QS=zeros(M,1);
WW=zeros(M,1);
PERC=zeros(M,1);

% Main ROUTINE
W=W_p*W_max;
PIOprec=0;
S=NaN;
Pcum=0;
IE=0;
for t=1:M
    if  (PIO(t)>0.0) && (PIOprec==0)
        S=Sr_coeff*(W_max-W);
        Ia=lambda*S;
    end
    if (PIO(t)==0) && (PIOprec==0) && not(isnan(S))
        S=NaN;
        Pcum=0;
    end
    
    if not(isnan(S))
        if Pcum<Ia
            IE=0;
            Pcum=Pcum+PIO(t);
        end
        if Pcum>=Ia
            Pcum1=Pcum+PIO(t);
            IE=((((Pcum1-Ia).^2)./(Pcum1+(1-lambda)*S))...
                -(((Pcum-Ia).^2)./(Pcum+(1-lambda)*S)));
            Pcum=Pcum1;
        end
    end
    
    E=EPOT(t)*W/W_max;
    PERC(t)=Nu*Ks*(W/W_max).^(m2);
    BF(t)=(1-Nu)*Ks*(W/W_max).^(m2);
    W=W+(PIO(t)-BF(t)-IE-PERC(t)-E);
    if W>=W_max
        SE=W-W_max;
        W=W_max;
    else
        SE=0;
    end
    QS(t)=IE+SE;
    WW(t)=W./W_max;
    if t>3,PIOprec=sum(PIO(t-3:t));end
end

% Convolution (GIUH)
IUH1=IUH_comp(gamma1,Ab,dt,delta_T)*dt;IUH1=IUH1./sum(IUH1);
IUH2=IUH_NASH(1,0.5*gamma1,Ab,dt,delta_T)*dt;IUH2=IUH2./sum(IUH2);
% IUH2=IUH_comp(gamma2,Ab,dt,delta_T);
QSint=interp1(1:M,QS,1:dt:M)';
BFint=interp1(1:M,BF,1:dt:M)';
temp1=conv(IUH1,QSint);
temp2=conv(IUH2,BFint);
Qsim1=temp2(1:round(1/dt):M*round(1/dt)).*(Ab.*1000./delta_T./3600);
Qsim=(temp1(1:round(1/dt):M*round(1/dt))+temp2(1:round(1/dt):M*round(1/dt)))...
    .*(Ab.*1000./delta_T./3600);

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

% Figure
if FIG==1
    % saturation degree
    figure,
    set(gcf,'paperpositionmode','manual','paperposition',[1 1 24 13])
    set(gcf,'position',[50   50   900   500])
    h(1) = axes('Position',[0.1 0.77 0.8 0.10]);
    s=(['NSE= ',num2str(NS),...
        ' ANSE= ',num2str(ANSE),...
        ' NSE(radQ)= ',num2str(NS_radQ),...        
        ' R^2= ',num2str(RQ)]);
    title(['\bf',s]); hold on
    plot(D,WW,'r');
    set(gca,'XAxisLocation','top')
    datetick('x',12,'keeplimits')
    ylabel({'saturation','[-]'});
    grid on, box on,
    axis([D(1) D(end) -0.05 1.05])
    %  observed rainfall
    h(2) = axes('Position',[0.1 0.65 0.8 0.10]);
    plot(D,PIO),hold on,
    datetick('x',12,'keeplimits')
    set(gca,'ydir','reverse')
    set(gca,'Xticklabel','')
    ylabel({'rainfall','[mm]'});
    grid on, box on, axis tight 
    % observed and simulated discharge
    h(3) = axes('Position',[0.1 0.1 0.8 0.54]);
    plot(D,Qobs,'-g','Linewidth',3);hold on
    plot(D,Qsim,'--r','Linewidth',1.0);
    legend('Q_o_b_s','Q_s_i_m',0);
    ylabel({'discharge',' [m^3/s]'})
    grid on, box on, axis tight
    datetick('x',11,'keeplimits')
    print(gcf,name,'-dpng','-r150')
%     save res
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
