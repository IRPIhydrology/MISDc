%% MISDc : "MODELLO IDROLOGICO SEMI-DISTRIBUITO IN CONTINUO"
%                   2-LAYER, LUMPED, SNOW VERSION
% -------------------------------------------------------------------
%           TUTORIAL for the use of the Matlab MISDc CODE
% -------------------------------------------------------------------
% 
% The MISDc model can be freely applied and used, just cite some of the references to the model reported below.
% The authors are highly interested to collaborate for the understanding of the model functioning and to improve its performance and applicability.
% 
% For any questions please do not hesitate to contact:
% luca.brocca@irpi.cnr.it
% 
% -------------------------------------------------------------------
% 
% The following files are distributed
% 
% MATLAB codes
% 1) "MISDc_WEB_2L_snow_IE.m": MISDc code at daily time scale
% 2) "MISDc_WEB_2L_snow_IE_h.m": MISDc code at hourly time scale
% 3) "cal_MISDc_WEB_2L_snow_IE.m": code for MISDc calibration (daily)
% 4) "cal_MISDc_WEB_2L_snow_IE_h.m": code for MISDc calibration (hourly)
%  
% note that the calibration codes require the optimization toolbox
%
% Auxiliary file:
% "IUH": contains the coordinates of the dimensionless IUH 
% (Instantaneous Unit Hydrograph),
% 
% INPUT file (example)
% "migi_0205_daily.txt": example file for Niccone at Migianella (137 km^2)
%      It contains 4 columns
%      a) date (in numeric Matlab format)
%      b) rainfall depth (in mm)
%      c) air temperature (°C)
%      d) observed discharge data (m^3/s)
% 
% OUTPUT file:
% "Qsim_MISDc_2L_daily_test.png": figure with the model output 
% (name of the figure defined in the function)
%      
% Upper panel: relative soil moisture for the two layers and rainfall time series
% Middle panel: rainfall
% Lower panel: Observed (green) vs simulated (red) river discharge
% 
%% Calibration of MISDc 2 layer daily 
clc,clear,close all
% [NS,ANSE,KGE,NSradQ,Xopt,WW,Qsim]=...
%     cal_MISDc_WEB_2L_snow_IE(data,X_ini,Ab,dir,name,name_suff,FIG);
%
% INPUT
% data: Matrix with input data (see above for the structure)
% X_ini: Initial set pf parameter values
% area: Basin area
% dirout: Directory where to save the figure output
% name: Fist part of the figure name
% name_suff: Second part of the figure name
% FIG: 1 for making the figure, otherwise no figure

data=load('Test_data\migi_0205_daily.txt'); % input data set
X_ini=ones(10,1)*.1;X_ini(1)=0.05; % (dimensionless)
area=137; % basin area in km^2
dirout=[cd,'\Results\']; % Results folder
name=['Qsim_MISDc_2L_daily'];
name_suff=['_test'];
FIG=1;

% Calibration function (default calibration wrt KGE, can be changed)
[NS,ANSE,KGE,NSradQ,Xopt,WW,Qsim]=...
    cal_MISDc_WEB_2L_snow_IE(data,X_ini,area,dirout,name,name_suff,FIG);

% OUTPUT 
% NS: Nash Sutcliffe Efficiency 
% ANSE: Nash Sutcliffe Efficiency for high flows
% KGE: Kling Gupta Efficiency
% NSradQ: Nash Sutcliffe Efficiency for low flows
% WW: Simulated soil moisture for the first layer
% Qsim: Simulated river discharge

%% Model run with optimized parameter set (Xopt)
[NS,ANSE,KGE,NS_radQ,Qsim,WW,WW2]=...
    MISDc_WEB_2L_snow_IE(data,Xopt,area,FIG,dirout,name,name_suff);
% Results are the same as in the calibration (of course)

%% A similar routine can be written at hourly time scale

% Calibration of MISDc 2 layer hourly 
% clc,clear,close all
% % [NS,ANSE,KGE,NSradQ,Xopt,WW,Qsim]=...
% %     cal_MISDc_WEB_2L_snow_IE_h(data,X_ini,Ab,dir,name,name_suff,FIG);
% %
% % INPUT
% % data: Matrix with input data (see above for the structure)
% % X_ini: Initial set pf parameter values
% % area: Basin area
% % dirout: Directory where to save the figure output
% % name: Fist part of the figure name
% % name_suff: Second part of the figure name
% % FIG: 1 for making the figure, otherwise no figure
% 
% data=load('Test_data\migi_0205.txt'); % input data set
% X_ini=ones(10,1)*.1;X_ini(1)=0.05; % (dimensionless)
% area=137; % basin area in km^2
% dirout=[cd,'\Results\']; % Results folder
% name=['Qsim_MISDc_2L_hourly'];
% name_suff=['_test'];
% FIG=1;
% 
% % Calibration function (default calibration wrt KGE, can be changed)
% [NS,ANSE,KGE,NSradQ,Xopt,WW,Qsim]=...
%     cal_MISDc_WEB_2L_snow_IE_h(data,X_ini,area,dirout,name,name_suff,FIG);
% 
% % OUTPUT 
% % NS: Nash Sutcliffe Efficiency 
% % ANSE: Nash Sutcliffe Efficiency for high flows
% % KGE: Kling Gupta Efficiency
% % NSradQ: Nash Sutcliffe Efficiency for low flows
% % WW: Simulated soil moisture for the first layer
% % Qsim: Simulated river discharge
% 
% %% Model run with optimized parameter set (Xopt)
% [NS,ANSE,KGE,NS_radQ,Qsim,WW,WW2]=...
%     MISDc_WEB_2L_snow_IE_h(data,Xopt,area,FIG,dirout,name,name_suff);
