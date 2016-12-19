% CESTNPOOLSEXAMPLE predicts the evolution of a n-pool system under chemical 
% exchange using least-squares regression
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Syntax
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Inputs
%
%       magField, Magnetic Field [Tesla]
%       satTime, Saturation Time [Seconds]
%       satPower, Saturation Power [Micro Tesla]
%       nPools, Number of Pools
%
%       x-Data, Saturation Offsets [ppm]                           
%       y-Data, Z-magnetization (Mz/Mo) [%]
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Authors
% Jaden Johnston                    Julio Cardenas-Rodriguez
% University of Arizona             University of Arizona
% jadenjohnston@email.arizona.edu   cardenaj@email.arizona.edu
%                       
% www.cardenaslab.org/resources
% v1.0 05/27/2016

% Code
%% Allocate variables
magField    = 7;    % Tesla
satTime     = 2;    % Seconds
satPower    = 1.4;  % MicroTeslas
nPools= 2;

SaturationOffset= linspace(-40,40,200);          % ppm

T1= [2.5 .62];                           % Seconds
T2= [.7 22E-2];                          % Seconds
Conc= [.1];                              % Moles/Lt
ExchangeRates= [500];                     % Hz
PoolOffsets= [0 -10];                     % ppm

ParamVec= [T1 T2 Conc ExchangeRates PoolOffsets];

if(size(ParamVec)~= 5*nPools-2)
    error('Please check size of ParamVec. Conc and ExchangeRates arrays do not contain values for water.');
    exit;
end

[Z,M]= CESTnPools(magField, satTime, satPower, SaturationOffset, ParamVec);
plot(SaturationOffset,100.*Z','LineWidth',2);
axis([-40 40 -1 101]); 
xlabel('Saturation Offset (ppm)');
ylabel('Mz/Mo (%)');

clearvars -except Z M SaturationOffset

%% Define Variables
magField    = 7;    % [Tesla]
satTime     = 2;    % [Seconds]
satPower    = 1.4;    % [Micro Tesla]

nPools      = 2;    % Number of pools

%% Curve Fitting Parameters
% Function
CESTfunc=@(ParamVec, SaturationOffset) CESTnPools(magField, satTime, satPower, SaturationOffset, ParamVec);

% x-Data [ppm]      
xData= SaturationOffset;                     

% y-Data (Mz/Mo) [%]
yData= Z; 

% Guess Initial Conditions
%     InitCond= 1XK vector: [T1 T2 Conc Rate PPM] where 
%         T1=   [ T1a  T1b  T1c ...  T1n]  Relaxation times [s]
%         T2=   [ T2a  T2b  T2c ...  T2n]  Relaxation times [s]
%         Conc= [       Cb   Cc ...   Cn]  Concentration of pool proton [M]
%         Rate= [      kba  kca ...  kna]  Exchange rates to water [Hz]
%         PPM=  [PPMa PPMb PPMc ... PPMn]  Parts per million of pools [ppm]
% 
% ## Water is always the first pool and the elements of Rate and Conc
% should start with the second pool. ##

T1= repmat(2.5,1,nPools);                          % [Seconds] 
T2= repmat(.7,1,nPools);                % [Seconds] 
Conc= repmat(.1,1,nPools-1);       % [Moles/Liter] 
Rate= repmat(500,1,nPools-1);                     % [Hertz] 
PPM= [0 -10];                      % [ppm]
InitCond= [T1 T2 Conc Rate PPM];   
clearvars T1 T2 Rate PPM;

% Lower bounds
T1LB= repmat(.1,1,nPools);
T2LB= repmat(.2,1,nPools);
ConcLB= zeros(1,nPools-1);
RateLB= zeros([1,nPools-1]);
PPMLB= repmat(min(xData),1,nPools);
LowerBounds= [T1LB T2LB ConcLB RateLB PPMLB];
clearvars T1LB T2LB ConcLB RateLB PPMLB;

% Upper bounds
T1UB= repmat(5,1,nPools);
T2UB= repmat(1,1,nPools);
ConcUB= repmat(110,1,nPools-1);
RateUB= Inf(1,nPools-1);
PPMUB= repmat(max(xData),1,nPools);
UpperBounds= [T1UB T2UB ConcUB RateUB PPMUB];
clearvars T1UB T2UB ConcUB RateUB PPMUB Conc;

if(size(InitCond)~= 5*nPools-2)
    error('Please check size of InitCond. Conc and Rate arrays do not contain values for water.');
    exit;
elseif(size(LowerBounds)~= 5*nPools-2)
    error('Please check size of LowerBounds. ConcLB and RateLB arrays do not contain values for water.');
    exit;
elseif(size(UpperBounds)~= 5*nPools-2)
    error('Please check size of UpperBounds. ConcUB and RateUB arrays do not contain values for water.');
    exit;
end

%% Curve fitting with confidence interval
[ParamVecPred,resnorm,residual,~,~,~,jacobian]=lsqcurvefit(CESTfunc, InitCond, xData, yData, LowerBounds, UpperBounds);
conf= nlparci(ParamVecPred,residual,'jacobian',jacobian);

%% Plot Results
yPred= CESTfunc(ParamVecPred, xData);

h=plot(xData,[yPred,yData]); 
h(1).LineWidth=2;   h(1).LineStyle='--'; 
h(2).LineStyle='none';  h(2).Marker='o';    h(2).MarkerSize=2;   
legend({'Predicted','Observed'},'FontSize',14,'Location','SouthWest');
xlabel('Offset (ppm)');
ylabel('Mz/Mzo (au)');