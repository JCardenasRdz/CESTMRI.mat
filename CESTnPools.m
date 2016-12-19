function [zSpectrum,Mstate] = CESTnPools(magField, satTime, satPower, PPM, ParamVec)
% CESTNPOOLS calculates the evolution of a n-pool system under chemical 
% exchange using a matrix exponential approximation to solve the equation,
% dM/dt = A*M
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Syntax
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% [zSpectrum] = CESTnPools(magField, satTime, satPower, PPM, ParamVec)
%
% Inputs
%
%     magField= 1X1 Scalar: magnetic field strength [Tesla]
%     satTime=  1X1 Scalar: Saturation time [seconds]
%     satPower= 1X1 Scalar: RF Saturation power [micro Tesla]
%     PPM=      KX1 vector: Saturation Offsets [ppm]
%     ParamVec= 1XK vector: [T1 T2 Conc Rate PPM] where 
%         T1=   [ T1a  T1b  T1c ...  T1n]  Relaxation times [s]
%         T2=   [ T1a  T2b  T2c ...  T2n]  Relaxation times [s]
%         Conc= [       Cb   Cc ...   Cn]  Concentration of pool proton [M]
%         Rate= [      kba  kca ...  kna]  Exchange rates to water [Hz]
%         PPM=  [PPMa PPMb PPMc ... PPMn]  Parts per million of pools [ppm]
% 
% ## Water is always the first pool and the elements of Rate and Conc
% should start with the second pool. ##
%
% Outputs
%
%     zSpectrum= KX1 vector of Mz component after water
%     Mstate=    [Mxa ... Mxp Mya ... Myp Mz0a ... Mz0p 1]' X K Matrix
%                of final state for each Saturation Offset;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% EXAMPLE
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% magField    = 7;    % Tesla
% satTime     = 2;    % Seconds
% satPower    = 1.4;  % MicroTeslas
% 
% SaturationOffset= linspace(-40,40,200);          % ppm
% 
% T1= [2.5 .62 .56 .770];                           % Seconds
% T2= [.7 22E-2 33E-3 33E-3];                       % Seconds
% Conc= [.1 .1 .050];                               % Moles/Lt
% 
% ExchangeRates= [500 500 500];                     % Hz
% PoolOffsets= [0 -10 20 -20];                     % ppm
% ParamVec= [T1 T2 Conc ExchangeRates PoolOffsets];
% 
% Z= CESTnPools(magField, satTime, satPower, SaturationOffset, ParamVec);
% plot(SaturationOffset,100.*Z','LineWidth',2);
% axis([-40 40 -1 101]); 
% xlabel('Saturation Offset (ppm)');
% ylabel('Mz/Mo (%)');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Reference
% Murase K, and Tanki Nobuyoshi. Magnetic Resonance Imaging 29 (2011) 126; 131
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Authors
% Jaden Johnston                    Julio Cardenas-Rodriguez
% University of Arizona             University of Arizona
% jadenjohnston@email.arizona.edu   cardenaj@email.arizona.edu
%                       
% www.cardenaslab.org/resources
% v1.0 05/27/2016


% Code
%% Properties of Hydrogen
gyroRatH= 42.576;            % Gyromagnetic Ratio [MHz/T]
concH= 110;                  % Molarity [M]

%% System Parameters
[~,cols]= size(ParamVec);
    nPools= (cols+2)/5;    % Number of pools
    sizeA= 3*nPools+ 1;    % Size of A (square matrix)

R1= 1./ ParamVec(1:nPools)';                                    % T1 Relaxation Rates [Hz]
R2= 1./ ParamVec(nPools+1:2*nPools)';                           % T2 Relaxation Rates [Hz]
fConc= zeros([nPools, 1]);                                      % Fractional concentration of the solute protons
    fConc(2:end)= ParamVec(2*nPools+1:3*nPools-1)./ concH;      
    fConc(1)= 1;                                                    % Water has ratio of 1
RateNA= zeros([nPools, 1]);                                     % Exchange rates from pools to water [Hz]
    RateNA(2:end)= ParamVec(3*nPools:4*nPools-2);                  
RateAN= RateNA.* fConc;                                             % Exchange rates from water to pools [Hz]
    RateNA(1)= sum(RateAN(2:end));                                  % Sum of exchange rates from water to pools [Hz]
PowerOffset= PPM.* gyroRatH.* magField.* 2.* pi;                 % Offset frequency of RF irradiation [rad/s]
    W= ParamVec(4*nPools-1:end)'.* gyroRatH.* magField.* 2.* pi;    % Larmor frequencies [rad/s]
    
w1= satPower* gyroRatH* 2* pi;  % Nutation rate of RF irradiation [rad/s]

% Initial Conditions
M0= zeros([sizeA 1]);
M0(2*nPools+1:3*nPools)= fConc;
M0(end)=1;

%% Allocate Propagation Matrix
%{ 
    Matrix A is in the form (see Fig 4 in Murase):
         A1     A5     0     0
        -A5     A1    A4     0
          0    -A4    A2    A3
          0      0     0     0
%}

CA= cell(4);

A1= -diag(R2+ RateNA);
    A1(1,2:end)= RateNA(2:end);
    A1(2:end,1)= RateAN(2:end);
A2= -diag(R1+ RateNA);
    A2(1,2:end)= RateNA(2:end);
    A2(2:end,1)= RateAN(2:end);
A3= R1.* fConc;
A4= w1*eye(nPools);

CA{1,1}= A1;
CA{2,2}= A1;
CA{3,3}= A2;
CA{3,4}= A3;
CA{2,3}= A4;
CA{3,2}= -A4;

CA{1,3}= zeros(nPools);
CA{1,4}= zeros(nPools,1);
CA{2,4}= zeros(nPools,1);
CA{3,1}= zeros(nPools);
CA{4,1}= zeros(1,nPools);
CA{4,2}= zeros(1,nPools);
CA{4,3}= zeros(1,nPools);
CA{4,4}= 0;

zSpectrum=zeros(length(PPM),1);
Mstate= zeros(nPools*3+1,length(PPM));

for i=1:length(PowerOffset)
    A5= diag(W- PowerOffset(i));
    
    CA{1,2}= A5;
    CA{2,1}= -A5;
    
    % Calculate M for each offset frequency
    A= cell2mat(CA);
    M= expm(A*satTime)*M0;
    Mstate(:,i)=M;
    zSpectrum(i,1)= M(2*nPools+ 1);
end    

end

