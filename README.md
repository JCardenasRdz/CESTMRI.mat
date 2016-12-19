# CESTMRI.mat
CEST.mat is MATLAB toolbox for simulating and fitting an arbitrary number of exchange pools in CEST MRI and NMR.

    [zSpectrum,Mstate]= CESTnPools(magField, satTime, satPower, PPM, ParamVec)
### SIMULATION OF CEST MRI DATA
In order to simulate CEST MRI data using the `CESTnPools` function, the following parameters must be defined:

_Single floating point numbers_

        magField    The magnetic field in Tesla
        satTime     The saturation time in seconds
        satPower    The saturation power in micro-Tesla

_Vectors_

        SaturationOffset    1 x K vector of saturation offsets in parts per million (ppm)
        ParamVec            1 x (5n-2) vector, where n is the number of CEST pools including water
        The vector is constructed by concatenating the following five vectors:

ParamVec= `[T1 T2 Conc ExchangeRates Pool_Offsets]`;

    T1              1 x n vector of T1 relaxation times in seconds
    T2              1 x n vector of T2 relaxation times in seconds
    Conc            1 x (n-1) vector of concentrations in moles/Lt, Water is 110 Moles/Lt 
    ExchangeRates   1 x (n-1) vector of exchange rates from each pool to water in Hz
    Pool_Offsets    1 x n vector of pool offsets in ppm (including water)
    Outputs

The `CESTnPools` function will solve the time-dependent Bloch equations using the parameters defined in the program.

    zSpectrum   Z-spectra of length K (same as saturation offsets) 
    Mstate      (3n+1) X K matrix of final magnetization components for each pool with the following format: 
    Mstate=     [Mx,A ... Mx,N My,A ... My,N Mz,A ... Mz,N 1]' x K

## Notation
Water should always be assigned to pool A. All other pools should be assigned to other letters, from pool B, pool C, all the way to pool N. This notation will be used in the following paragraphs when constructing the five vectors. The T1, T2, and pool saturation vectors are constructed in the following form:

T1=     [T1A T1B T1C ... T1N]; 
T2=     [T2A T2B T2C ... T2N];
PPM=    [ppmA ppmB ppmC ... ppmN];
T1A, T2A, and ppmA are all associated with pool A (water); T1B, T2B, and ppmB are all associated with pool B; and so on. The concentrations and exchange rates vectors are constructed in another form:

Conc=           [CB CC ... CN];
ExchangeRates=  [kBA kCA ... kNA];
These are the exchange rates from one pool to water. Notice that the components of these two vectors begin with pool B, not pool A. Like before, CB and kBA are associated with pool B, CC and kCA are associated with pool C, and so on.

FITTING OF CEST MRI DATA
In order to fit CEST MRI data using the CESTnPools function, some known floating-point numbers must be defined. These include the magnetic field in Tesla’s (magField), the saturation time in seconds (satTime), the saturation power in micro-Tesla’s, and the number of pools (nPools).

Secondly, some curve fitting parameters must be defined. These include the curve fitting function (CESTfunc), the x-data (xData), the y-data (yData), initial conditions (InitCond), and lower (LowerBounds) and upper bounds (UpperBounds).

The curve fitting function is defined by the following, where more details on the variables can be found in the description of simulating CEST MRI data.

CESTfunc= @(ParamVec, SaturationOffset) CESTnPools(magField, satTime, satPower, SaturationOffset, ParamVec);
The x- and y-data come from imported data. The x-data is the saturation offset vector in parts-per-million and the y-data is the associated Z-spectrum in arbitrary units.

The initial conditions and lower and upper bounds for the CEST parameters follow the same format as the parameter vector (ParamVec) in the CEST simulation section. Initial conditions is the parameter vector for expected results, the lower bounds vector is the parameter vector of lower bounds, and the upper bounds vector is the parameter vector of upper bounds. The more precisely these three vectors are defined, the more accurate the fitting will be.

After all parameters are defined, curve fitting can commence. The lsqcurvefit function is used with the following syntax, where resnorm, residual, and jacobian are defined by MATLAB’s help features:

[ParamVecPred,resnorm,residual,~,~,~,jacobian]= lsqcurvefit(CESTfunc, InitCond, xData, yData, LowerBounds, UpperBounds);
ParamVecPred contains the fitting results in the format of ParamVec. The confidence intervals (conf) are determined using the nlparci function with the following syntax:

conf= nlparci(ParamVecPred,residual,'jacobian',jacobian);
The interpretation of conf is located within MATLAB’s help features. In order to plot the results of the CEST fitting, CESTfunc is used again to determine the predicted fit (yPred) with the following syntax:

yPred= CESTfunc(ParamVecPred, xData);
Plotting yPred versus xData using the plot function will result in the desired graph for the fitted curve.
