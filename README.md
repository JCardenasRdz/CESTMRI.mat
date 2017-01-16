*CEST.mat* is MATLAB toolbox for simulating and fitting an arbitrary number of exchange pools in CEST MRI and NMR.
```Matlab
    [zSpectrum,Mstate]= CESTnPools(magField, satTime, satPower, PPM, ParamVec)
```

## Notation
Water should always be assigned to pool A. All other pools should be assigned to other letters, from pool B, pool C, all the way to pool N. This notation will be used in the following paragraphs when constructing the five vectors needed to simulate data. The `T1`, `T2`, and pool saturation vectors are constructed in the following form:
```Matlab
    T1=     [T1A T1B T1C ... T1N]; 
    T2=     [T2A T2B T2C ... T2N];
    PPM=    [ppmA ppmB ppmC ... ppmN];
```
`T1A`, `T2A`, and `ppmA` are all associated with pool A (water); `T1B`, `T2B`, and `ppmB` are all associated with pool B; and so on. The concentrations and exchange rates vectors are constructed in another form:
```Matlab
    Conc=           [CB CC ... CN];
    ExchangeRates=  [kBA kCA ... kNA];
```
These are the exchange rates from one pool to water. Notice that the components of these two vectors begin with pool B, not pool A. Like before, `CB` and `kBA` are associated with pool B, `CC` and `kCA` are associated with pool C, and so on.

# SIMULATION OF CEST MRI DATA
In order to simulate CEST MRI data using the `CESTnPools` function, the following parameters must be defined:
_Single floating point numbers_
```
magField    The magnetic field in Tesla
satTime     The saturation time in seconds
satPower    The saturation power in micro-Tesla
```
_Vectors_
```
SaturationOffset    1 x K vector of saturation offsets in parts per million (ppm)
ParamVec            1 x (5n-2) vector, where n is the number of CEST pools including water
```
The vector is constructed by concatenating the following five vectors:
```
ParamVec= `[T1 T2 Conc ExchangeRates Pool_Offsets]`;
```
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


# FITTING OF CEST MRI DATA
In order to fit CEST MRI data using the `CESTnPools` function, some of the variables used in the simulation but me fixed using an anonymous function, otherwise the fitting algorithm will try to estiamte them. The two variables thats are always knows are:

       magField = magnetic field in Teslas
       satTime = saturatin time in seconds

If you are very confident on the saturation power applied to the sample you can also fix `satPower`. If you are want to estimate the actual power, just include it as one of the variables. Once you define the anonymous function, you can use it and `lsqcurvefit` (or other minimizer) to estimate the paramters.

_Define anonymous function_
```
CESTfunc= @(ParamVec, SaturationOffset) CESTnPools(magField, satTime, satPower, SaturationOffset, ParamVec);
```

## Perform Non-Linear Leastsquare curve fitting_

    [ParamVecPred,resnorm,residual,~,~,~,jacobian]= lsqcurvefit(CESTfunc, InitCond, SaturationOffset, yData, LowerBounds, UpperBounds);
where:

    InitCond =              Initial guess for all the parameters to be estimated
    SaturationOffset=       Experimental saturation offsets in ppm
    yData=                  Experimental Z spectra
    LowerBounds=            Lower limit for all the parameters to be estimated
    UpperBounds=            Upper limit for all the parameters to be estimated

`ParamVecPred` contains the fitting results in the format of ParamVec. The confidence intervals (conf) can be determined using the `nlparci` function with the following syntax:
```
conf= nlparci(ParamVecPred,residual,'jacobian',jacobian);
```
