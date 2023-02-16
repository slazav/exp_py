## fit_res -- Python library for fitting resonance data

Usage: `import fit_res001 as fit_res`

### Functions (v001)

#### `fit` -- fit frequensy sweep

Usage: `res = fit(data, coord=0, npars=6, do_fit=1)`

Function is very similar to fit_res command-line program
(`https://github.com/slazav/fit_res`).

Fit data (X,Y vs F) with one of functions:
* `(X + 1j*Y) = (A + 1j*B) + (C + 1j*D)/(f0**2 - F**2 + 1j*F*df) + (F + 1j*E)*(F-f0)` -- coordivnate fit
* `(X + 1j*Y) = (A + 1j*B) + 1j*F* (C + 1j*D)/(f0**2 - F**2 + 1j*F*df) + (F + 1j*E)*(F-f0)` -- velocity fit
* same with `E = F = 0`

Parameters:
* `data` -- Nx5 numpy array with usual columns T-F-X-Y-D. Can be obtained with `f4wire.get_sweep_*` functions.
* `coord` -- boolean value, do coordinate (1) or velosity (0) fit.
* `npars` -- use 6 or 8 parameters (constant or linear offset)
* `do_fit` -- if 0 then only estimation of initial parameter is done, no fitting.

Return value: A member of `fit_res_t` class with folloowing fields:
* `par` -- parameters, list with 8-values: A,B,C,D,f0,df,E,F
* `err` -- parameter uncertainties, list with 8 values
* `time` -- Mean time (from initial data)
* `drive` -- Mean drive (from initial data)
* `e` -- RMS difference between data and model
* `coord` -- coord/velocity fit
* `npars` -- numpber of free parameters, 6 or 8
* `(A,B,C,D,f0,df,E,F)` -- individual fit parameters
* `__init__(time,drive,e,par,err, npars, coord)` -- constructor, used to set all values
* `func(f)` -- Access to model function.
* `dbfmt()` -- Data list in the same form as in fit databases (19 columns: time, drive, e, A, Ae, B, Be, ...)

#### `fitfunc` -- function used for fitting

Usage: `v = fitfunc(par,coord,F)`

