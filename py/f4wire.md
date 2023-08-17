## Python library for getting vibrating wire data from graphene database

Usage: `import f4wire001 as f4wire`

### Changelog

* v001 -- 2023.02.14, V.Z

### Functions (v001)

#### `get_data` -- get data without background, converted to real voltage/current

Usage: `get_data(name, t1, t2, use_bg=1, cnv_drive=1, cnv_volt=1, cache="")`

Works for both tracking and frequency sweeping modes.

Arguments:

* `name` -- wire name, such as 'w1bt'
* `t1`, `t2` -- timestamps, same as in `graphene.py` library
* `use_bg` -- get background from `<name>_pars:f2` database, subtract it from data
* `cnv_drive` -- convert drive from Vpp to Arms (using drive box settings)
* `cnv_volt` -- convert voltage to real voltage on the wire (using cold transformer gain)
* `cache` -- if cache file exists load data from there, if not get data from database and save to cache

Return value: empty or 2D numpy array (Nx5). Columns: time, frequency(Hz), X(Vrms), Y(Vrms), D(Vpp or Arms).

#### `get_sweep_*` -- get list of frequency sweeps

* `get_sweep_prev(name, t1, sweep_dir=None, nskip=0, nsweeps=1, ...)` -- get sweep starting before t1
* `get_sweep_next(name, t1, sweep_dir=None, nskip=0, nsweeps=1, ...)` -- get sweep starting after t1
* `get_sweep_range(name, t1, t2, sweep_dir=None, ...)` -- get all sweeps starting between t1 and t2
* `get_sweep(name, t1, sweep_dir=None, ...)` -- get sweep starting before t1 and ending after t1
* `get_sweep_list(name, tlist, sweep_dir=None, ...)` -- similar to `get_sweep()` but use list of timestamps to get list of sweeps

Arguments:

* `sweep_dir` -- is -1 or 1, then select only one sweep direction.
* `nskip`    -- number of sweeps to skip in get_sweep_prev/get_sweep_next.
* `nsweeps`  -- number of sweeps to get in get_sweep_prev/get_sweep_next.
* other arguments are same as in `get_data()` function
* if `cache` parameter is not empty then all sweeps are saved in a single file

Return value: Python list of sweeps - numpy arrays from get_data() function

#### `merge_sweeps` -- merge sweeps

Only merge sweeps with same drive (same_drive=1) or merge all sweeps with amplitude restcaling.

Usage: `sweeps = merge_sweeps(sweeps, same_drive=1)`

#### `track_res_lin` -- Process tracking mode data, linear resonance

Usage: `(f0, df) = track_res_lin(data, fit)`

Process tracking mode data. Assuming that resonance is linear, and
A,B,C,D,E,F parameters are proportional to drive find f0 and df
parameters.

Arguments:
* `data` -- data from `<name>_sweep` database (output of `f4wire.get_data()`)
* `fit`  -- output of `fit_res`, fit of a frequency sweep which is used for processing

Return values: resonance frequency and width, `f0` and `df`, numpy arrays, same length as `data[:,0]`

#### `track_heat` -- Process tracking mode data, heating

Usage: `pwr = track_heat(data, fit)`

Process tracking mode data. Calculate power dissipated by wire. Power in Watts if
input voltage in Volts(rms) and drive is Amps(rms)

Arguments: same as for `track_res_lin()`

Return value: dissipated power `pwr` -- numpy array, same length as `data[:,0]`

#### `get_track` -- high-level function for processing tracking mode data.

Usage: `(time, f0, df, Vpar, Vperp, vel, PWR) =
 get_track(name, t1, t2, get=0, cache="", plot="", nsweeps=1, nskip=0, prev_sweeps=1, fit_coord=0, fit_npars=6, use_bphase=0, bphase=None, verb=0, osc="")`

Arguments:
* `name` -- wire name
* `t1`, `t2` -- timestamps
* `cache` -- prefix for cache file (will be saved in .npz format)
* `get`   -- always get data from database even if cache file exists (it will be updated)
* `plot`  -- prefix for plot (.png)
* `nsweeps`, `nskip`, `prev_sweeps` -- Use sweeps before (`prev_sweeps=1`) or after the data range.
   skip `nskip` sweeps and merge next `nsweeps` (with amplitude rescaled to drive).
* `fit_coord` -- use coordinate fit instead of velocity fit
* `fit_npars` -- number of parameters for Lorenzian fit
* `use_bphase` -- Flag for B-phase corrections. Correction parameters are
   taken from wire_info_t class for the given wire name.
   Corrections are used for fitting initial sweep and for calculating corrected resonance width `dF0`.
* `bphase` -- explicetely set bphase correction parameters (if there is a need to use non-standard ones)
* `verb` - verbosity flag
* `osc` -- scope binary file

Return values:
* `TT` -- Unix time [s]
* `F0`, `dF` -- resonance frequency and width [Hz]
* `dF` -- uncorrected resonance width [Hz]
* `dF0` -- corrected resonance width [Hz]
* `FF` -- array of sweep frequencies [Hz]
* `DD` -- drive current [A]
* `sweep` -- sweep points
* `fit` -- sweep fit results
* `VX`, `VY` -- components of voltage projected to direction of driving force [V].
* `vel` -- velocity amplitude [m/s]. Field from `demag_pc:f2` database is used.
* `PWR` -- dissipated power (drive * vpar) [W].
* `ttc` -- T/T_c
* `wire` -- name of the wire
* `press` -- pressure [bar]
* `field` -- magnetic field [T]

def get_track(name, t1, t2,
     cache="", plot="", nsweeps=1, nskip=0, prev_sweeps=1, fit_coord=0, fit_npars=6):
