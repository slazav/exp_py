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

* `get_sweep_prev(name, t1, sweep_dir=None, ...)` -- get sweep starting before t1
* `get_sweep_next(name, t1, sweep_dir=None, ...)` -- get sweep starting after t1
* `get_sweep_range(name, t1, t2, sweep_dir=None, ...)` -- get all sweeps starting between t1 and t2
* `get_sweep(name, t1, sweep_dir=None, ...)` -- get sweep starting before t1 and ending after t1
* `get_sweep_list(name, tlist, sweep_dir=None, ...)` -- similar to `get_sweep()` but use list of timestamps to get list of sweeps

Arguments:

* `sweep_dir` -- is -1 or 1, then select only one sweep direction.
* other arguments are same as in `get_data()` function
* if `cache` parameter is not empty then all sweeps are saved in a single file

Return value: Python list of sweeps - numpy arrays from get_data() function

#### `merge_sweeps` -- merge sweeps with same drive

Usage: `sweeps = merge_sweeps(sweeps)`
