## Python library for getting data from graphene database

Usage: `import graphene002 as graphene`

### Changelog

* 001 - first version, until 13.02.2023, V.Z.
* 002 - 13.02.2003 - ...
        more strict cache files (name should be specified)

### Functions (v002)

#### `set_source` -- Configure access to graphene database.

Usage `set_source(a)`

Argument: list of arguments for executing a command, or pre-defined name:

* "local":   ('device_c', 'ask', 'db') -- should be used on a computer with
configured device `db` for `device_c` program (f4a, f4b).

* "ssh_f4b": ('ssh', 'f4b', 'device_c', 'ask', 'db') -- should be used if ssh access
to f4a or f4b is available.

* "xyz_f4":  ('wget', 'http://slazav.xyz:8095/', '-O', '-', '-o', '/dev/null') -- http access
* "xyz_f2":  ('wget', 'http://slazav.xyz:8091/', '-O', '-', '-o', '/dev/null')
* "xyz_dd":  ('wget', 'http://slazav.xyz:8085/', '-O', '-', '-o', '/dev/null')

For very fast access one can copy database file locally and use
`set_source(('graphene', '-d', '.', '-E', 'none'))`.


#### `graphene_cmd`, `get_range`, `get_wrange`, `get_prev`, `get_next`, `get` -- Get data from graphene.

Usage: `graphene_cmd (cmd, name, t1="0", t2='inf', dt=0, raw=False, usecols=None, unpack=False, cache="")`

A separate function is available for each command:
*`get_range(name, t1, t2, **kwargs)`
*`get_wrange(name, t1, t2, **kwargs)`
*`get_prev(name, t, **kwargs)`
*`get_next(name, t, **kwargs)`
*`get(name, t, **kwargs)`

See `https://github.com/slazav/graphene` for command description.

Arguments:

* `cmd` -- Command name: "get_range", "get_wrange", "get_prev", "get_next", or "get".
* `name` -- Database name. Can contain column name "<name>:<n>" or filter "<name>:f<n>"
* `t1`, `t2` -- Timestamps. Same form as for `timeconv` function.
* `dt` -- Time step for `get_range` and `get_wrange` commands.
* `numpy` -- Return python list with original text data instead of numpy array.
* `usecols` -- Use specified columns from the database. Works for both numpy and list output.
   Values larger then number of data columns are allowed (nan values on output).
* `unpack` -- Transpore numpy arrays (same as in `numpy.loadtxt`).
* `cache` -- Filename for using as cache. Data will be saved there and used next time from the file.
Default value "" means no caching.

Return value:

* 2D or empty numpy array (if numpy==True)
* list of lists or empty list (if numpy==False)

#### `graphene_load -- load data from stream

Usage: `graphene_load(ff, unpack=False, usecols=None, numpy=True)`

Function used in `graphene_cmd` for loading data.

#### `graphene_run -- run graphene command and get data

Usage: `graphene_run(cmd, name, t1="0", t2='inf', dt=0)`

Function used in `graphene_cmd` for running graphene program.

#### `timeconv` -- Convert time from human-readable form to unix timestamp.

Usage: `graphene.timeconv(t, fmt="")`

Arguments:

* `t` -- Time in string form. Could be "now" (for current time), "now_s"
(for current time truncated to seconds), "0" or "inf" (for smallest and
largest graphene timestamp), unix timestamp (as string), or time in
human-readable form matching fmt argument.

* `fmt` -- Time format. If empty (default value), a few standard formats
are tried: '%Y-%m-%d %H:%M:%S', '%Y-%m-%dT%H:%M:%S', '%Y-%m-%d %H:%M',
'%Y-%m-%dT%H:%M', '%Y-%m-%d %H', '%Y-%m-%dT%H', '%Y-%m-%d'.

Return value: timestamp suitable for graphene in string form.

Note: graphene uses <seconds>.<nanoseconds> format which is too wide
for standard floating-point numbers. If exact timestamp is needed it's
better to use string for it. In almost all real cases microsecond precision
is used, then it is possible to use float values, converting them to string
with "%.6f" format.

Function used in `graphene_cmd` for converting input timestamps.
