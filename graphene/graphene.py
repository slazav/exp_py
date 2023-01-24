import subprocess
import warnings
import numpy
import re
import os.path

sources = {
  "local":   ('device_c', 'ask', 'db'),
  "ssh_f4b": ('ssh', 'f4b', 'device_c', 'ask', 'db'),
  "xyz_f4":  ('wget', 'http://slazav.xyz:8095/', '-O', '-', '-o', '/dev/null'),
  "xyz_f2":  ('wget', 'http://slazav.xyz:8091/', '-O', '-', '-o', '/dev/null'),
  "xyz_dd":  ('wget', 'http://slazav.xyz:8085/', '-O', '-', '-o', '/dev/null'),
}

gr_args = ['device_c', 'ask', 'db']
cache_dir = 'data'

### Set program for accessing graphene
### - list of arguments or name from sources table
def set_source(a):
  global gr_args
  if type(a) == list or type(a) == tuple:
    gr_args = list(a)
  else:
    gr_args = sources[a]

### Set cache folder (empty for no caching)
def set_cache(a):
  global cache_dir
  cache_dir = a


### Do arbitrary graphene command
def graphene_cmd(cmd, name, t1=0, t2='inf', unpack=False, usecols=None, fname=""):

  ### do we want to use cache?
  if cache_dir != "":
    if not os.path.isdir(cache_dir): os.mkdir(cache_dir)
    if fname=="":
      fname = name + "_" + cmd + "_" + str(t1) + "_" + str(t2)
    fname = cache_dir + "/" + fname
    if os.path.isfile(fname):
      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = numpy.loadtxt(fname, comments="#", usecols=usecols, unpack=unpack)
      return data

  # convert timestammps to strings if needed
  if type(t1) != str: t1='%f'%(t1)
  if type(t2) != str: t2='%f'%(t2)

  ### build command: gr_args + <command> + <db name> + <times>
  args = list(gr_args)
  if len(gr_args)>1 and gr_args[0] == 'wget':
    if args[1][-1] != '/': args[1] += '/'
    args[1] += "%s?name=%s"%(cmd,name)
    if cmd=='get_wrange': args[1]+="&t1=%s&t2=%s"%(t1,t2)
    if cmd=='get_range':  args[1]+="&t1=%s&t2=%s"%(t1,t2)
    if cmd=='get_next':   args[1]+="&t1=%s"%(t1)
    if cmd=='get_prev':   args[1]+="&t2=%s"%(t2)
    if cmd=='get':        args[1]+="&t2=%s"%(t2)
  else:
    args.extend((cmd, name))
    if cmd=='get_wrange': args.extend((t1,t2))
    if cmd=='get_range': args.extend((t1,t2))
    if cmd=='get_next':  args.append(t1)
    if cmd=='get_prev':  args.append(t2)
    if cmd=='get':       args.append(t2)

  ### run the command and load values
  with subprocess.Popen(args,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True) as proc:
    if cache_dir != "":
      print(proc.stdout.read(), file=open(fname, "w"))
    else:
      fname = proc.stdout
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      data = numpy.loadtxt(fname, comments="#", usecols=usecols, unpack=unpack)
    rc = proc.wait()
    if rc:
      print('> Graphene error:', proc.stderr.read())
      exit(1)
  return data


###

def get_range(name, t1, t2, **kwargs):
  return graphene_cmd('get_range', name, t1=t1, t2=t2, **kwargs)

def get_wrange(name, t1, t2, **kwargs):
  return graphene_cmd('get_wrange', name, t1=t1, t2=t2, **kwargs)

def get_prev(name, t, **kwargs):
  return graphene_cmd('get_prev', name, t2=t, **kwargs)

def get_next(name, t, **kwargs):
  return graphene_cmd('get_next', name, t1=t, **kwargs),

def get(name, t, **kwargs):
  return graphene_cmd('get', name, t2=t, **kwargs)

###

# Convert time from human-readable to graphene-readable form
def timeconv(t):
  if t=="now" or t=="now_s" or t=="inf": return t
  if re.fullmatch('[0-9]+[+-]?', t): return t
  return subprocess.check_output(['date', '+%s', '-d', t]).decode('utf-8')[0:-1]
