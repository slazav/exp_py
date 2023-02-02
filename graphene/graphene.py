import subprocess
import warnings
import numpy
import re
import tempfile
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

### Convert time from human-readable to graphene-readable form
def timeconv(t):
  if t=="now" or t=="now_s" or t=="inf": return t
  if re.fullmatch('[0-9\.]+[+-]?', t): return t
  return subprocess.check_output(['date', '+%s', '-d', t]).decode('utf-8')[0:-1]

###############################################################
# Communication with graphene, read-only operations.
# No caching, no data parsing.
# t1 and t2 should be strings.
def graphene_read(cmd, name, t1="0", t2='inf', dt=0):

  ### build command: gr_args + <command> + <db name> + <times>
  args = list(gr_args)
  if len(gr_args)>1 and gr_args[0] == 'wget':
    if args[1][-1] != '/': args[1] += '/'
    args[1] += "%s?name=%s"%(cmd,name)
    if cmd=='get_wrange': args[1]+="&t1=%s&t2=%s&dt=%f"%(t1,t2,dt)
    if cmd=='get_range':  args[1]+="&t1=%s&t2=%s&dt=%f"%(t1,t2,dt)
    if cmd=='get_next':   args[1]+="&t1=%s"%(t1)
    if cmd=='get_prev':   args[1]+="&t2=%s"%(t2)
    if cmd=='get':        args[1]+="&t2=%s"%(t2)
  else:
    args.extend((cmd, name))
    if cmd=='get_wrange': args.extend((t1,t2,"%f"%(dt)))
    if cmd=='get_range': args.extend((t1,t2,"%f"%(dt)))
    if cmd=='get_next':  args.append(t1)
    if cmd=='get_prev':  args.append(t2)
    if cmd=='get':       args.append(t2)

  ### run the command and load values
  proc = subprocess.Popen(args,
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True)
  data = proc.stdout.read()
  if proc.wait():
    print('> Graphene error:', proc.stderr.read())
    exit(1)
  return data

###############################################################
## Load data (with numpy.loadtxt).
## Suppress warnings on empty data.
## Problem with variable number of columns.
def graphene_load(ff, unpack=False, usecols=None):
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    data = numpy.loadtxt(ff, comments="#", usecols=usecols, unpack=unpack)
  return data

###############################################################
## Load data (without numpy.loadtxt)
## No problems with variable data length
def graphene_load2(ff, unpack=False, usecols=None):
  data=[]
  for x in ff.readlines():
    line = x.split()
    if len(line)==0: continue
    data.append(line)

  mlen=0
  # calculate max number of columns
  for x in data:
    if mlen<len(x): mlen=len(x)
  if usecols!=None and mlen < max(usecols)+1:
    mlen = max(usecols)+1

  # pad data to mlen
  for x in data:
    if len(x)<mlen: x += float('nan') * (mlen - len(x))

  # convert to numpy
  data = numpy.array(data, dtype='float')

  # take only needed indices
  if usecols!=None:
    data = numpy.take(data, usecols)

  # transpose if needed
  if usecols: data = numpy.transpose(data)
  return data

###############################################################
### Do arbitrary graphene command, read output, cache data
def graphene_cmd(cmd, name, t1="0", t2='inf', dt=0, unpack=False, usecols=None, fname=""):

  # convert timestammps to strings if needed
  if type(t1) != str: t1='%f'%(t1)
  if type(t2) != str: t2='%f'%(t2)

  t1 = timeconv(t1)
  t2 = timeconv(t2)

  ### do we want to use cache?
  if cache_dir != "":
    if not os.path.isdir(cache_dir): os.mkdir(cache_dir)
    if fname=="":
      fname = name + "_" + cmd + "_" + t1 + "_" + t2
    fname = cache_dir + "/" + fname
    if os.path.isfile(fname):
      ff = open(fname)
      return graphene_load2(ff, usecols=usecols, unpack=unpack)

  if cache_dir != "": ff = open(fname, "w+")
  else: ff = tempfile.TemporaryFile(mode="w+")

  data = graphene_read(cmd, name, t1, t2, dt)
  print(data, file=ff)
  ff.seek(0)

  return graphene_load2(ff, usecols=usecols, unpack=unpack)

###############################################################

def get_range(name, t1, t2, **kwargs):
  return graphene_cmd('get_range', name, t1=t1, t2=t2, **kwargs)

def get_wrange(name, t1, t2, **kwargs):
  return graphene_cmd('get_wrange', name, t1=t1, t2=t2, **kwargs)

def get_prev(name, t, **kwargs):
  return graphene_cmd('get_prev', name, t2=t, **kwargs)

def get_next(name, t, **kwargs):
  return graphene_cmd('get_next', name, t1=t, **kwargs)

def get(name, t, **kwargs):
  return graphene_cmd('get', name, t2=t, **kwargs)


