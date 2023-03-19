import subprocess
import warnings
import numpy
import re
import tempfile
import os.path
import time
import math
import sys

sources = {
  "local":   ('device_c', 'ask', 'db'),
  "ssh_f4b": ('ssh', 'f4b', 'device_c', 'ask', 'db'),
  "xyz_f4":  ('wget', 'http://slazav.xyz:8095/', '-O', '-', '-o', '/dev/null'),
  "xyz_f2":  ('wget', 'http://slazav.xyz:8091/', '-O', '-', '-o', '/dev/null'),
  "xyz_dd":  ('wget', 'http://slazav.xyz:8085/', '-O', '-', '-o', '/dev/null'),
}

gr_args = ['device_c', 'ask', 'db']

### Set program for accessing graphene
### - list of arguments or name from sources table
def set_source(a):
  global gr_args
  if type(a) == list or type(a) == tuple:
    gr_args = list(a)
  else:
    gr_args = sources[a]


###############################################################
### Convert time in human-readable form for graphene
def timeconv(t, fmt=""):

  # convert timestamp to string if needed
  if type(t) != str: t='%f'%(t)

  if t=="now" or t=="now_s" or t=="inf": return t
  if re.fullmatch('[0-9\.]+[+-]?', t): return t

  if fmt!="":
    t = time.mktime(time.strptime(t, fmt))
    return "%.6f"%(t)

  for fmt in ('%Y-%m-%d %H:%M:%S',
              '%Y-%m-%dT%H:%M:%S',
              '%Y-%m-%d %H:%M',
              '%Y-%m-%dT%H:%M',
              '%Y-%m-%d %H',
              '%Y-%m-%dT%H',
              '%Y-%m-%d'):
    try:
      t = time.mktime(time.strptime(t, fmt))
      return "%.6f"%(t)
    except ValueError:
      pass
  raise ValueError('no valid date format found for ', t)


###############################################################
# Communication with graphene, read-only operations.
# No caching, no data parsing.
# t1 and t2 should be strings.
def graphene_run(cmd, name, t1="0", t2='inf', dt=0, verb=1):

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
  if verb>0:
    print("Running command: ", " ".join(args), file=sys.stderr)
  proc = subprocess.Popen(args,
      stdout=subprocess.PIPE,
      stderr=subprocess.STDOUT,
      text=True)
  data = proc.stdout.read()
  if proc.wait():
    raise Exception('> Graphene error:', data)

  return data

###############################################################
## Load data (without numpy.loadtxt)
## No problems with variable data length
## raw -- return python list of lists with raw string data instead if numpy array
def graphene_load(ff, unpack=False, usecols=None, raw=0):
  data=[]
  if usecols != None and\
     not isinstance(usecols, list) and\
     not isinstance(usecols, tuple): usecols=(usecols,)

  for x in ff.readlines():
    line = x.split()
    if len(line)==0:
      continue

    if usecols == None:
      data.append(line)
      continue

    d=[]
    for c in usecols:
      if c<len(line): d.append(line[c])
      else: d.append("nan")
    data.append(d)

  if raw: return data


  mlen=0
  # calculate max number of columns
  for x in data:
    if mlen<len(x): mlen=len(x)

  # pad data to mlen
  for x in data:
    if len(x)<mlen: x += [float('nan')] * (mlen - len(x))

  # convert to numpy
  data = numpy.array(data, dtype='float')
  if data.size==0: return data

  # transpose if needed
  if unpack: data = numpy.transpose(data)
  return data

###############################################################
### Do arbitrary graphene command, read output, cache data
def graphene_cmd(cmd, name, t1="0", t2='inf',
                 dt=0, unpack=False, usecols=None, raw=False, cache="", verb=1):

  t1 = timeconv(t1)
  t2 = timeconv(t2)

  ### do we want to use cache?
  if cache != "":
    if os.path.isfile(cache):
      if verb>0:
        print("Using cache: ", cache, file=sys.stderr)
      ff = open(cache)
      return graphene_load(ff, usecols=usecols, unpack=unpack, raw=raw)
    ff = open(cache, "w+")
  else: ff = tempfile.TemporaryFile(mode="w+")

  data = graphene_run(cmd, name, t1, t2, dt, verb=verb)
  print(data, file=ff)
  ff.seek(0)

  return graphene_load(ff, usecols=usecols, unpack=unpack, raw=raw)

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


