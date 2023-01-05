import subprocess
import numpy
import os.path

gr_args = ['device_c', 'ask', 'db']
cache_dir = ''

### Set program for accessing graphene, as separate arguments:
def set_args(a):
  global gr_args
  gr_args = list(a)

### Set cache folder (empty for no caching)
def set_cache(a):
  global cache_dir
  cache_dir = a

### Do arbitrary graphene command
def graphene_cmd(cmd, name, t1=0, t2='+inf', cols=(0,1), fname=""):

  ### do we want to use cache?
  if cache_dir != "":
    if not os.path.isdir(cache_dir): os.mkdir(cache_dir)
    if fname=="":
      fname = name + "_" + cmd + "_" + str(t1) + "_" + str(t2)
    fname = cache_dir + "/" + fname
    if os.path.isfile(fname):
      data = numpy.loadtxt(fname, comments="#", usecols=cols, unpack=True)
      return data

  ### build command: gr_args + <command> + <db name> + <times>
  args = list(gr_args)
  args.extend((cmd, name))

  if type(t1) != str: t1='%f'%(t1)
  if type(t2) != str: t2='%f'%(t2)

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
    data = numpy.loadtxt(fname, comments="#", usecols=cols, unpack=True)
    rc = proc.wait()
    if rc:
      print('> Graphene error:', proc.stderr.read())
      exit(1)
  return data


###

def get_range(name, t1, t2, cols=(0,1), fname=""):
  return graphene_cmd('get_range', name, t1=t1, t2=t2, cols=cols, fname=fname)

def get_prev(name, t, cols=(0,1), fname=""):
  return graphene_cmd('get_prev', name, t2=t, cols=cols, fname=fname)

def get_next(name, t, cols=(0,1), fname=""):
  return graphene_cmd('get_next', name, t1=t, cols=cols, fname=fname)

def get(name, t, cols=(0,1), fname=""):
  return graphene_cmd('get', name, t2=t, cols=cols, fname=fname)

