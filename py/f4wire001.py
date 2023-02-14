# processing vibrating wire data

import numpy
import sys
import os
import graphene002 as graphene

###########################################################
# Calculate background using standard 12-parameter model.
#
def calc_bg(bg, x, im=0):
  if bg.size==0:
    return numpy.zeros_like(x)
  if bg.size!=12:
    print("ERROR: background contains < 12 numbers:\n", bg, file=sys.stderr)
    return numpy.zeros_like(x)
  if im: sh = 6
  else: sh = 0
  return 1e-6*(bg[sh+0] + bg[sh+1]*x/1000 + bg[sh+2]*(x/1000)**2 + bg[sh+3]*(x/1000)**3)/\
              (((x/1000)**2-bg[sh+4]**2)**2 + (bg[sh+5]*x/1000)**2)

###########################################################
# Base function for getting data from *_sweeps databases,
# subtracting background, converting voltage and current.
# Works for both frequency sweeps and tracking data.
#
def get_data(name, t1, t2, use_bg=1, cnv_drive=1, cnv_volt=1, cache=""):

  if cache != "" and os.path.isfile(cache):
    return numpy.loadtxt(cache)

  # data: T,F,X,Y,D
  data = graphene.get_range(name + "_sweeps", t1, t2)

  if use_bg:
    bg = graphene.get_prev(name + "_dbox:f2", t1)
    if bg.size>0:
      data[:,2] -= calc_bg(bg[0][1:], data[:,1], 0) * data[:,4]
      data[:,3] -= calc_bg(bg[0][1:], data[:,1], 1) * data[:,4]
      #print("bg: ", bg[0])

  if cnv_drive:
    v2i = graphene.get_prev(name + "_dbox:f1", t1, usecols=1)
    if v2i.size>0:
      data[:,4] *= v2i[0]
      #print("v2i: ", v2i[0])

  if cnv_volt:
    v2v = graphene.get_prev(name + "_meas:f1", t1, usecols=1)
    if v2v.size>0:
      data[:,4] /= v2v[0]
      #print("v2v: ", v2v[0])

  if cache != "": numpy.savetxt(cache, data)
  return data

###########################################################
# Get sweeps using parameters from <name>_pars database
def get_sweeps_(name, pars, sweep_dir=None, cache="", **kwargs):
  sweeps=[]
  # no sweeps
  if pars.size == 0: return sweeps

  # load cache if needed
  if cache != "" and os.path.isfile(cache):
    ff=open(cache, "r")
    s = []
    while 1:
      l = ff.readline()
      if not l: break
      l = l.split()
      if len(l) == 0:
        if len(s)>0:
          sweeps.append(numpy.array(s, dtype=float))
          s=[]
      elif len(l) == 5:
        s.append(l)
      else:
        print("ERROR: wrong data in cache file: ", cache, l)
    if len(s)>0:
      sweeps.append(numpy.array(s))
      s=[]
    return sweeps

  # single sweep
  if len(pars.shape) == 1:
    pars = numpy.reshape(pars, (1,-1))

  for i in range(pars.shape[0]):
    (t1, dt) = pars[i,0:2]
    if sweep_dir ==  1 and pars[i,7] !=  1: continue
    if sweep_dir == -1 and pars[i,7] != -1: continue
    # Get sweep points
    sweep = get_data(name, t1, t1+dt, **kwargs)
    sweeps.append(sweep)

  # save cache if needed
  if cache != "":
    ff=open(cache, "w")
    for s in sweeps:
      for i in range(s.shape[0]):
        print("%.6f %e %e %e %e"%(*s[i,:],), file=ff)
      print("", file=ff)

  return sweeps

###########################################################
# Get sweeps

def get_sweep_prev(name, t1, **kwargs):
  pars = graphene.get_prev(name + '_pars', t1)
  return get_sweeps_(name, pars, **kwargs)

def get_sweep_next(name, t1, **kwargs):
  pars = graphene.get_prev(name + '_pars', t1)
  return get_sweeps_(name, pars, **kwargs)

def get_sweep_range(name, t1, t2, **kwargs):
  pars = graphene.get_range(name + '_pars', t1, t2)
  return get_sweeps_(name, pars, **kwargs)

# sweep covering time t1
def get_sweep(name, t1, **kwargs):
  t1 = graphene.timeconv(t1)
  pars = graphene.get_prev(name + '_pars', t1)
  if pars.size==0: return ()
  pars=pars[0]
  # Note: special timestamps (inf, now, <t>+) will fail here!
  if pars.size==0 or pars[0] + pars[1] < float(t1): return (); # too old
  return get_sweeps_(name, pars, **kwargs )

# list of sweeps
def get_sweep_list(name, tlist, **kwargs):
  pars = numpy.empty((0,10))
  for t in tlist:
    t = graphene.timeconv(t)
    p = graphene.get_prev(name + '_pars', t)
    # Note: special timestamps (inf, now, <t>+) will fail here!
    if p.size==0: continue
    p=p[0]
    if p[0] + p[1] < float(t): continue; # too old
    pars = numpy.append(pars, numpy.reshape(p,(1,10)), axis=0)
  return get_sweeps_(name, pars, **kwargs)
