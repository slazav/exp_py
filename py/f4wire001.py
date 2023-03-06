# processing vibrating wire data

import numpy
import scipy.optimize
import sys
import os
import graphene002 as graphene

###########################################################
# vibrating wire information

wire_info_tab = {
  # Cell 2020
  #         diameter, cm   length, cm    density, g/cm^3  cm/s
  'w1ta2': {'D': 127e-4,   'L': 3.61e-1, 'rho': 16.7}, # length measured on photo
  'w2ta2': {'D': 127e-4,   'L': 5.16e-1, 'rho': 16.7}, # length measured on photo
  'w1bh':  {'D': 13.5e-4,  'L': 2.74e-1, 'rho': 6.05}, # length measured on photo
  'w2bh':  {'D': 13.5e-4,  'L': 2.58e-1, 'rho': 6.05}, # length measured on photo
  # thin wires of classical design
  'w1bt':  {'D': 4.5e-4,   'L': 1.49e-1, 'rho': 6.05}, # length measured on photo
  'w2bt':  {'D': 4.5e-4,   'L': 1.43e-1, 'rho': 6.05}, # length measured on photo
  'w1a':   {'D': 4.5e-4,   'L': 1.4e-1,  'rho': 6.05}, # unknown length
  'w2a':   {'D': 4.5e-4,   'L': 1.4e-1,  'rho': 6.05}, # unknown length
  'w1b':   {'D': 0.315e-4, 'L': 1e-1,    'rho': 6.05}, # unknown length
  'w2b':   {'D': 0.180e-4, 'L': 1e-1,    'rho': 6.05}, # unknown length
  # wires on PCB
  'w1c':   {'D': 0.390e-4, 'L': 0.5e-1,  'rho': 6.05}, # unknown length
  'w2c':   {'D': 0.315e-4, 'L': 0.5e-1,  'rho': 6.05}, # unknown length
  'w1d':   {'D': 0.180e-4, 'L': 0.5e-1,  'rho': 6.05}, # unknown length
  'w2d':   {'D': 0.180e-4, 'L': 0.5e-1,  'rho': 6.05}, # unknown length
  # Cell 2020
  'w0ta':  {'D': 127e-4,   'L': 5.0e-1,  'rho': 16.7},  # unknown length
  'w1ta':  {'D': 127e-4,   'L': 4.87e-1, 'rho': 16.7},  # length measured on photo
  'w2ta':  {'D': 127e-4,   'L': 5.13e-1, 'rho': 16.7},  # length measured on photo
  'w0um':  {'D': 4.5e-4,   'L': 1.4e-1,  'rho': 6.05}, # unknown length
  'w1um':  {'D': 4.5e-4,   'L': 1.4e-1,  'rho': 6.05}, # unknown length
  'w2um':  {'D': 4.5e-4,   'L': 1.4e-1,  'rho': 6.05}, # unknown length
  'w1nm':  {'D': 0.45e-4,  'L': 1e-1,    'rho': 6.05}, # unknown length
  'w2nm':  {'D': 0.45e-4,  'L': 1e-1,    'rho': 6.05}, # unknown length
  #
  'mcta':  {'D': 127e-4,   'L': 5e-1,    'rho': 16.7},  # unknown length
}

# Intrinsic width measurements in vacuum, 2023-01-07:
# w1bh -2.789404*B^2 + 2.098990
# w2bh -4.708766*B^2 + 2.480987
# w1a 14.815948*B^2 + 0.767711
# w2a 4.852521*B^2 + 0.494533
# w1bt 9.220223*B^2 + 0.636275
# w2bt 33.285247*B^2 + 0.791410
# w1b 3.182436*B^2 + 0.358630
# w2d 101.593025*B^2 + 0.436765


## add B-phase non-linear parameters
####                                 cm/s          Hz             Hz/T^2
wire_info_tab['w1a'].update(  {'vmax': 0.14, 'dfi0': 0.292,  'dfi2': 14.82, 'S0': 1, 'S1': 6.7859, 'S2': -3.7249})
wire_info_tab['w1b'].update(  {'vmax': 0.72, 'dfi0': 0.0507, 'dfi2': 3.182, 'S0': 1, 'S1': 0.4699, 'S2':  0.0151})
wire_info_tab['w2a'].update(  {'vmax': 0.39, 'dfi0': 0.0810, 'dfi2': 4.853, 'S0': 1, 'S1': 1.0609, 'S2':  0.1434})
wire_info_tab['w1bt'].update( {'vmax': 0.45, 'dfi0': 0.1664, 'dfi2': 9.220, 'S0': 1, 'S1': 1.2103, 'S2': 0.1638})
wire_info_tab['w2bt'].update( {'vmax': 0.19, 'dfi0': 0.2492, 'dfi2': 33.29, 'S0': 1, 'S1': 3.4182, 'S2': 0.1825})
wire_info_tab['w1bh'].update( {'vmax': 0.32, 'dfi0': 0.3314, 'dfi2': 0.000, 'S0': 1, 'S1': 3.2977, 'S2': -0.8180})
wire_info_tab['w2bh'].update( {'vmax': 0.29, 'dfi0': 0.4377, 'dfi2': 0.000, 'S0': 1, 'S1': 2.2225, 'S2': -0.3893})
wire_info_tab['w1ta2'].update({'vmax': 0.49, 'dfi0': 0.0000, 'dfi2': 0.000, 'S0': 1, 'S1': 0.4508, 'S2': 0.1005})
wire_info_tab['w2ta2'].update({'vmax': 0.55, 'dfi0': 0.1966, 'dfi2': 0.000, 'S0': 1, 'S1': 0.5634, 'S2': 0.0383})

class wire_info_t:
  D = 0    # diameter [cm]
  L = 0    # length (projection to plane perpendicular to B) [cm]
  rho = 0  # material density [g/cm^3].

  # B-phase non-linear parameters:
  vmax = 0 # max velocity for the non-linear model [cm/s]
  dfi0 = 0 # intrinsic width at zero field [Hz]
  dfi2 = 0 # field-dependent part of intrinsic width [Hz/T^2]
  S0 = 1   # S-function parameters
  S1 = 0   #
  S2 = 0   #

  # 4-th power polyfit of  pF^2 * vF * 2N(0)/2 [sgs] vs P [bar]
  pp_d = (-1.6016e-03, 1.2582e-01, -3.9976e+00, 8.8950e+01, 2.0300e+03)
  # 4-th power polyfit of  kTc/pF [cm/s] vs P [bar]
  pp_v = (-2.9402e-06, 2.5419e-04, -9.4237e-03, 2.0209e-01, 1.5584e+00)
  # 4-th power polyfit of gap/kTc vs P [bar]
  pp_g = (-1.8629e-07, 1.4784e-05, -4.6497e-04, 8.8839e-03, 1.7725e+00)

  #####################
  # correction function S
  def sfunc(self, P, B, ttc, vel=None, volt=None):
    if type(volt)!= type(None):
      vel = volt/B/self.L * 1e4  # V/T/cm -> cm/s
    if type(vel) == type(None):
      raise Exception("wire_info_t::delta: missing vel or volt parameter")
    vv = vel / (numpy.polyval(self.pp_v, P)*ttc)  # v/v0
    return self.S0/(1 + self.S1*vv + self.S2*vv**2)

  def ttc_to_delta0(self, P, ttc):
    if ttc <= 0: return 0
    gap = numpy.polyval(self.pp_g, P)
    d0  = 2*numpy.polyval(self.pp_d, P)/self.rho/self.D / (2*numpy.pi) #[sgs]
    return d0 * numpy.exp(-gap/ttc) 

  def delta0_to_ttc(self, P, delta0):
    if delta0 <= 0: return 0
    gap = numpy.polyval(self.pp_g, P)
    d0  = 2*numpy.polyval(self.pp_d, P)/self.rho/self.D / (2*numpy.pi) #[sgs]
    return -gap/numpy.log(delta0/d0)

  #####################
  # calibration delta(P,B,V,ttc)
  def ttc_to_delta(self, P, B, ttc, vel=None, volt=None):
    dfi = self.dfi0 + self.dfi2 * B**2
    delta0 = self.ttc_to_delta0(P, ttc)
    S = self.sfunc(P,B, ttc, vel=vel, volt=volt)
    return dfi + delta0*S

  #####################
  # Inversed calibration, ttc(P,B,V,delta)
  def delta_to_ttc(self, P, B, delta, vel=None, volt=None):
    dfi = self.dfi0 + self.dfi2 * B**2
    gap = numpy.polyval(self.pp_g, P)
    d0 = 2*numpy.polyval(self.pp_d, P)/self.rho/self.D / (2*numpy.pi) #[sgs]

    ttcp = 0.2 # starting point for iterations
    for i in range(20):
      if delta<=dfi or ttcp==0:
        S = 1
        ttc = 0
      else:
        S = self.sfunc(P,B,ttcp,vel=vel,volt=volt)
        ttc = -gap/numpy.log((delta - dfi)/d0/S)
      if (abs(ttc-ttcp) < 1e-6): break
      else: ttcp = ttc
    return ttc

  #####################
  # Correction, delta0(P,B,V,delta)
  def delta0(self, P, B, delta, vel=None, volt=None):
    ttc = self.ttc(P, B, delta, vel=None, volt=None)
    dfi = self.dfi0 + self.dfi2 * B**2
    return (delta-dfi)/self.sfunc(P,B,ttc,vel=vel,volt=volt)

  #####################
  def __init__(self, name):
    if name not in wire_info_tab:
      raise Exception("'ERROR: f4wire: unknown name: ' + name")
    w = wire_info_tab[name]
    if 'D'   in w: self.D   = w['D']
    if 'L'   in w: self.L   = w['L']
    if 'rho' in w: self.rho = w['rho']
    if 'dfi0' in w: self.dfi0 = w['dfi0']
    if 'dfi2' in w: self.dfi2 = w['dfi2']
    if 'S0'   in w: self.S0   = w['S0']
    if 'S1'   in w: self.S1   = w['S1']
    if 'S2'   in w: self.S2   = w['S2']
    return

###########################################################
# Wire thickness and length (projection to plane perpendicular to B), mm
def wire_dim(name):
  if name not in wire_info_tab:
    raise Exception("'ERROR: f4wire: unknown name: ' + name")
  w = wire_info_tab[name]
  return (w['D']*10, w['L']*10)

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
      data[:,2] /= v2v[0]
      data[:,3] /= v2v[0]
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
      sweeps.append(numpy.array(s, dtype=float))
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

def get_sweep_prev(name, t1, nsweeps=1, nskip=0, **kwargs):
  sweeps = []
  for i in range(nskip):
    pars = graphene.get_prev(name + '_pars', t1)
    t1 = pars[0][0]-1e-6
  for i in range(nsweeps):
    pars = graphene.get_prev(name + '_pars', t1)
    sweeps.append(get_sweeps_(name, pars, **kwargs)[0])
    t1 = pars[0][0]-1e-6
  return sweeps

def get_sweep_next(name, t1, nsweeps=1, nskip=0, **kwargs):
  sweeps = []
  for i in range(nskip):
    pars = graphene.get_next(name + '_pars', t1)
    t1 = pars[0][0]-1e-6
  for i in range(nsweeps):
    pars = graphene.get_next(name + '_pars', t1)
    sweeps.append(get_sweeps_(name, pars, **kwargs)[0])
    t1 = pars[0][0]+1e-6
  return sweeps

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

###########################################################
# merge sweeps (only merge same drive or merge all with amplitude rescaling)
def merge_sweeps(sweeps, same_drive=1):

  if same_drive==0:
    ret = numpy.row_stack(sweeps)
    ret[:,2] *= ret[0,4] / ret[:,4]
    ret[:,3] *= ret[0,4] / ret[:,4]
    return [ret,]

  ret = []
  for s in sweeps:
    if len(ret)>0 and ret[-1][0,4] == s[0,4]:
      ret[-1] = numpy.row_stack((ret[-1], s))
    else:
      ret.append(s)
  return ret

###########################################################
# Process tracking mode data. Assuming that resonance is linear,
# and A,B,C,D,E,F parameters are not changing (or scaled with drive)
# find f0 and df parameters.
def track_res_lin(data, fit):
  # scale to original drive, subtract offset
  FF = data[:,1]
  XX = data[:,2]*fit.drive/data[:,4] - fit.A - fit.E*(FF-fit.f0)
  YY = data[:,3]*fit.drive/data[:,4] - fit.B - fit.F*(FF-fit.f0)

  # coord:     X + i*Y = (C + i*D) / (f0^2-f^2 + i*f*df)
  # velocity:  X + i*Y = 1j*F*(C + i*D) / (f0^2-f^2 + i*f*df)
  # find V = f0^2-f^2 + i*f*df:
  VV = (fit.C + 1j*fit.D)/(XX + 1j*YY)
  if not fit.coord: VV *= 1j*FF

  F0 = numpy.sqrt(numpy.real(VV) + FF**2)
  dF = numpy.imag(VV)/FF
  return (F0, dF)

###########################################################
# Process tracking mode data.
# Calculate dissipated power: drive current multiplied by in-phase voltage.
def track_heat(data, fit):

  # subtract offset scaled to new drive
  FF = data[:,1]
  XX = data[:,2] - (fit.A + fit.E*(FF-fit.f0))*data[:,4]/fit.drive
  YY = data[:,3] - (fit.B + fit.F*(FF-fit.f0))*data[:,4]/fit.drive

  # Project (X,Y) to (C,D) in coord mode, (-D,C) in velocity mode.
  # C and D are not scaled to drive, but absolute value is not important.
  # We use drive for that
  if fit.coord:
#    Vpar = (fit.C*XX + fit.D*YY)/numpy.hypot(fit.C,fit.D)
    Vperp = (-fit.C*YY + fit.D*XX)/numpy.hypot(fit.C,fit.D)
  else:
#    Vpar = (fit.C*YY - fit.D*XX)/numpy.hypot(fit.C,fit.D)
    Vperp = (fit.C*XX + fit.D*YY)/numpy.hypot(fit.C,fit.D)

  return data[:,4]*Vperp



###########################################################
# Process tracking mode data.
def get_track(name, t1, t2,
     get=0, cache="", plot="", nsweeps=1, nskip=0, prev_sweeps=1,
     fit_coord=0, fit_npars=6, bphase=None):
  import fit_res002 as fit_res

  if cache != "" and get==0 and os.path.isfile(cache+".npz"):
    data = numpy.load(cache+".npz")

    # old cache files have less data
    if 'arr_7' in data:
      return (data["arr_0"], data["arr_1"], data["arr_2"],
              data["arr_3"], data["arr_4"], data["arr_5"],
              data["arr_6"], data["arr_7"])

  # Get data with correct current and voltage
  data = get_data(name, t1, t2)

  # Get field
  field = graphene.get_prev("demag_pc:f2", t1, usecols=1)[0][0]

  # Wire dimensions, mm (projection to plane perpendicular to B)
  (wd, wl) = wire_dim(name)

  # Get previous frequency sweep for thermometer and heater:
  if prev_sweeps:
    sweep = get_sweep_prev(name, t1, nsweeps=nsweeps, nskip=nskip)
    sweep = merge_sweeps(sweep, same_drive=0)[0]
  else:
    sweep = get_sweep_next(name, t2, nsweeps=nsweeps, nskip=nskip)
    sweep = merge_sweeps(sweep, same_drive=0)[0]

  # Fit the sweep
  fit = fit_res.fit(sweep, coord=fit_coord, npars=fit_npars, bphase=bphase)

  # Scale offset and amplitude to new drive
  TT = data[:,0]
  FF = data[:,1]
  DD = data[:,4]
  C = fit.C * DD
  D = fit.D * DD
  XX = data[:,2] - (fit.A + fit.E*(FF-fit.f0))*DD
  YY = data[:,3] - (fit.B + fit.F*(FF-fit.f0))*DD

  # Resonance frequency and width:
  # coord:     X + i*Y = (C + i*D) / (f0^2-f^2 + i*f*df)
  # velocity:  X + i*Y = 1j*F*(C + i*D) / (f0^2-f^2 + i*f*df)
  # find V = f0^2-f^2 + i*f*df:
  VV = (C + 1j*D)/(XX + 1j*YY)
  if not fit.coord: VV *= 1j*FF
  F0 = numpy.sqrt(numpy.real(VV) + FF**2)
  dF = numpy.imag(VV)/FF

  # Project (X,Y) to (C,D) in coord mode, (-D,C) in velocity mode.
  if fit.coord:
    Vpar = (C*XX + D*YY)/numpy.hypot(C,D)
    Vperp = (-C*YY + D*XX)/numpy.hypot(C,D)
  else:
    Vpar = (C*YY - D*XX)/numpy.hypot(C,D)
    Vperp = (C*XX + D*YY)/numpy.hypot(C,D)

  # Power
  PWR = DD*Vperp

  # Velocity
  vel=numpy.hypot(Vpar, Vperp)/field/(wl*1e-3)

  # do plot if needed
  if plot!="":
    import matplotlib.pyplot as plt
    (fig, ax) = plt.subplots(2,2)

    # sweep
    a=ax[0,0]
    a.plot(sweep[:,1], sweep[:,2], 'r.', label="X")
    a.plot(sweep[:,1], sweep[:,3], 'b.', label="Y")
    ff=numpy.linspace(min(sweep[:,1]), max(sweep[:,1]), 100)
    vv1=fit.func(ff, numpy.min(sweep[:,4]))
    vv2=fit.func(ff, numpy.max(sweep[:,4]))
    a.plot(ff, numpy.real(vv1), 'k-', linewidth=1)
    a.plot(ff, numpy.imag(vv2), 'k-', linewidth=1)
    a.plot(ff, numpy.real(vv1), 'k-', linewidth=1)
    a.plot(ff, numpy.imag(vv2), 'k-', linewidth=1)
    a.set_xlabel("freq, Hz")
    a.set_ylabel("volt, Vrms")
    a.set_title("frequency sweep")

    t0 = TT[0]
    # width and frequency in tracking mode
    a1=ax[0,1]
    a2=a1.twinx()
    a2.plot(TT-t0, F0, 'b.-', label="f0rack")
    a2.plot(TT-t0, FF, 'g.-', label="f_meas")
    a1.plot(TT-t0, dF, 'r.-', label="dfrack")
    xx=[0, TT[-1]-t0]
    a1.plot(xx, [fit.df]*2, 'm-', label='df_fit')
    a2.plot(xx, [fit.f0]*2, 'c-', label='f0_fit')
    a1.set_xlabel("time, s")
    a1.set_ylabel("df, Hz")
    a2.set_ylabel("f0, Hz")
    a1.legend()
    a2.legend()
    a1.set_title("freq, width")

    # voltage components
    a=ax[1,0]
    a.plot(TT-t0, Vpar,  'r.-', label="Vpar")
    a.plot(TT-t0, Vperp, 'b.-', label="Vperp")
    a.set_xlabel("time, s")
    a.set_ylabel("Voltage")
    a.legend()
    a.set_title("X,Y")

    # power
    a=ax[1,1]
    a.semilogy(TT-t0, PWR,  'r.-')
    a.set_xlabel("time, s")
    a.set_ylabel("Power, W")
    a.set_title("power")

    # save plot
    plt.gcf().set_size_inches(12, 12)
    plt.savefig(plot, dpi=100)
    plt.close()

  if cache != "": numpy.savez(cache, TT, F0, dF, Vpar, Vperp, vel, PWR, field)
  return (TT, F0, dF, Vpar, Vperp, vel, PWR, field)

###########################################################
# Process background sweep.
def get_bg(nn, t1, fmin=0, fmax=1e6, fminres=0, fmaxres=0, cache="", plot=""):

  # get data for sweep, no conversions!
  sweeps = get_sweep(nn, t1, cache=cache, use_bg=0, cnv_volt=0, cnv_drive=0)

  if len(sweeps) < 1: return
  s = sweeps[0]
  print("\n", nn)

  F = s[:,1] *1e-3; # kHz
  X = s[:,2]/s[:,4] * 1e6; # uV/V drive
  Y = s[:,3]/s[:,4] * 1e6; # uV/V drive

  ii=numpy.logical_and(F>=fmin/1000, F<=fmax/1000)
  if fminres!=0 and fmaxres!=0:
    jj=numpy.logical_or(F<fminres/1000, F>fmaxres/1000)
    ii=logical_and(ii,jj)

  # 6-prameter function compatible with background setting in
  # forc_cw interface
  px=(1000,1000,1,1,5,5)
  py=(1000,1000,1,1,5,5)
  def fitfunc(x, a,b,c,d, f,g):
    return (a+b*x+c*x**2+d*x**3)/((x**2-f**2)**2 + (g*x)**2)

  px = scipy.optimize.curve_fit(fitfunc, F[ii], X[ii], px, maxfev=200000)[0]
  py = scipy.optimize.curve_fit(fitfunc, F[ii], Y[ii], py, maxfev=200000)[0]

  print("%f %f %f %f %f %f"%(*px,))
  print("%f %f %f %f %f %f"%(*py,))

  if plot!="":
    import matplotlib.pyplot as plt
    xx=numpy.linspace(0.5,10, 200)
    plt.clf()
    plt.plot(F,X, 'r.')
    plt.plot(F,Y, 'b.')
    plt.plot(xx, fitfunc(xx, *px), 'k-')
    plt.plot(xx, fitfunc(xx, *py), 'k-')

    plt.xlabel('freq, kHz')
    plt.ylabel('X/D, Y/D [uV/V]')

    fig = plt.gcf()
    fig.set_size_inches(8, 6)
    plt.savefig(plot + ".png", dpi=100)
    plt.close()
  return (*px, *py)
