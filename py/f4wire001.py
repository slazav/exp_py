# processing vibrating wire data

import numpy
import scipy.optimize
import sys
import os
import graphene002 as graphene

###########################################################
# Wire thickness and length (projection to plane perpendicular to B), mm
def wire_dim(name):
  # Cell 2022
  if (name == 'w1ta2'): return (0.127, 3.61)  # measured on photo
  if (name == 'w1ta2'): return (0.127, 5.16)  # measured on photo
  if (name == 'w1bh'): return (13.5e-3, 2.74) # measured on photo
  if (name == 'w2bh'): return (13.5e-3, 2.58) # measured on photo
  if (name == 'w1bt'): return (4.5e-3, 1.49)  # measured on photo
  if (name == 'w2bt'): return (4.5e-3, 1.43) # measured on photo
  # thin wires of classical design
  if (name == 'w1a'): return (4.5e-3, 1.4) # unknown length
  if (name == 'w2a'): return (4.5e-3, 1.4) # unknown length
  if (name == 'w1b'): return (0.315e-3, 1) # unknown length
  if (name == 'w2b'): return (0.180e-3, 1) # unknown length
  # wires on PCB
  if (name == 'w1c'): return (0.390e-3, 1) # unknown length
  if (name == 'w2c'): return (0.315e-3, 1) # unknown length
  if (name == 'w1d'): return (0.180e-3, 1) # unknown length
  if (name == 'w2d'): return (0.180e-3, 1) # unknown length
  # Cell 2020
  if (name == 'w0ta'): return (0.127, 5.0)  # unknown length
  if (name == 'w1ta'): return (0.127, 4.87)  # measured on photo
  if (name == 'w2ta'): return (0.127, 5.13)  # measured on photo
  if (name == 'w0um'): return (4.5e-3, 1.4) # unknown length
  if (name == 'w1um'): return (4.5e-3, 1.4) # unknown length
  if (name == 'w2um'): return (4.5e-3, 1.4) # unknown length
  if (name == 'w1nm'): return (0.45e-3, 1) # unknown length
  if (name == 'w2nm'): return (0.45e-3, 1) # unknown length
  #
  if (name == 'mcta'): return (0.127, 5.0)  # unknown length
  print('ERROR: f4wire.wire_len: unknown name: ', name, file=sys.stderr)
  exit(1)

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
     cache="", plot="", nsweeps=1, nskip=0, prev_sweeps=1, fit_coord=0, fit_npars=6):
  import fit_res001 as fit_res

  if cache != "" and os.path.isfile(cache+".npz"):
    data = numpy.load(cache+".npz")
    return (data["arr_0"], data["arr_1"], data["arr_2"],
            data["arr_3"], data["arr_4"], data["arr_5"], data["arr_6"])

  # Get data with correct current and voltage
  data = get_data(name, t1, t2)

  # Get previous frequency sweep for thermometer and heater:
  if prev_sweeps:
    sweep = get_sweep_prev(name, t1, nsweeps=nsweeps, nskip=nskip)
    sweep = merge_sweeps(sweep, same_drive=0)[0]
  else:
    sweep = get_sweep_next(name, t2, nsweeps=nsweeps, nskip=nskip)
    sweep = merge_sweeps(sweep, same_drive=0)[0]

  # Fit the sweep
  fit = fit_res.fit(sweep, coord=fit_coord, npars=fit_npars)

  # Scale offset and amplitude to new drive
  TT = data[:,0]
  FF = data[:,1]
  DD = data[:,4]
  C = fit.C * DD/fit.drive
  D = fit.D * DD/fit.drive
  XX = data[:,2] - (fit.A + fit.E*(FF-fit.f0))*DD/fit.drive
  YY = data[:,3] - (fit.B + fit.F*(FF-fit.f0))*DD/fit.drive

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

  # Get field
  field = graphene.get_prev("demag_pc:f2", t1, usecols=1)[0][0]

  # Wire dimensions, mm (projection to plane perpendicular to B)
  (wd, wl) = wire_dim(name)

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
    vv=fit.func(ff)
    a.plot(ff, numpy.real(vv), 'k-', linewidth=1)
    a.plot(ff, numpy.imag(vv), 'k-', linewidth=1)
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

  if cache != "": numpy.savez(cache, TT, F0, dF, Vpar, Vperp, vel, PWR)
  return (TT, F0, dF, Vpar, Vperp, vel, PWR)

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
