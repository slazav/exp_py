# Processing tracking mode data.
# This is a new version which uses arbitrary D-function
# delta(delta0, |vel|)
# Old version is still in f4wire001.py
# V.Zavjalov, 14.06.2023

import numpy
import scipy.optimize
import sys
import os
import graphene002 as graphene
import f4wire001 as f4wire
import fit_res003 as fit_res


def get_track(name, t1, t2,
     get=0, cache="", plot="", nsweeps=1, nskip=0, prev_sweeps=1,
     fit_coord=0, fit_npars=6, dfunc=None, verb=0, osc=""):

  ############################
  # get data from DB or cache
  if cache != "" and get==0 and os.path.isfile(cache+".npz"):
    if verb: print("get_track: load cache: ", cache)
    d = numpy.load(cache+".npz")
    (data, sweep, field, press) = (d["arr_0"], d["arr_1"], d["arr_2"], d["arr_3"])

  else:
    if verb: print("get_track: get data from db")
    # Get data with correct current and voltage
    if osc=="":
      data = f4wire.get_data(name, t1, t2, use_bg=1, cnv_volt=1, cnv_drive=1)
    else:
      data = f4wire.get_data_osc(name, osc, use_bg=1, cnv_volt=1, cnv_drive=1)

      # if t1, t2 = None, use full time span of the oscilloscope file
      # if not, use them as relative values calculated from beginning of the file
      if not t1: t1 = data[0,0]
      else:  t1 = data[0,0] + t1

      if not t2: t2 = data[-1,0]
      else: t2 = data[0,0] + t2

      ii = numpy.logical_and(data[:,0]>=t1, data[:,0]<=t2)
      data=data[ii,:]

    # Get previous frequency sweep for thermometer and heater:
    if prev_sweeps:
      sweep = f4wire.get_sweep_prev(name, t1, nsweeps=nsweeps, nskip=nskip)
      sweep = numpy.row_stack(sweep)
    else:
      sweep = f4wire.get_sweep_next(name, t2, nsweeps=nsweeps, nskip=nskip)
      sweep = numpy.row_stack(sweep)

    # Get field
    field = graphene.get_prev("demag_pc:f2", t1, usecols=1)[0][0]

    # Get pressure
    press = graphene.get_prev("cell_press", t1, usecols=1)[0][0]

    if cache != "":
      numpy.savez(cache, data, sweep, field, press)

  ############################

  if data.size == 0: return None

  # Wire information
  wire = f4wire.wire_info_t(name)

  # Fit the sweep
  fit = fit_res.fit(sweep, coord=fit_coord, npars=fit_npars, dfunc=dfunc)

  # Scale amplitude to new drive, remove offset
  TT = data[:,0]
  FF = data[:,1]
  DD = data[:,4]
  fit = fit
  C = fit.C * DD
  D = fit.D * DD
  XX = data[:,2] - (fit.A + fit.E*(FF-fit.f0))*DD
  YY = data[:,3] - (fit.B + fit.F*(FF-fit.f0))*DD

  # Return object
  class ret_t: pass
  ret = ret_t()
  ret.TT = TT
  ret.FF = FF
  ret.DD = DD
  ret.field=field
  ret.press=press
  ret.sweep=sweep
  ret.fit = fit

  # Resonance frequency and width:
  # coord:     X + i*Y = (C + i*D) / (f0^2-f^2 + i*f*df)
  # velocity:  X + i*Y = 1j*F*(C + i*D) / (f0^2-f^2 + i*f*df)
  # find V = f0^2-f^2 + i*f*df:
  VV = (C + 1j*D)/(XX + 1j*YY)
  if not fit.coord: VV *= 1j*FF
  ret.F0 = numpy.sqrt(numpy.real(VV) + FF**2)
  ret.dF = numpy.imag(VV)/FF

  # Project (X,Y) to (C,D) in coord mode, (-D,C) in velocity mode.
  CD = numpy.hypot(C,D)
  if fit.coord:
    ret.VX = (C*XX + D*YY)/CD
    ret.VY = (-C*YY + D*XX)/CD
  else:
    ret.VX = (C*YY - D*XX)/CD
    ret.VY = (C*XX + D*YY)/CD

  # Power
  ret.PWR = DD*ret.VY

  # Velocity [cm/s]
  ret.volt = numpy.hypot(XX,YY)
  ret.volt2vel = 1/field/(wire.L*1e-2) * 100
  ret.vel = ret.volt*ret.volt2vel

  # do plot if needed
  if plot!="":
    import matplotlib.pyplot as plt
    (fig, ax) = plt.subplots(2,2)

    # sweep
    a=ax[0,0]
    fit_res.plot(a, a, sweep, fit, npts=200, xlabel="X", ylabel="Y")
    a.set_ylabel("volt, Vrms")
    a.set_title("frequency sweep")

    t0 = TT[0]
    # width and frequency in tracking mode
    a1=ax[0,1]
    a2=a1.twinx()
    a2.plot(TT-t0, ret.F0, 'b-', label="f0 track")
    a2.plot(TT-t0, FF, 'g-', label="f meas")

    a1.plot(TT-t0, ret.dF, 'r-', label="df track")
    if bphase:
      a1.plot(TT-t0, ret.dF0, 'm-', label="dF corr")

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
    a.plot(TT-t0, ret.VX, 'r.-', label="Vpar")
    a.plot(TT-t0, ret.VY, 'b.-', label="Vperp")
    a.set_xlabel("time, s")
    a.set_ylabel("Voltage")
    a.legend()
    a.set_title("X,Y")

    # power
    a=ax[1,1]
    a.semilogy(TT-t0, ret.PWR,  'r.-')
    a.set_xlabel("time, s")
    a.set_ylabel("Power, W")
    a.set_title("power")

    # save plot
    plt.gcf().set_size_inches(12, 12)
    plt.savefig(plot, dpi=100)
    plt.close()

  return ret
