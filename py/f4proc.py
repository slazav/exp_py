# random functions for data processing

import numpy
import scipy.optimize
import matplotlib.pyplot as plt

import os
import sys
import f4wire001 as f4wire
import fit_res002 as fit_res

# Get data and do basic plot for full sweep measurement
def plot_drive_sweeps(t1,t2,names, shift=0.2, nmin=None, ncnt=None):
  os.makedirs('data', exist_ok=1)
  for nn in names:
    sweeps = f4wire.get_sweep_range(nn, t1, t2, cache='data/%s.dat'%(nn))

    if nmin: sweeps = sweeps[nmin:]
    if ncnt: sweeps = sweeps[:ncnt]

    if len(sweeps)<1: continue

    # get max amplitude
    m0=0
    for s in sweeps:
      m = numpy.max(abs(s[:,2]/s[:,4]))
      if m0<m: m0=m
      m = numpy.max(abs(s[:,3]/s[:,4]))
      if m0<m: m0=m
    print(m0)

    # fit first sweep to get phase
    fit = fit_res.fit(sweeps[0], coord=0)
    ph = numpy.arctan2(fit.D, fit.C)
    d0 = sweeps[0][0,4]
    print(nn, ' ph: ', ph)

    # plot sweeps
    (fig, ax) = plt.subplots(1,3)
    sh=0
    dp=None
    for s in sweeps:
      d = s[0,4]
      vv = (s[:,2] + 1j*s[:,3])/d - (fit.A + 1j*fit.B)*d0/d
      vv *= numpy.exp(-1j*ph)

      ax[0].plot(s[:,1], sh+numpy.real(vv), '-')
      ax[1].plot(s[:,1], sh+numpy.imag(vv), '-')
      ax[2].plot(s[:,1], sh+numpy.abs(vv), '-')

      if d!=dp: sh+=m0*shift
      dp = d

    # plot fit
    ll=plt.xlim()
    ff=numpy.linspace(ll[0], ll[1], 100)
    vv=(fit.func(ff, 1) - (fit.A + 1j*fit.B))*numpy.exp(-1j*ph)
    ax[0].plot(ff, numpy.real(vv), 'k-', linewidth=0.7)
    ax[1].plot(ff, numpy.imag(vv), 'k-', linewidth=0.7)
    ax[2].plot(ff, numpy.abs(vv), 'k-', linewidth=0.7)

    plt.gcf().set_size_inches(10, 8)
    plt.savefig(nn + ".png", dpi=100)
    plt.close()

######################################################################

# Process bolometer heating measurements (steps)
# - take osc file, do usual tracking mode processing, get width or corrected width, or temperature
# - take heater data, find power
# - find positions where power is switching, find time constants
# - find temperature vs power dependence,
#
# Arguments:
# - fname    -- oscilloscope file name (without .osc extension)
# - name     -- thermometer wire name, 'w1bt' or 'w2bt'
# - t1,t2    -- get this time range from the oscilloscope file
# - dtm      -- time before switching heater to be used in fitting
# - dtp      -- time after switching heater to be used in fitting
# - htr_thr  -- threshold power for the heater
def process_bolo_heat(fname, name, t1 = 20, t2 = 320, dtm = 10, dtp = 15, htr_thr = 1e-14):
  os.makedirs('data_bolo', exist_ok=1)
  # get thermometer data
  data1 = f4wire.get_track(name, t1, t2, use_bphase=1,
    osc=fname + '.sig', cache="data_bolo/" + fname, plot="data_bolo/" + fname + '.png')

  if name == 'w1bt': heater='w1bh'
  if name == 'w2bt': heater='w2bh'
  # get heater data
  data2 = f4wire.get_track(heater, data1.TT[0], data1.TT[-1], use_bphase=1,
    cache="data_bolo/" + fname + '_htr', plot="data_bolo/" + fname + '_htr.png')

  # get indices of heater changes
  ihtr=[]
  for i in range(1,data2.PWR.size):
    if data2.PWR[i] >= htr_thr and data2.PWR[i-1] < htr_thr: ihtr.append(i-1)
    if data2.PWR[i] < htr_thr and data2.PWR[i-1] >= htr_thr: ihtr.append(i-1)


  def fitfunc(t, t0, T1,T2,tau):
    ret = numpy.zeros_like(t)
    ret[t<t0]  = T1
    ret[t>=t0] = T2 + (T1-T2)*numpy.exp(-(t[t>=t0]-t0)/tau)
    return ret

  (fig, ax) = plt.subplots(2,2)
  ttcx=[]
  ttcy=[]
  taux=[]
  tauy=[]
  taue=[]
  for i in ihtr:
    t0 = data2.TT[i]
    ii1 = numpy.logical_and(data1.TT >= t0 - dtm, data1.TT <= t0 + dtp)
    ii2 = numpy.logical_and(data2.TT >= t0 - dtm, data2.TT <= t0 + dtp)
    time1 = data1.TT[ii1] - t0
    time2 = data2.TT[ii2] - t0
    ttc=data1.ttc[ii1]
    pwr=data2.PWR[ii2] * 1e12 # power, pW
    pwr1 = numpy.mean(pwr[time2<-2])
    pwr2 = numpy.mean(pwr[time2>2])

    # choose only heater down
    if data2.PWR[i] < data2.PWR[i+1]: continue

    ax[0,0].plot(time1, ttc, '-')
    ax[1,0].plot(time2, pwr, '*-')

    # fit data
    par=[1, ttc[0], ttc[-1], 1]
    res = scipy.optimize.curve_fit(
       fitfunc, time1, ttc, par)
    par = res[0]
    err = numpy.sqrt(numpy.diag(res[1]))

    ttcx.append(pwr1)
    ttcx.append(pwr2)
    ttcy.append(par[1])
    ttcy.append(par[2])
    taux.append(pwr1)
    tauy.append(par[3])
    taue.append(err[3])

    # plot fit
    tt=numpy.linspace(-dtm, +dtp, 200)
    ax[0,0].plot(tt, fitfunc(tt,*par), 'k-', linewidth=0.7)
    ax[0,0].text(par[0]+2, par[1], 'tau: %.3f s\nT/Tc: %.4f %.4f\npwr: %.3f pW'%(par[3], par[1], par[2], pwr1))
#  ax[1,0].text(par[0]+1, pwr1,   'pwr: %.3f pW'%(pwr1*1e12))

  ax[0,1].plot(ttcx, ttcy, 'b*')
  pp=numpy.polyfit(ttcx, ttcy, 3)
  xx=numpy.linspace(0, numpy.max(ttcx) + 1)
  ax[0,1].plot(xx, numpy.polyval(pp,xx), 'k-')

  ax[1,1].errorbar(taux, tauy, taue, fmt='b*')

  ax[0,0].set_xlabel('time, s')
  ax[0,0].set_ylabel('T/Tc')
  ax[1,0].set_xlabel('time, s')
  ax[1,0].set_ylabel('power, pW')
  ax[0,1].set_xlabel('power, pW')
  ax[0,1].set_ylabel('T/Tc')
  plt.gcf().set_size_inches(12, 12)
  plt.savefig(fname + '.png', dpi=100)
  plt.close()

######################################################################
