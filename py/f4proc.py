# random functions for data processing

import numpy
import scipy.optimize
import matplotlib.pyplot as plt

import os
import re
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

