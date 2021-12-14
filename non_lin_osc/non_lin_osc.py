#!/usr/bin/python3

import numpy
import math
import scipy.optimize
import scipy.integrate

#######################################################

# Equilibriun function of a non-linear oscillator
# at periodic drive.
#
# u,v - van der Pol coordinates (x,x' rotated by w*t)
# w - frequency
# func - second derivative function, x'' = -func(x, x')
#
# Coordinated are averaged over rotation period (phase 0..2*pi)
# Output is u' and v',  they should be zero in equilibrium.
#
# See text in http://www.scholarpedia.org/article/Duffing_oscillator
#

#######################################################
# Find (\dot u,\dot v) as a function of (u,v) for period-averaged motion of
# a periodically-driven non-linear oscillator x'' + func(x,x') = F*cos(w*t)
def nonlin_osc_eq(uv,w, F, func, fpars):
  p = 2*math.pi*numpy.linspace(0,1,100) # phase for integration
  sp = numpy.sin(p)
  cp = numpy.cos(p)

  x = uv[0]*cp-uv[1]*sp;
  y = w*(-uv[0]*sp-uv[1]*cp);
  z = (F*cp - func(x,y,fpars) + w**2*x)/w;

  return [-scipy.integrate.trapz(z*sp, p),\
          -scipy.integrate.trapz(z*cp, p)];

#######################################################

# Find equilibrium (zero of nonlin_osc_eq function).
# This is enough for simple small non-linearities,
# for duffing oscillator it does not work properly
# (one should integrate trajectories in u-v space instead)

def nonlin_osc_solve0(w, F, func, fpars, uv0):
  return scipy.optimize.fsolve(nonlin_osc_eq, uv0, (w,F,func,fpars))

#######################################################

## harmonic osc
def osc_harm(x,dx, pars):
  w0  = pars[0] # resonant frequency
  tau = pars[1] # relaxation time at low velocuty
  return w0**2*x + dx/tau

## pseudoplastic osc N1
def osc_pseudopl1(x,dx, pars):
  w0  = pars[0] # resonant frequency
  tau = pars[1] # relaxation time at low velocuty
  vc  = pars[2] # critical velocity
  k   = pars[3] # relaxation factor for large velocity
  return w0**2*x + dx/tau * (k - (1-k) * vc/numpy.sqrt(vc**2 + dx**2))

## pseudoplastic osc N2
def osc_pseudopl2(x,dx, pars):
  w0  = pars[0] # resonant frequency
  tau = pars[1] # relaxation time at low velocuty
  vc  = pars[2] # critical velocity
  return w0**2*x + dx/tau * vc/numpy.sqrt(vc**2 + dx**2)

## duffing osc
def osc_duffing(x,dx, pars):
  w0  = pars[0] # resonant frequency
  tau = pars[1] # relaxation time at low velocuty
  a   = pars[2] # non-linear parameter
  return w0**2*x + dx/tau + a*x**3

#######################################################
# Equilibriun function and solver for Duffung oscillator at periodic drive.
# see http://www.scholarpedia.org/article/Duffing_oscillator
# Note that solver does not work when hysteresis appears.

def nonlin_osc_eq_duff(uv, w, F, w0, tau, a):
  return [(-w*uv[0]/tau + (w**2-w0**2)*uv[1] - 3/4.0*a*(uv[1]**2+uv[0]**2)*uv[1])/2/w,\
          (-w*uv[1]/tau - (w**2-w0**2)*uv[0] + 3/4.0*a*(uv[1]**2+uv[0]**2)*uv[0] - F)/2/w]

def nonlin_osc_solve_duff(w, F, w0, tau, a, uv0):
  return scipy.optimize.fsolve(nonlin_osc_eq_duff, uv0, (w,F,w0, tau, a))
