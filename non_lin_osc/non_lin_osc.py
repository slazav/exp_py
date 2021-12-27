#!/usr/bin/python3

import numpy
import math
import scipy.optimize
from scipy.integrate import solve_bvp,trapz,quad

#######################################################
#######################################################

############################
# Solve ODE of a non-linear oscillator driven by periodic force:
#   x'' + func(x, x', pars) = F*cos(w*t)
# on one period (2pi/w) with periodic BC, x(0)=x(T), x'(0)=x'(T).
# Returns full solution as a cubic spline (see solve_bvp()).
# Arguments:
#   func - function
#   pars - parameters for the function
#   F - drive amplitude
#   w - drive frequency
# Return:
#   result from solve_bvp()
#
def osc_solve_per_func(func, pars, F, w):

  def bc(Xa, Xb):
    return numpy.array([Xa[0]-Xb[0], Xa[1]-Xb[1]])

  def rhs_func(t, X):
    RHS = -func(X[0], X[1], pars) - F*numpy.cos(w*t)
    return numpy.vstack((X[1], RHS))

  # initial mesh:
  t0 = numpy.linspace(0, 2*math.pi/w, 5)
  x0 = numpy.zeros((2, t0.size))
  return solve_bvp(rhs_func, bc, t0, x0)

############################
# Calculate N-th harmonic of a function returned by
# osc_solve_per()
#
def osc_solve_per_harm(res, N):
  T=res.x[-1]
  W=2*math.pi/T * N
  fx = lambda t: numpy.cos(W*t)*res.sol(t)[0]
  fy = lambda t: numpy.sin(W*t)*res.sol(t)[0]
  x = quad(fx, 0, T, limit=200)[0]
  y = quad(fy, 0, T, limit=200)[0]
  return (2*x/T, 2*y/T)

############################
# Calculate response of a non-linear oscillator on the driving frequency.
# Combination of osc_solve_per_func() and osc_solve_per_harm() calls.
#
def osc_solve_per(func, pars, F, w):
  res = osc_solve_per_func(func, pars, F, w)
  return osc_solve_per_harm(res,1)


#######################################################
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

  return [-trapz(z*sp, p),\
          -trapz(z*cp, p)];

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
