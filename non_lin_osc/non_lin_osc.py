#!/usr/bin/python3

import numpy
import math
import scipy.optimize
from scipy.integrate import solve_bvp,trapz,quad

#######################################################
#######################################################
# Here we are working with a non-linear oscillator
# driven by periodic force:
#     x'' + func(x, x', pars) = F*cos(w*t)
#
# Examples of oscillatiors:

## harmonic osc
def osc_harm(x,dx, pars):
  (w0,tau)  = pars
  return w0**2*x + dx/tau

## duffing osc
def osc_duffing(x,dx, pars):
  (w0,tau, a) = pars
  return w0**2*x + dx/tau + a*x**3

## pseudoplastic osc N1
def osc_pseudopl1(x,dx, pars):
  (w0,tau,vc,k) = pars
  return w0**2*x + dx/tau * (k - (1-k) * vc/numpy.sqrt(vc**2 + dx**2))

## pseudoplastic osc N2
def osc_pseudopl2(x,dx, pars):
  (w0,tau,vc) = pars
  return w0**2*x + dx/tau * vc/numpy.sqrt(vc**2 + dx**2)


############################
# Solve ODE of the non-linear oscillator driven by periodic force
# on one period (2pi/w) with periodic BC, x(0)=x(T), x'(0)=x'(T).
# Returns full solution as a cubic spline (see solve_bvp() return value).
# Arguments:
#   func - function
#   pars - parameters for the function
#   F - drive amplitude
#   w - drive frequency
# Return:
#   result from solve_bvp()
#
def osc_solve_per_func(func, pars, F, w, a0=0, p0=0, nper=1):

  def bc(Xa, Xb):
    return numpy.array([Xa[0]-Xb[0], Xa[1]-Xb[1]])

  # Two equations: x' = v; v' = RHS(x,v)
  def rhs_func(t, X):
    RHS = -func(X[0], X[1], pars) - F*numpy.cos(w*t)
    return numpy.vstack((X[1], RHS))

  # initial mesh:
  t0 = numpy.linspace(0, nper*2*numpy.pi/w, 50)

  # initial guess
  x0 = numpy.zeros((2, t0.size))
  if a0 != 0:
    x0[0,:] =  a0*numpy.cos(w*t0+p0)
    x0[1,:] = -a0*w*numpy.sin(w*t0+p0)
  res = solve_bvp(rhs_func, bc, t0, x0, tol=1e-4, max_nodes=10000)
  res.nper = nper
  res.w = w

  return res

############################
# Calculate N-th harmonic of a function returned by
# osc_solve_per()
#
def osc_solve_per_harm(res, N):
  if res.status>0: return (numpy.nan,numpy.nan)
  T=res.x[-1]
  W=res.w*N

#  fx = numpy.cos(W*res.x)*res.y[0]
#  fy = numpy.sin(W*res.x)*res.y[0]
#  print('2>', res.x.size, fy[0], fy[-1])
#  x = trapz(fx, res.x)
#  y = trapz(fy, res.x)

  fx = lambda t: numpy.cos(W*t)*res.sol(t)[0]
  fy = lambda t: numpy.sin(W*t)*res.sol(t)[0]
  x = quad(fx, 0, T, limit=200)[0]
  y = quad(fy, 0, T, limit=200)[0]
  return (2*x/T, 2*y/T)


############################
# Calculate response of a non-linear oscillator on the driving frequency.
# Combination of osc_solve_per_func() and osc_solve_per_harm() calls.
#
def osc_solve_per(func, pars, F, w, a0=0, p0=0):
  res = osc_solve_per_func(func, pars, F, w, a0, p0)
  return osc_solve_per_harm(res,1)


#######################################################
#######################################################

# Equilibriun function of a non-linear oscillator
# at periodic drive.
#
# u,v - van der Pol coordinates (x,x'/w rotated by w*t)
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
#
def osc_solve_vdp_eq(uv, func, pars, F, w):
  p = 2*math.pi*numpy.linspace(0,1,100) # phase for integration
  sp = numpy.sin(p)
  cp = numpy.cos(p)

  x = uv[0]*cp-uv[1]*sp;
  y = w*(-uv[0]*sp-uv[1]*cp);
  z = (F*cp - func(x,y,pars) + w**2*x)/w;

  return [-trapz(z*sp, p), -trapz(z*cp, p)];

#######################################################

# Find equilibrium (zero of nonlin_osc_eq function).
# This is enough for simple small non-linearities,
# for duffing oscillator it does not work properly
# (one should integrate trajectories in u-v space instead)
#
def osc_solve_vdp(func, pars, F, w, a0=0, p0=0):
  uv0 = [a0*numpy.cos(p0), a0*numpy.sin(p0)]
  (uv0, info, stat, msg) = scipy.optimize.fsolve(osc_solve_vdp_eq, uv0, (func,pars,F,w,), full_output=True)
  if stat!=1: return (numpy.nan, numpy.nan)

  # Equilibrium is unstable if Jacobian
  # has an eigenvalue with positive real part.
  # or if (trJ<0 || detJ>0)
  # see https://stackoverflow.com/questions/52309635/check-the-stability-of-differential-equation-using-fsolve
  #     http://www.scholarpedia.org/article/Equilibrium
  # does not work yet?
#  R = numpy.zeros((2,2))
#  R[numpy.triu_indices(2)] = info["r"]
#  J = (info["fjac"].T).dot(R)
#  ee = numpy.real(numpy.linalg.eigvals(J))
#  if ee[0]>0 or ee[1]>0: return (numpy.nan, numpy.nan)

  return uv0

#######################################################
#######################################################
# Equilibriun function and solver for Duffung oscillator at periodic drive.
# see http://www.scholarpedia.org/article/Duffing_oscillator
# Note that solver does not work when hysteresis appears.

def osc_solve_vdp_duff_eq(uv, w0, tau, a, F, w):
  return [(-w*uv[0]/tau + (w**2-w0**2)*uv[1] - 3/4.0*a*(uv[1]**2+uv[0]**2)*uv[1])/2/w,\
          (-w*uv[1]/tau - (w**2-w0**2)*uv[0] + 3/4.0*a*(uv[1]**2+uv[0]**2)*uv[0] - F)/2/w]

def osc_solve_vdp_duff(w0, tau, a, F, w, a0=0, p0=0):
  uv0 = [a0*numpy.cos(p0), a0*numpy.sin(p0)]
  (uv0, info, stat, msg) = scipy.optimize.fsolve(osc_solve_vdp_duff_eq, uv0, (w0, tau, a, F, w), full_output=True)
  if stat!=1: return (numpy.nan, numpy.nan)
  return uv0

# Analytical expression for amplitude of Duffing oscillator.
#  ((w**2 - w0**2 - 3/4.0*a * A**2)**2 + (w/tau)**2)*A**2 = F**2
# This is a qubic equation for A**2.
# Returns up to 3 numbers
#
def osc_duff_amp(w,F,w0,tau, a):
  p = [9/16.0*a**2,  - 3/2.0*(w**2-w0**2)*a, (w**2-w0**2)**2 + (w/tau)**2, -F**2]
  rr = numpy.sqrt(numpy.roots(p))
  return numpy.real(rr[numpy.imag(rr)==0])
