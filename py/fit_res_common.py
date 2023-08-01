import numpy
import math

###############################################################
# Initial conditions.
def init_pars(FF,XS,YS,coord):

  # points with min/max freq (note that frequency could be non-monotonic)
  ifmin = numpy.argmin(FF)
  ifmax = numpy.argmax(FF)

  xmin=XS[ifmin]; xmax=XS[ifmax]
  ymin=YS[ifmin]; ymax=YS[ifmax]
  fmin=FF[ifmin]; fmax=FF[ifmax]

  # A,B - in the middle between these points:
  A = (xmin+xmax)/2
  B = (ymin+ymax)/2

  # E,F - slope of the line connecting first and line points
  E = (xmax-xmin)/(fmax-fmin);
  F = (ymax-ymin)/(fmax-fmin);

  # furthest point from line connecting min and max point
  # it should be near resonance:
  dist = numpy.hypot(XS - xmin - (FF-fmin)*E, YS - ymin - (FF-fmin)*F)
  ires = numpy.argmax(dist)
  F0 = FF[ires]

  # min/max freq where distance > dmax/sqrt(2),
  # this is resonance width:
  ii = dist > dist[ires]/math.sqrt(2)
  dF = numpy.max(FF[ii]) - numpy.min(FF[ii])
  if dF == 0:
    dF = abs(FF[max(0,ires-1)] - FF[min(ires+1,FF.size-1)])

  # amplitudes:
  if coord:
    C = -F0*dF*(YS[ires]-B);
    D =  F0*dF*(XS[ires]-A);
  else:
    C = dF*(XS[ires]-A);
    D = dF*(YS[ires]-B);

  return (A,B,C,D,F0,dF,E,F)

###############################################################
# Normal lorentzian function:

# complex function used for fitting
def fitfunc_lor(par,coord,F,D=1):
  V = (par[2] + 1j*par[3])*D/(par[4]**2 - F**2 + 1j*F*par[5])
  if not coord: V *= 1j*F  # velocity fit
  V += (par[0] + 1j*par[1])*D
  if len(par)==8: V += (par[6] + 1j*par[7])*(F-par[4])*D
  return V

# function for minimization
def minfunc_lor(par, coord,F,X,Y,D=1):
  V = fitfunc_lor(par, coord,F,D)
  print(par)
  return numpy.linalg.norm(X + 1j*Y - V)

  #/numpy.linalg.norm(X + 1j*Y)

###############################################################
# Plot fit.
# Parameters:
# - ax,ay   - axis for X and Y (could be same)
# - sweep   - sweep data
# - fit_res - fit result
# - npts    - number of points for fit plot
# - sh      - y-shift
# - sc      - y-scale

def plot(ax,ay, sweep, fit, npts=100, sh=0, sc=1, xlabel=None, ylabel=None, plot_bg=1, plot_lin=1):
  import matplotlib.pyplot as plt

  # sweep data
  ax.plot(sweep[:,1], sh+sc*sweep[:,2], 'r.', label=xlabel)
  ay.plot(sweep[:,1], sh+sc*sweep[:,3], 'b.', label=ylabel)

  # note that drive could be non-constant!
  drive=numpy.mean(sweep[:,4])
  ff=numpy.linspace(min(sweep[:,1]), max(sweep[:,1]), npts)

  # fit result
  if fit != None:
    vv = sh*(1 + 1j) + sc*fit.func(ff, drive)
    ax.plot(ff, numpy.real(vv), 'k-', linewidth=1)
    ay.plot(ff, numpy.imag(vv), 'k-', linewidth=1)

  # plot also linear Lorentzian
  if plot_lin and fit.func_lin:
    vv = sh*(1 + 1j) + sc*fit.func_lin(ff, drive)
    ax.plot(ff, numpy.real(vv), 'k--', linewidth=0.7)
    ay.plot(ff, numpy.imag(vv), 'k--', linewidth=0.7)

  # plot background
  if plot_bg:
    vv = sh*(1 + 1j) + sc*drive*(fit.par[0]+1j*fit.par[1] + (fit.par[-2]+1j*fit.par[-1])*(ff-fit.par[4]))
    ax.plot(ff, numpy.real(vv), 'k--', linewidth=0.7)
    ay.plot(ff, numpy.imag(vv), 'k--', linewidth=0.7)

  ax.set_xlabel("freq, Hz")
  if ax!=ay:
    ay.set_xlabel("freq, Hz")

