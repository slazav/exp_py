import numpy
import math
import scipy.optimize

###############################################################
# complex function used for fitting
def fitfunc(par,coord,F):

  V = (par[2] + 1j*par[3])/(par[4]**2 - F**2 + 1j*F*par[5])

  if not coord: V *= 1j*F  # velocity fit
  V += par[0] + 1j*par[1]

  if len(par)==8: V += (par[6] + 1j*par[7])*(F-par[4])
  return V

# function for minimization
def minfunc(par, coord,F,X,Y):
  V = fitfunc(par, coord,F)
  return numpy.linalg.norm(X + 1j*Y - V)

###############################################################
class fit_res_t:

  def __init__(self, time,drive,e,par,err, npars, coord):
    if len(par)!=8 or len(err)!=8:
      raise Exception("ERROR: fit_res_t: list of size 8 expected")
    self.time=time   # Mean time
    self.drive=drive # Mean drive
    self.e=e         # RMS difference between data and model
    self.coord=coord # coord/velocity fit
    self.npars=npars # 6/8 pars
    self.par=par     # parameters (8 members)
    self.err=err     # parameter uncertainties (8 members)
    self.A=par[0]
    self.B=par[1]
    self.C=par[2]
    self.D=par[3]
    self.f0=par[4]
    self.df=par[5]
    self.E=par[6]
    self.F=par[7]
    self.A_e=err[0]
    self.B_e=err[1]
    self.C_e=err[2]
    self.D_e=err[3]
    self.f0_e=err[4]
    self.df_e=err[5]
    self.E_e=err[6]
    self.F_e=err[7]
    self.amp=numpy.hypot(par[2], par[3])/par[5]
    if coord: self.amp/=par[4]

  # function
  def func(self, f):
    return fitfunc(self.par, self.coord, f)

  # Data array as in the database format (19 columns)
  # time, drive, e, A, Ae, B, Be, ...
  def dbfmt(self):
    return [self.time, self.drive, self.e,
      self.par[0], self.err[0], self.par[1], self.err[1],
      self.par[2], self.err[2], self.par[3], self.err[3],
      self.par[4], self.err[4], self.par[5], self.err[5],
      self.par[6], self.err[6], self.par[7], self.err[7]]

###############################################################
# Fit frequency sweep.
# Similar to fit_res program (https://github.com/slazav/fit_res)
# data is Nx5 numpy array with usual columns:
#   T-F-X-Y-D
# Return is similar to fit databases (19 columns):
#   time, drive, err, A, Aerr, B, Berr, C, Cerr,
#   D, Derr, f0, f0err, df, dferr, E, Eerr, F, Ferr

def fit(data, coord=0, npars=6, do_fit=1):

  if npars!=6 and npars!=8:
    raise Exception("ERROR: npars should be 6 or 8")

  k = numpy.max(numpy.hypot(data[:,2], data[:,3]))
  FF = data[:,1]
  XX = data[:,2]/k
  YY = data[:,3]/k
  time  = numpy.mean(data[:,0])
  drive = numpy.mean(data[:,4])

  ##########################
  # Initial conditions.
  # points with min/max freq
  ifmin = numpy.argmin(FF)
  ifmax = numpy.argmax(FF)

  # A,B - in the middle between these points:
  A = (XX[ifmin] + XX[ifmax])/2
  B = (YY[ifmin] + YY[ifmax])/2

  # furthest point from (A,B),
  # it should be near resonance:
  ires = numpy.argmax(numpy.hypot(XX-A, YY-B))
  F0 = FF[ires]

  # min/max freq where distance > dmax/sqrt(2),
  # this is resonance width:
  ii = numpy.hypot(XX-A, YY-B) > numpy.hypot(XX[ires]-A, YY[ires]-B)/math.sqrt(2)
  dF = numpy.max(FF[ii]) - numpy.min(FF[ii])
  if dF == 0:
    dF = abs(FF[max(0,ires-1)] - FF[min(ires+1,FF.size-1)])

  # amplitudes:
  if coord:
    C = -F0*dF*(YY[ires]-B);
    D =  F0*dF*(XX[ires]-A);
  else:
    C = dF*(XX[ires]-A);
    D = dF*(YY[ires]-B);

  # E,F - slope of the line connecting first and line points
  E = (XX[-1] - XX[0])/(FF[-1] - FF[0]);
  F = (YY[-1] - YY[0])/(FF[-1] - FF[0]);

  ##########################

  par = [A,B,C,D,F0,dF]
  err = [0,0,0,0,0,0]
  e   = 0
  if npars==8:
    par.extend((E,F))
    err.extend((0,0))

  if do_fit:
    res = scipy.optimize.minimize(minfunc, par, (coord, FF,XX,YY),
      options={'disp': False, 'maxiter': 1000})

    # Parameter uncertainty which corresponds to res.fun
    # which is relative RMS difference between function and the model.
    # df = d2f/dx2 dx2 -> dx = dqrt(0.5*df*H^-1)
    err = numpy.sqrt(0.5*res.fun*numpy.diag(res.hess_inv)).tolist()
    par = res.x.tolist()
    e = res.fun/FF.size

  if len(par)<8: par = par + [0,]*(8-len(par))
  if len(err)<8: err = err + [0,]*(8-len(err))

  e*=k
  for i in (0,1,2,3,6,7):
    par[i]*=k
    err[i]*=k

  ret = fit_res_t(time, drive, e, par, err, npars, coord)

  return ret

