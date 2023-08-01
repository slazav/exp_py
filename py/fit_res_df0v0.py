import numpy
import math
import scipy.optimize
import fit_res_common

# Version of fit_res with additional parameter v0.
# A simple B-phase damping model with
# delta(delta0, v0) = delta0/S(v/v0)
# with v0 as a fitting parameter.
# 
# V.Zavjalov, 15.06.2023

# Fit of theoretical S function
A0=0.486387040089067
B0=1.14243469051564

###############################################################
# Complex function used for fitting
# (A+i*B) + i*f*(C+i*D)/()
def fitfunc(par, F,D):
  AB=(par[0] + 1j*par[1])
  CD=(par[2] + 1j*par[3])
  F0  = par[4]
  dF0 = par[5]
  V0  = par[6]

  if par[2]==0 and par[3]==0:
    return numpy.zeros_like(F)
  if par[5]==0:
    return numpy.zeros_like(F)

  Vx = 0
  E0 = 2
  while 1:
    VV=abs(Vx)/V0
    dF = dF0/(1 + A0*abs(VV)**B0)

    V = 1j*F * CD*D/(F0**2-F**2 + 1j*F*dF)
    V[numpy.isnan(V)] = 0
    V[numpy.isinf(V)] = 0
    E = numpy.max(abs(V-Vx)/abs(V))
    if E < 1e-6: break
    if E > E0: break # avoid infinite growth of V
    Vx = V
    Ex = E

  V += AB*D
  if len(par)==9: V += (par[7] + 1j*par[8])*(F-F0)*D
  return V

# function for minimization
def minfunc(par, F,X,Y,D):
  V = fitfunc(par, F,D)
  return numpy.linalg.norm(X + 1j*Y - V)/numpy.linalg.norm(X + 1j*Y)

###############################################################
class fit_res_t:
  def __init__(self, time,e,par,err, npars):
    if len(par)!=9 or len(err)!=9:
      raise Exception("ERROR: fit_res_t: list of size 9 expected")

    self.time=time   # Mean time
    self.e=e         # RMS difference between data and model
    self.npars=npars # 7/9 pars
    self.coord=0     #
    self.par=par     # parameters (9 members)
    self.err=err     # parameter uncertainties (9 members)
    self.A=par[0]
    self.B=par[1]
    self.C=par[2]
    self.D=par[3]
    self.f0=par[4]
    self.df=par[5]
    self.v0=par[6]
    self.E=par[7]
    self.F=par[8]
    self.A_e=err[0]
    self.B_e=err[1]
    self.C_e=err[2]
    self.D_e=err[3]
    self.f0_e=err[4]
    self.df_e=err[5]
    self.v0_e=err[6]
    self.E_e=err[7]
    self.F_e=err[8]
    if par[5]!=0 and par[4]!=0:
      self.VX=par[3]/par[5]
      self.VY=-par[2]/par[5]
      self.amp=numpy.hypot(self.VX,self.VY)
    else:
      self.VX=float('nan')
      self.VY=float('nan')
      self.amp=float('nan')

  # function
  def func(self, f,d):
    return fitfunc(self.par, f,d)

  def func_lin(self, f,d):
    return fit_res_common.fitfunc_lor(
      numpy.take(self.par, (0,1,2,3,4,5,7,8)), self.coord, f,d)

  def dfunc(self, df0, v):
    VV=v/self.v0
    return df0/(1 + A0*abs(VV)**B0)

###############################################################
# Fit frequency sweeps.
# Similar to fit_res program (https://github.com/slazav/fit_res)
# data is Nx5 numpy array with usual columns:
#   T-F-X-Y-D

def fit(data, npars=7, do_fit=1):

  if npars!=7 and npars!=9:
    raise Exception("ERROR: npars should be 7 or 9")

  kv = numpy.max(numpy.hypot(data[:,2], data[:,3]))
  kd = numpy.max(abs(data[:,4]))
  time  = numpy.mean(data[:,0])
  FF = data[:,1]
  XX = data[:,2]/kv
  YY = data[:,3]/kv
  DD = data[:,4]/kd

  (A,B,C,D,F0,dF,E,F) = fit_res_common.init_pars(FF, XX/DD, YY/DD, 0)
  V0 = numpy.hypot(C,D)/dF

  par = [A,B,C,D,F0,dF,V0]
  err = [0,0,0,0,0,0,0]
  e   = 0
  if npars==9:
    par.extend((E,F))
    err.extend((0,0))

  if do_fit:
    res = scipy.optimize.minimize(minfunc, par, (FF,XX,YY,DD),
        options={'disp': False, 'maxiter': 10000, 'eps': 1e-9})

    # Parameter uncertainty which corresponds to res.fun
    # which is relative RMS difference between function and the model.
    # df = d2f/dx2 dx2 -> dx = dqrt(0.5*df*H^-1)
    err = numpy.sqrt(0.5*res.fun*numpy.diag(res.hess_inv)).tolist()
    par = res.x.tolist()
    e = res.fun/FF.size

  if len(par)<9: par = par + [0,]*(9-len(par))
  if len(err)<9: err = err + [0,]*(9-len(err))

  e*=kv
  for i in (0,1,2,3,7,8):
    par[i]*=kv/kd
    err[i]*=kv/kd

  for i in (6,):
    par[i]*=kv
    err[i]*=kv

  return fit_res_t(time, e, par, err, npars)

###############################################################
# Plot fit.

def plot(ax,ay, sweep, fit, **kargs):
  fit_res_common.plot(ax,ay, sweep, fit, **kargs)
