import numpy
import math
import scipy.optimize
import fit_res_common

# Changes between fit_res001 and fit_res002:
# - Use drive in fitting, parameters A,B,C,D,E,F should be scaled with drive.
#   Data with different drives can be fitted together.
# - Remove fit_res_t.dbfmt
# - Non-linear function for low-temperature B-phase
#
# Changes between fit_res002 and fit_res003:
# - Use arbitrary "D-function" (delta(delta0, |vel|)) instead
#   of B-phase corrections.
#   No knowledge about field, pressure, wire parameters.

###############################################################
# Complex function used for fitting
# Loerentzian with S-function
def fitfuncD(par, dfunc, kv, F,D):
  if par[2]==0 and par[3]==0:
    return numpy.zeros_like(F)
  if par[5]==0:
    return numpy.zeros_like(F)
  V0 = 0
  E0 = 2
  while 1:
    delta = dfunc(par[5], abs(V0)*kv)
    V = 1j*F * (par[2] + 1j*par[3])*D/(par[4]**2 - F**2 + 1j*F*delta)
    V[numpy.isnan(V)] = 0
    V[numpy.isinf(V)] = 0
    E = numpy.max(abs(V-V0)/abs(V))
    if E < 1e-6: break
    if E > E0: break # avoid infinite growth of V
    V0 = V
    E0 = E

  V += (par[0] + 1j*par[1])*D
  if len(par)==8: V += (par[6] + 1j*par[7])*(F-par[4])*D
  return V

# function for minimization
def minfuncD(par, dfunc, kv, F,X,Y,D):
  V = fitfuncD(par, dfunc, kv, F,D)
  return numpy.linalg.norm(X + 1j*Y - V)/numpy.linalg.norm(X + 1j*Y)


###############################################################
class fit_res_t:
  def __init__(self, time,e,par,err, npars, coord, dfunc):
    if len(par)!=8 or len(err)!=8:
      raise Exception("ERROR: fit_res_t: list of size 8 expected")

    if coord!=0 and dfunc!=None:
      raise Exception("ERROR: in dfunc mode coord parameter should be 0")

    self.time=time   # Mean time
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
    self.dfunc=dfunc
    if par[5]!=0 and par[4]!=0:
      if coord:
        self.VX=-par[2]/par[5]/par[4]
        self.VY=par[3]/par[5]/par[4]
      else:
        self.VX=par[3]/par[5]
        self.VY=-par[2]/par[5]
      self.amp=numpy.hypot(self.VX,self.VY)
    else:
      self.VX=float('nan')
      self.VY=float('nan')
      self.amp=float('nan')

  # function
  def func(self, f,d):
    if self.dfunc!=None:
      return fitfuncD(self.par, self.dfunc, 1, f,d)
    else:
      return fit_res_common.fitfunc_lor(self.par, self.coord, f,d)

  def func_lin(self, f,d):
    return fit_res_common.fitfunc_lor(self.par, self.coord, f,d)


###############################################################
# Fit frequency sweeps.
# Similar to fit_res program (https://github.com/slazav/fit_res)
# data is Nx5 numpy array with usual columns:
#   T-F-X-Y-D

def fit(data, coord=0, npars=6, dfunc=None, do_fit=1):

  if npars!=6 and npars!=8:
    raise Exception("ERROR: npars should be 6 or 8")

  if coord!=0 and dfunc!=None:
    raise Exception("ERROR: in dfunc mode coord parameter should be 0")

  kv = numpy.max(numpy.hypot(data[:,2], data[:,3]))
  kd = numpy.max(abs(data[:,4]))
  time  = numpy.mean(data[:,0])
  FF = data[:,1]
  XX = data[:,2]/kv
  YY = data[:,3]/kv
  DD = data[:,4]/kd

  (A,B,C,D,F0,dF,E,F) = fit_res_common.init_pars(FF, XX/DD, YY/DD, coord)
  par = [A,B,C,D,F0,dF]
  err = [0,0,0,0,0,0]
  e   = 0
  if npars==8:
    par.extend((E,F))
    err.extend((0,0))

  if do_fit:
    if dfunc!=None:
      res = scipy.optimize.minimize(minfuncD, par, (dfunc, kv, FF,XX,YY,DD),
        options={'disp': False, 'maxiter': 1000})
    else:
      res = scipy.optimize.minimize(fit_res_common.fitfunc_lor, par, (coord, FF,XX,YY,DD),
        options={'disp': False, 'maxiter': 1000})

    # Parameter uncertainty which corresponds to res.fun
    # which is relative RMS difference between function and the model.
    # df = d2f/dx2 dx2 -> dx = dqrt(0.5*df*H^-1)
    err = numpy.sqrt(0.5*res.fun*numpy.diag(res.hess_inv)).tolist()
    par = res.x.tolist()
    e = res.fun/FF.size

  if len(par)<8: par = par + [0,]*(8-len(par))
  if len(err)<8: err = err + [0,]*(8-len(err))

  e*=kv
  for i in (0,1,2,3,6,7):
    par[i]*=kv/kd
    err[i]*=kv/kd

  return fit_res_t(time, e, par, err, npars, coord, dfunc)

###############################################################
# Plot fit

def plot(ax,ay, sweep, fit, **kargs):
  fit_res_common.plot(ax,ay, sweep, fit, **kargs)
