import numpy
import re
import graphene

# Get a few sweeps for wire <name> starting at time range <t1>...<t2>.
# If <t2>==None get a single sweep covering time <t1>
# Return list of dictionaries with sweep parameters.
#
# 1. Sweep parameters (including exact time span) are extracted
#    from <W>_pars database: t0, t1-t0, F1, F2, N, dt, dtf, dir, drive volt, drive ph
# 2. Fit extracted from <W> database:
#    Tcnt,drive volt,Err,A,Ae,B,Be,C,Ce,D,De,F,Fe,dF,dFe
# 3. F-sweep extracted from <W>_sweeps database:
#    T,F,X,Y

def get_sweep(name, t1):
  t1 = graphene.timeconv(t1)
  pars = graphene.get_prev(name + '_pars', t1, unpack=False, cols=None)
  # Note: special timestamps (inf, now, <t>+) will fail here!
  if pars[0] + pars[1] < float(t1): return (); # too old
  return get_sweeps_(name, pars)

def get_sweep_range(name, t1, t2):
  t1 = graphene.timeconv(t1)
  t2 = graphene.timeconv(t2)
  pars = graphene.get_range(name + '_pars', t1, t2, unpack=False, cols=None)
  return get_sweeps_(name, pars)

def get_sweep_list(name, tlist):
  pars = numpy.empty((0,10))
  for t in tlist:
    t = graphene.timeconv(t)
    p = graphene.get_prev(name + '_pars', t, unpack=False, cols=None)
    # Note: special timestamps (inf, now, <t>+) will fail here!
    if p[0] + p[1] < float(t): continue; # too old
    pars = numpy.append(pars, numpy.reshape(p,(1,10)), axis=0)
  return get_sweeps_(name, pars)


#######################
# Internal function for getting sweeps using
# sweep parameters from <name>_pars database.

def get_sweeps_(name, pars):
  sweeps=[]

  # no sweeps
  if pars.size == 0: return sweeps

  # single sweep
  if len(pars.shape) == 1:
    pars = numpy.reshape(pars, (1,-1))

  for i in range(pars.shape[0]):
    (t1p, dtp, f1, f2, n, s, s0, d, vpp, ph) = pars[i,:]
    t2p = t1p+dtp

    tcent = (t2p+t1p)/2
    fcent = (f2+f1)/2
    tspan = abs(t2p-t1p)
    fspan = abs(f2-f1)

    # Get sweep points
    (ts,fs,xs,ys) = graphene.get_range(name + '_sweeps', t1p, t2p, cols=None)

    # Get fit
    # (tcnt,dr,Err,A,Ae,B,Be,C,Ce,D,De,f0,f0e,df,dfe,[E,Ee,F,Fe])
    fit = graphene.get_range(name, t1p, t2p, cols=None, unpack=False).copy()
    if fit.size<19: fit.resize(19)

    amp = numpy.hypot(fit[7],fit[9])/fit[11]/fit[13]

    # get demag field
    Bdemag = graphene.get_prev('demag_pc:f2', tcent, cols=1)

    sweep = {
      't1': t1p, 't2': t2p,
      'tspan': tspan, 'tcent': tcent,
      'f1':  f1, 'f2':  f2, 'n': n,
      'fspan': fspan, 'fcent': fcent,
      'step': s, 'step0': s0, 'dir': d,
      'drive': vpp, 'phase': ph,
      'T': ts, 'F': fs, 'X': xs, 'Y': ys,
      'fit_err': fit[2],
      'fit_A': fit[3],   'fit_A_err': fit[4],
      'fit_B': fit[5],   'fit_B_err': fit[6],
      'fit_C': fit[7],   'fit_C_err': fit[8],
      'fit_D': fit[9],   'fit_D_err': fit[10],
      'fit_f0': fit[11], 'fit_f0_err': fit[12],
      'fit_df': fit[13], 'fit_df_err': fit[14],
      'fit_E': fit[15],  'fit_E_err': fit[16],
      'fit_F': fit[17],  'fit_F_err': fit[18],
      'fit_amp': amp,
      'demag': Bdemag
    }
    sweeps.append(sweep)

  return sweeps

