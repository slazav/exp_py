# functions for finding signal frequency

import numpy
from math import pi,hypot,ceil,floor

###################################

# fft maximum
def find_freq_fft(T,A, fmin=-1, fmax=-1):
  if T.size!=A.size:
    raise ValueError("T and A arrays have different number of points")

  N = T.size - 1
  tstep = (T[-1]-T[0])/(T.size-1)
  df = 1/tstep/(N+1)

  imin=0
  imax=ceil(N/2)
  if fmin>0: imin = floor(fmin/df)
  if fmax>0: imax = ceil(fmax/df)

  if (imax <= imin):
    raise ValueError("Wrong fmin/fmax setting: %e %e" %(fmin, fmax))

  fft = numpy.fft.fft(A)
  fft[0] = 0  # remove constant component
  f_ind = numpy.argmax(numpy.abs(fft[imin:imax])) + imin

  return f_ind*df, fft, f_ind

###################################

# fit fft maximum by 3-point parabola
def find_freq_fftfit_q(T,A, fmin=-1, fmax=-1):
  freq, fft, f_ind = find_freq_fft(T,A, fmin, fmax)
  v1 = abs(fft[f_ind-1])
  v2 = abs(fft[f_ind])
  v3 = abs(fft[f_ind+1])
  di = -(v3-v1)/(v1 - 2*v2 + v3)/2
  return freq*(1 + di/f_ind)

###################################

# Fit fft near maximum.
# This is similar method which I use for fitting
# decaying signals in
# https://github.com/slazav/pico_osc/tree/master/modules/fit_signal

def find_freq_fftfit_l(T,A, fmin=-1, fmax=-1):

  freq1, fft1, f_ind1 = find_freq_fft(T,A, fmin, fmax)

  # Now we should solve problem with too "pure"
  # signals where we have maxinun and zeros near it.

  # we want to remove k points to shift
  # fft node by 1/2 near the resonance i:
  #  i - i (n-k)/n = 1/2,
  #  k = 1/2 n/i,
  N = T.size
  k = int(N/2/f_ind1)

  freq2, fft2, f_ind2 = find_freq_fft(T[0:-k],A[0:-k], fmin, fmax)

  # choose signal with smaller maximum
  # larger side points:
  if abs(fft1[f_ind1]) > abs(fft2[f_ind2]):
    freq = freq2
    fft = fft2
    f_ind = f_ind2
  else:
    freq = freq1
    fft = fft1
    f_ind = f_ind1

  # fft is proportional to 1/(f-f0)
  # We want to fit 1/fft it with linear function A*x+B
  # and find zero crossing (-B/A)
  # As usual, we minimize  S = sum (a*x +b - y)^2
  # Differentiating by A and B we have
  # A*sum(x*x) + B*sum(x) - sum(y*x) =  0
  # A*sum(x) + B*sum(w) - sum(y) =  0

  # Solving this we have:
  # B = [sum(y*x)*sum(x) - sum(y)*sum(x*x)]/[sum(x)*sum(x) - sum(1)*sum(x*x)]
  # A = [sum(y)*sum(x) - sum(y*x)*sum(1)]/[sum(x)*sum(x) - sum(1)*sum(x*x)]
  # x0 = -B/A = (sum(y*x)*sum(x) - sum(y)*sum(x*x)) / (sum(y*x)*sum(1) - sum(y)*sum(x))

  # it's good to fit in a very narrow window, because 1/fft is very noisy far from f0
  ind_win = 2 #
  sn = 0
  sx = 0
  sy = 0
  sxx = 0
  sxy = 0
  for x in range(-ind_win, +ind_win+1):
    if fft[f_ind + x]==0: continue
    y = numpy.real(1/fft[f_ind + x])
    w = 1/y**2 # Very strong weighting function
    sn  += w
    sx  += x*w
    sy  += y*w
    sxx += x*x*w
    sxy += x*y*w

  di = (sxy*sx-sy*sxx)/(sxy*sn-sy*sx)
  # Another approach, much less accurate:
  # v1 = fft[f_ind-1]
  # v2 = fft[f_ind]
  # v3 = fft[f_ind+1]
  # f_indm = f_ind - (v3-v1)/(v1 - 2*v2 + v3)/2

  return freq*(1 + di/f_ind)

###################################
# calculate fourier component of the signal
# at frequency F and its harmonics

from scipy.integrate import trapz

def calc_fourier(freq, T,A, harm=[1]):

  if freq==0:
    raise ValueError("calc_fourier: freq=0")

  if 1/freq > T[-1]-T[0]:
    raise ValueError("calc_fourier: signal length is smaller then one period")

  # We want to do integration over integer number of periods.
  # Let's split it into integral over the signal plus some tail which
  # extends from signal end to the nearest end of period.
  # fill the tail with data using periodicity

  tspan = T[-1]-T[0]
  tend = T[0] + ceil(tspan*freq)/freq
  ii = numpy.logical_and(T+1/freq > T[-1], T+1/freq < tend)
  T1 = numpy.array([T[-1]] + list(T[ii]+1/freq) + [tend])
  A1 = numpy.array([A[-1]] + list(A[ii]) + [A[0]])

  ret = numpy.zeros(len(harm), dtype=complex)
  for i in range(len(harm)):
    if harm[i] == 0:
      ret[i] = numpy.mean(A)
      continue

    E  = numpy.exp(1j*harm[i]*2*pi*freq*T)
    E1 = numpy.exp(1j*harm[i]*2*pi*freq*T1)

    ret[i] = 2*(trapz(E*A, T) + trapz(E1*A1, T1))/(tend-T[0])

  return ret

###################################

# one more approach: we calculate "slow" Furier at some frequency F and 
# then maximize its value as a function of F

from scipy.optimize import minimize

def find_freq_fmax(T,A, fmin=-1, fmax=-1):

  # initial frequency value:
  freq, fft, f_ind = find_freq_fft(T,A, fmin, fmax)
  df = 1/(T[-1]-T[0]) # frequency resolution
  def minfunc(freq, T,A):
    return -numpy.abs(calc_fourier(freq[0], T,A))

  res = minimize(minfunc, freq, (T,A),
    bounds=[(freq-df,freq+df)],
    options={'disp': False, 'maxiter': 1000})
  return res.x[0]

