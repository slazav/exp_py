# B-phase correction with proc2 *.res2 files obtained with proc2 function,
# See data_f4_wire_b repository
# V.Zavjalov, 14.06.2023

import numpy

############################################################
# Convert vel+delta to delta0.
# (Extrapolation of velocity-dependent width to zero velocity).
# Useful for processing tracking mode.
#
def res2_delta_to_delta0(vel, delta, res2file):
  # name press field   A B s2 dfi   err(A, B, s2, dfi)
  (A,B,s2,dfi) = numpy.loadtxt(res2file, usecols=(3,4,5,6), unpack=1)

  A0=0.439023879201735
  B0=0.045824276758537
  C0=-0.00192214696297599

  # v0 = B/(A - numpy.log(abs(delta0)))
  # delta = delta0/(1 + A0*(v/v0) + s2*B0*(v/v0)**2) + dfi

  delta0=delta
  while 1:
    v0 = B/(A - numpy.log(abs(delta0)))
    tmp=delta0
    delta0 = (delta - dfi)*(1 + A0*(vel/v0) + s2*B0*(vel/v0)**2)
    print(delta0[0], v0[0])
    if all(abs((delta0-tmp)/delta0) < 1e-6): break

  return delta0

############################################################
# Return S-function suitable for fit_res003 method
#
def res2_dfunc(res2file):
  # name press field   A B s2 dfi   err(A, B, s2, dfi)
  (A,B,s2,dfi) = numpy.loadtxt(res2file, usecols=(3,4,5,6), unpack=1)

  def dfunc(delta0):
    v0 = B/(A - numpy.log(abs(delta0)))
    delta = delta0/(1 + A0*(vel/v0) + s2*B0*(vel/v0)**2) + dfi

  return sfunc
