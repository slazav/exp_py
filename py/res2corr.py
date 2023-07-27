# B-phase corrections
# See data_f4_wire_b repository
# V.Zavjalov, 14.06.2023

# Data: see data_f4_wire_b/fit_one_press
data = {
'w1a_00':  [17.404997, 0.313306, 0.486618, 13.839619, 0.872453],
'w1a_02':  [16.005874, 0.305199, 0.161843, 26.812043, 0.836122],
'w1a_05':  [19.618373, 0.778049, 0.065260, 30.108355, 0.958040],
'w1b_00':  [10.514851, 1.605444, 0.080513, 7.9653690, 1.047847],
'w1b_02':  [12.160643, 2.557796, 0.000179, 13.847249, 1.018225],
'w1b_05':  [12.485458, 3.748643, -0.02095, 18.841482, 1.050879],
'w2a_00':  [10.483333, 1.068978, 0.188191, -0.342545, 1.109421],
'w2a_02':  [11.748810, 1.663515, 0.058767, 8.2609610, 1.103366],
'w2a_05':  [12.361161, 2.202467, 0.058756, 8.7126320, 1.092822],
'w1bt_00': [11.326197, 1.014530, 0.021374, 22.320859, 1.127634],
'w1bt_02': [11.259790, 1.195596, 0.071972, 15.704778, 1.074271],
'w1bt_05': [9.8686580, 1.358947, 0.039289, 18.763331, 1.100349],
'w2bt_00': [9.7123850, 0.350739, 0.391738, 19.987295, 1.101966],
'w2bt_02': [11.559636, 0.534391, 0.200483, 26.139833, 1.045476],
'w2bt_05': [10.548072, 0.589810, 0.178895, 27.613787, 1.063568],
}

import numpy

############################################################
def get_data(name, P):
  key = "%s_%02.0f"%(name, round(P))
  print(key)
  if key not in data:
    raise Exception("unknown wire name or pressure: ", name, P)
  return data[key][:5]


############################################################
# Convert vel+delta to delta0.
# (Extrapolation of velocity-dependent width to zero velocity).
# Useful for processing tracking mode.
#
def delta_to_delta0(vel, delta, name, P, field):
  (A,B,dfi,dff,p0) = get_data(name, P)
  A0=0.486387040089067
  delta0=delta
  for i in range(100):
    v0 = B/(A - numpy.log(abs(delta0)))
    tmp=delta0.copy()
    delta0 = (delta - dfi - dff*field**2)*(1 + A0*(vel/v0)**p0)
    if all(abs((delta0-tmp)/delta0) < 1e-6): break

  return delta0

############################################################
def delta_to_ttc(vel, delta, name, P, field):
  (A,B,dfi,dff,p0) = get_data(name, P)
  A0=0.486387040089067

  delta0=delta
  for i in range(100):
    v0 = B/(A - numpy.log(abs(delta0)))
    tmp=delta0
    delta0 = (delta - dfi - dff*field**2)*(1 + A0*(vel/v0)**p0)
    print(delta0[0], v0[0])
    if all(abs((delta0-tmp)/delta0) < 1e-6): break

  # 4-th power polyfit of gap/kTc vs P [bar]
  gap_p = (-1.8629e-07, 1.4784e-05, -4.6497e-04, 8.8839e-03, 1.7725e+00)
  gap = numpy.polyval(gap_p, P)
  return gap/(A - numpy.log(abs(delta0)))

############################################################
# Return S-function suitable for fit_res003 method
#
def dfunc(delta, name, P, field):
  (A,B,dfi,dff,p0) = get_data(name, P)
  A0=0.486387040089067
  def dfunc(delta0):
    v0 = B/(A - numpy.log(abs(delta0)))
    delta = delta0/(1 + A0*(vel/v0)**p0) + dfi + dff*field**2

  return dfunc
