## Response of a non-linear oscillator driven by periodic force

```
x'' + func(x,x',pars) = F*cos(w*t)
```

#### Examples of oscillators

* harmonic osc, `osc_harm(x,dx, (w0,tau))`

Parameters: 
- `w0` -- resonant frequency
- `tau` -- relaxation time

Formula: `func(x,v) = w0**2*x + v/tau`

* duffing osc, `osc_duffing(x,dx, (w0,tau,a))`

- `w0`  - resonant frequency
- `tau` - relaxation time at low drive
- `a`   - non-linear parameter

Formula: `func(x,v) = w0**2*x + v/tau + a*x**3`

* pseudoplastic osc N1, `osc_pseudopl1(x,dx, (w0,tau,vc,k))`

Formula: `func(x,v) = w0**2*x + v/tau * (k - (1-k) * vc/numpy.sqrt(vc**2 + v**2))`

* pseudoplastic osc N2, `osc_pseudopl2(x,dx, (w0,tau,vc))`

Formula: `func(x,v) = w0**2*x + v/tau * vc/numpy.sqrt(vc**2 + v**2)`



#### Solving DE with periodic BC

This is the most general approach for finding response of a non-linear oscillator.
It is slow and in some cases not accurate.

* `res = osc_solve_per_func(func, pars, F, w)`

Solve the oscillator equation on one period using periodic BC (`x(0)=x(T)`, `x'(0)=x'(T)`), return
result from solve_bvp(). Use `res.sol(t)[0]` to access smooth function, `res.x` and `res.y` for
final mesh nodes and function values.


* `osc_solve_per_harm(res, N)`

Calculate N-th harmonic of a function returned by osc_solve_per().


* `osc_solve_per(func, pars, F, w)`

Calculate response of a non-linear oscillator on the driving frequency.
It's just a combination of osc_solve_per_func() and osc_solve_per_harm() calls.


#### Averaging motion in Van-der-Pol coordinates

This is recommended way of solving non-linear oscillators. The idea is to
use van der Pol coordinates (u,v) = (x, x'/w rotated by w*t), and then
numerically average equations over rotation period (phase 0..2*pi) And
then find equilibrium, (u',v')=0.

See texts:
- http://www.scholarpedia.org/article/Duffing_oscillator
- http://www.scholarpedia.org/article/Equilibrium

* `osc_solve_vdp_eq(uv, func, pars, F, w)`

Find (\dot u,\dot v) as a function of (u,v) for period-averaged motion of
the periodically-driven non-linear oscillator.

* `osc_solve_vdp(func, pars, F, w, a0=0, p0=0)`

Find equilibrium (zero of nonlin_osc_eq function).
This is enough for simple small non-linearities,
for duffing oscillator it does not work properly
(one should integrate trajectories in u-v space instead)

#### Duffing oscillator

Formulas for Duffing oscillator. They should be fully equivalent to `osc_solve_vdp_*`,
but faster.

* `osc_solve_vdp_duff_eq(uv, w0, tau, a, F, w)` -- Equivalent of `osc_solve_vdp_eq`.

* `osc_solve_vdp_duff(w0, tau, a, F, w, a0=0, p0=0)` -- Equivalent of `osc_solve_vdp`

* `osc_duff_amp(w,F,w0,tau, a)` -- Analytical expression for amplitude of Duffing oscillator.
(Solve qubic equation `((w**2 - w0**2 - 3/4.0*a * A**2)**2 + (w/tau)**2)*A**2 = F**2`,
return only real roots.)

