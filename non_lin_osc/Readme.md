## Response of a non-linear oscillator driven by periodic force

```
x'' + func(x,x',pars) = F*cos(w*t)
```

#### Solving DE with periodic BC

This is the most general approach for finding response of a non-linear oscillator.

* `res = osc_solve_per_func(func, pars, F, w)`

Solve the oscillator equation on one period using periodic BC (`x(0)=x(T)`, `x'(0)=x'(T)`), return
result from solve_bvp(). Use `res.sol(t)[0]` to access smooth function, `res.x` and `res.y` for
final mesh nodes and function values.


* `osc_solve_per_harm(res, N)`

Calculate N-th harmonic of a function returned by osc_solve_per().


* `osc_solve_per(func, pars, F, w)`

Calculate response of a non-linear oscillator on the driving frequency.
It's just a combination of osc_solve_per_func() and osc_solve_per_harm() calls.


