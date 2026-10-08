# SPEC — EPB_AP.jl (first- and second-order schemes)

Source: Arun & Ghorai, *A collocated finite volume scheme for the Euler–Poisson–Boltzmann system*
(referee version 1). Equation numbers refer to that manuscript; "[1]" is Arun & Ghorai, CAMWA 185 (2025).

## 1. Equations (1)–(3)

```math
\partial_t \rho + \partial_x(\rho u) = 0, \qquad
\partial_t(\rho u) + \partial_x(\rho u^2) = -\rho \partial_x \phi, \qquad
-\lambda^2 \partial_{xx} \phi + e^{\phi} = \rho .
```

$\lambda$ is the scale parameter (Debye length) of the Poisson equation; the manuscript writes $\varepsilon$, and the code
keyword is `lambda`. Interval $[a,b]$, $N$ uniform cells, centres $x_i = a + (i-\tfrac12)\Delta x$, all unknowns cell
centred (6). All problems are periodic except the Riemann problem (walls, §7).

## 2. Discrete framework (§2)

* Edge $e = i+\tfrac12$, indices mod $N$. $(\partial_E q)_e = (q_{i+1}-q_i)/\Delta x$,
  $(\mathrm{div}_M v)_i = (v_{i+1/2}-v_{i-1/2})/\Delta x$, $\Delta_M = \mathrm{div}_M \circ \partial_E$ (8).
* Log mean $\mathrm{LM}(a,b) = (b-a)/(\log b - \log a)$, $\mathrm{LM}(a,a) = a$ (9).
  Harmonic edge mean $1/\rho^h_e = \tfrac12(1/\rho_i + 1/\rho_{i+1})$ (11).
* $\lbrace u\rbrace_e = \tfrac12(u_i+u_{i+1})$, $\Delta u_e = u_{i+1}-u_i$, $a^\pm = \tfrac12(a \pm \lvert a\rvert)$ (12).
* Energy (14):
```math
E_h = \sum_i \Delta x \tfrac12 \rho_i u_i^2 + \sum_i \Delta x e^{\phi_i}(\phi_i - 1)
      + \frac{\lambda^2}{2} \sum_e \Delta x (\partial_E \phi)_e^2 .
```
* Source (15): $S_i = -\tfrac12\big[\bar\rho_{i+1/2}(\partial_E\phi)_{i+1/2} + \bar\rho_{i-1/2}(\partial_E\phi)_{i-1/2}\big]$.

## 3. Fluxes (17)–(19)

$\bar\rho_e = \mathrm{LM}(\rho_i, \rho_{i+1})$; $F_e = \bar\rho_e \lbrace u\rbrace_e - Q_e$;
$G_e = u^{up}_e F_e$ with $u^{up}_e = u_i$ if $F_e \ge 0$, else $u_{i+1}$.

## 4. Scheme A — §3.3 (Theorem 3.5), `scheme = :shift`

Step (27)–(29):
```math
\rho^{n+1} = \rho^n - \Delta t \mathrm{div}_M F, \qquad
(\rho u)^{n+1} = (\rho u)^n - \Delta t \mathrm{div}_M G + \Delta t S, \qquad
-\lambda^2 \Delta_M \phi^{n+1} + e^{\phi^{n+1}} = \rho^{n+1},
```
with $F, G$ built from $\bar\rho^n, u^n$ and $Q_e = \eta_e \Delta t (\partial_E \phi^{n+1})_e$; $S$ from $\bar\rho^n, \phi^{n+1}$.
* $\eta_e = \tfrac54 (\bar\rho^n_e)^2/\rho^{n,h}_e$ (43), edge by edge (local) in both schemes. The global variant of §5.2
  ($\max_e$ of (43) at every edge) is not implemented.
* Time step (44), with $X = \Delta t/\Delta x$ and the minimum over $i$:
```math
X \le \frac{\tfrac25 \rho_i}{A_i + \sqrt{A_i^2 + \tfrac45 C_i \rho_i}}, \qquad
A_i = \sum_{e \in \mathcal N(i)} \bar\rho_e \lvert\lbrace u\rbrace_e\rvert, \qquad
C_i = \Delta x \sum_{e \in \mathcal N(i)} \eta_e \lvert(\partial_E \phi^n)_e\rvert
```
  ($C_i$ with $\Delta x$: see §9 item 1).
* A posteriori (CFL): $\frac{\Delta t}{\Delta x}\big(\lvert F_{i+1/2}\rvert + \lvert F_{i-1/2}\rvert\big) \le \rho^n_i/5$
  with the computed $F$; if violated, repeat the step with $\Delta t$ multiplied by the smallest ratio over $i$.

## 5. Scheme B — §3.4, linearised variant (52), `scheme = :dissipation`

Same as A except (46)
```math
Q_e = \mu_e (\partial_E \phi^{n+1})_e - \kappa \lambda^2 (\partial_E \Delta_M \phi^{n+1})_e, \qquad
\mu_e = \eta_e \Delta t + \kappa \hat\rho^n_e, \qquad
\hat\rho^n_e = \mathrm{LM}\big(e^{\phi^n_i}, e^{\phi^n_{i+1}}\big), \qquad \nu_e = 0 .
```
* $\theta = \tfrac12$; local $\eta_e = (\bar\rho^n_e)^2 / \big(2(1-\theta)\rho^{n,h}_e\big)$ (Cor. 3.17(i), edge by edge as in §4);
  $\kappa^n = \tfrac{\Delta x}{2} \max_i \lvert u^n_i\rvert$ (Remark 3.23).
  For §5.6 "choice of $\kappa$": $\kappa^n = c_\kappa \Delta x (\max\lvert u^n\rvert + s_\kappa)$,
  $(c_\kappa, s_\kappa) \in \lbrace(\tfrac12,1), (\tfrac12,0), (\tfrac14,0), (\tfrac1{10},0)\rbrace$.
* Time step (Remark 3.22): $c_e = \bar\rho^n_e \lbrace u^n\rbrace_e - \kappa(\partial_E\rho^n)_e$; outward
  $c_{i,e} = +c_e$ ($e = i+\tfrac12$), $-c_e$ ($e = i-\tfrac12$);
```math
A_i = \sum_e c^+_{i,e}, \qquad C_i = \Delta x \sum_e \eta_e \lvert(\partial_E\phi^n)_e\rvert, \qquad
X_i = \frac{2\theta\rho_i}{A_i + \sqrt{A_i^2 + 4\theta C_i \rho_i}}, \qquad
\Delta t = 0.9 \min_i X_i \Delta x .
```
  No other bound on $\Delta t$ (see §9 item 10).
  A posteriori (CFLθ): $b_i = \frac{\Delta t}{\Delta x} \sum_e F^+_{i,e} \le \theta \rho^n_i$; retry as in §4.

## 6. Implicit solve (Remark 3.22), both schemes

$\rho^* = \rho^n - \Delta t \mathrm{div}_M(\bar\rho^n \lbrace u^n\rbrace)$. Solve for $\phi = \phi^{n+1}$ (scheme A: $\mu = \eta\Delta t$, $\kappa = 0$):
```math
e^{\phi} - \lambda^2 \Delta_M \phi - \Delta t \mathrm{div}_M(\mu \partial_E \phi) + \kappa\lambda^2 \Delta t \Delta_M^2 \phi = \rho^* \qquad (52)
```
*Derivation.* The discrete Poisson equation (29), the discrete form of (3), holds at level $n+1$:
$-\lambda^2\Delta_M\phi^{n+1} + e^{\phi^{n+1}} = \rho^{n+1}$. Insert (27) with the flux (18) and $Q$ from (46):
```math
\rho^{n+1} = \rho^n - \Delta t \mathrm{div}_M\big(\bar\rho^n\lbrace u^n\rbrace - \mu\partial_E\phi^{n+1} + \kappa\lambda^2\partial_E\Delta_M\phi^{n+1}\big)
            = \rho^* + \Delta t \mathrm{div}_M(\mu\partial_E\phi^{n+1}) - \kappa\lambda^2\Delta t \Delta_M^2\phi^{n+1},
```
using $\mathrm{div}_M\partial_E = \Delta_M$. Substituting this for $\rho^{n+1}$ in (29) and moving the $\phi^{n+1}$ terms to the left gives (52).
Thus $\rho^{n+1}$ is eliminated and (52) is one equation for $\phi^{n+1}$ alone; the explicit part of the mass flux is in $\rho^*$.

Newton with Jacobian $\mathrm{diag}(e^\phi) + \lambda^2 L + \Delta t L_\mu + \kappa\lambda^2\Delta t L^2$, $L = -\Delta_M$,
periodic corners included (walls: rows of the boundary cells without the outer edge); sparse SPD solve. Then $F$ from (18), $\rho^{n+1}$ from (27), so (29) holds by construction;
$u^{n+1} = (\rho u)^{n+1}/\rho^{n+1}$. $\phi^0$ from $\rho^0$ by (29) (same solver, $\Delta t = 0$).

## 7. Test cases ($\phi^0$ from (29))

The parameters of each case (domain, boundary, $N$, $\lambda$, $T$, scheme, time step, $n_r$, $c$) are read from
`params/<case>.toml`, not from the source; the table gives their values.

Periodic, except the Riemann problem, which uses the boundary conditions of [1, §7.4] at both ends: extrapolation for
$\rho$, no-slip for $u$, homogeneous Neumann for $\phi$. Discretely, the two boundary edges carry $F = G = 0$ and
$(\partial_E\phi) = 0$ (so also $Q = 0$ and no source contribution there); with $u = 0$ at the wall the extrapolated
$\rho$ enters no flux. Mass and the energy identity are kept. Up to $T = 50$ no wave reaches a wall (rarefaction head
at $x = -t \ge -50$, shock at $u_s t \le 59$), so the walls only see exponentially small tails.

| case | domain | $\rho^0$, $u^0$ | $\lambda$ | $T$ | $N$ | $\Delta t$ | reference |
|---|---|---|---|---|---|---|---|
| smooth (EulerAP.jl, not §5.2) | $[-10,10]$ | $1-0.3\mathrm{sech}(2x)$, $U = 0.2$ | 1 | 1 | 50; order 100…800 | (44) | fine grid $N=3200$, cell-averaged |
| steep (§5.2) | $[0,1]$ | $1+0.8\sin^3 2\pi x$, $0.25\cos 2\pi x$ | 1, 1e-3 | 0.5 | 200 | (44) | none |
| near vacuum (§5.2) | $[0,1]$ | $10^{-4}+(1-10^{-4})e^{-60(x-1/2)^2}$, 0 | 1, 1e-3 | 0.5 | 200 | (44) | none |
| bump (§5.6) | $[0,1]$ | $1+\tfrac12 e^{-(x-1/2)^2/0.005}$, 1 | 1, 0.1 | 2 | 200 | scheme rule | none |
| step (§5.6) | $[0,1]$ | 1 on $(\tfrac14,\tfrac34)$, $\tfrac12$ else; 0 | 0.1 | 1 | 200, 800, 3200 | scheme rule | $\min\rho$ stays $\tfrac12$ |
| expansion (§5.6) | $[0,1]$ | 1, $\sin 2\pi x$ | 1, 0.1 | 0.6 | 200 (400, 1600 at $\lambda=0.1$) | scheme rule | $\rho(0,T) \to 1/(1+2\pi T)$ |
| soliton (§5.6) | $[0,40]$ | travelling wave $c=1.3$, peak at $x=20$ | 1 | order 2; robustness $L/c$ | 400; order 100…1600 | rule, $\le 0.4\Delta x$ (order) | exact (translation) |
| ap (§5.6) | $[0,1]$ | $1+0.2\sin 2\pi x$, $0.2\cos 2\pi x$ | 1 … 1e-6 | 0.2 | 200 | fixed $0.25\Delta x$ | $\lVert e^\phi-\rho\rVert_\infty \propto \lambda^2$ |
| simple wave (§5.6) | $[0,1]$ | $e^{u^0}$, $0.2\sin 2\pi x$ | 1e-4 | 0.3 | 100…1600 | rule, $\le 0.4\Delta x$ | exact (characteristics) |
| riemann ([1, §7.4]) | $[-80,100]$, walls | 1 for $x<0$, $n_r$ for $x \ge 0$; $u=0$; $n_r \in \lbrace 0.5, 0.75, 0.95\rbrace$ | 1e-4 | 50 | 400 (9000 in [1]) | (44) | ICE solution [1, (7.8)–(7.10)] |
| riemann_speed (§5.5) | as above | as above, $n_r = 0.5$ | 1e-4 | 8 (snapshot 4) | 2000, 4000, 8000 | (44) | $u_s$ from [1, (7.10)] |

* Data files: `data/<case>_<scheme>_N<N>_lambda<λ>_T<T>.csv` (x, ρ, u, φ and exact columns if any), `…_t<t>.csv`
  snapshots, `…_history.csv` (t, Δt, $E_h$, min ρ, mass, $\lVert e^\phi-\rho\rVert_\infty$, CFL ratio); header lines `# key = value`.
* Soliton: $u = c - \sqrt{c^2-2\phi}$, $\rho = c/\sqrt{c^2-2\phi}$, $\lambda^2(\phi')^2 = 2V(\phi)$,
  $V = e^\phi - 1 + c\big(\sqrt{c^2-2\phi} - c\big)$; $\phi_{max}$ from $V(\phi_{max}) = 0$. Computed by integrating
  $\lambda^2\phi'' = e^\phi - \rho(\phi)$ from $(\phi_{max}, 0)$ with RK4 ($h = 10^{-3}$) and Hermite interpolation;
  the exact solution at $t$ is the profile at $\xi = x - 20 - ct$ (periodic).
* Simple wave: $u(x,t) = u^0(\xi)$, $x = \xi + (u^0(\xi)+1)t$ (Newton in $\xi$), $\rho = e^u$, $\phi = \log\rho$.
* smooth: the at-rest dip of EulerAP.jl moved with the uniform drift $U$ (Galilean invariance on the periodic domain), so
  $u > 0$ everywhere: (44) gives $\Delta t \propto \Delta x$ (from rest, (90): $\propto \sqrt{\Delta x}$) and there is no stagnation point.
* Riemann: $u_m$ solves $(1 - n_r e^{u_m})(u_m^2 - 2u_m - 2\log n_r) - 2u_m^2 = 0$, $u_s = u_m/(1 - n_r e^{u_m})$.
  Shock speed $(x_s(8) - x_s(4))/4$, $x_s$ by equal areas over $[x_a, b]$, $x_a = \tfrac12(u_m - 1 + u_s)t$.
* Unit-test states (no time integration, §5.3 order-1 column): smooth, steep ($u^0 = 0.5\cos$), near vacuum,
  shock-like ($1, 0.4 \mid 0.3, -0.4$ at $x = \tfrac12$), random rough ($\rho = 0.5+r$, $u = 2s-1$, `Random.seed!(20260923)`).

## 8. Properties checked by tests

Mass conservation to round-off; positivity $\rho^{n+1} \ge \tfrac45\rho^n$ (A) / $(1-\theta)\rho^n$ (B);
$E_h^{n+1} \le E_h^n$ per step (tolerance $10^3\epsilon_{mach}\max(1,\lvert E^0\rvert)$); semidiscrete rate (93) equals (23)
with $D \equiv 0$; one-step random states (Theorem 3.16); $\lVert e^\phi - \rho\rVert_\infty$ falls ≈100× per decade of
$\lambda$; observed order ≈ 1.

## 9. Unclear items and the reading used

1. (44) prints $C_i$ without $\Delta x$; the derivation, Remark 3.22 and (90) need
   $C_i = \Delta x \sum_e \eta_e \lvert\partial_E\phi\rvert$. **Used: with $\Delta x$** (author's decision). Note: the
   printed form reproduces the §5.2 tables to all digits, so those step counts and energy changes differ here;
   all properties (positivity, energy decrease, AP, order) hold with either form.
2. Soliton robustness $T$ is missing: **$T = L/c = 40/1.3$** (one domain traversal, Degond et al. $t_L$, as instructed).
3. Riemann: **all data from [1, §7.4]** to reproduce [1, Fig. 8]: domain $[-80,100]$, 9000 cells (default 400 for
   quick runs; set `N` in `params/riemann.toml`), $\lambda = 10^{-4}$,
   $n_r \in \lbrace 0.5, 0.75, 0.95\rbrace$, $T = 50$, the boundary conditions of [1] (§7), log mean only, scheme A.
   The manuscript's §5.5 shock-speed study ($n_r = 0.5$, $N$ = 2000–8000, $t$ = 4, 8) uses the same problem with the same
   walls. Shock positions at $T = 50$: 58.9, 53.6, 50.6 for $n_r$ = 0.5, 0.75, 0.95.
4. §5.6 uses the linearised variant (52); the literal variant (45) is not implemented (same within 3–4 digits per the
   manuscript). $\nu_e = 0$.
5. Soliton: the text writes $\lambda^2(\phi')^2 = V$ but the quadrature $\xi = \lambda\int d\phi/\sqrt{2V}$ and a
   direct derivation give $\lambda^2(\phi')^2 = 2V$; used $2V$.
6. Newton stop: $\lVert\delta\phi\rVert_\infty \le 10^{-12}(1+\lVert\phi\rVert_\infty)$, or
   $\lVert R\rVert_\infty \le 10^{-14} s (1+\lVert\phi\rVert_\infty)$ with $s = \max \mathrm{diag} J$ (the stencil
   scale of Remark 4.13); max 50 iterations; each update damped to $\lVert\delta\phi\rVert_\infty \le 1$.
7. Collapse: $\Delta t < 10^{-4}\Delta x$ stops the run (§5.6). The last step is shortened to hit $T$ (and snapshot times).
8. Random one-step states (§5.1) are loosely specified: 48 cells, $\log\rho \sim N(0,\sigma^2)$, $\sigma \in (0,3]$,
   optionally smoothed by a 3-point average; $u = A\cdot N(0,1)$, $A \in (0,3]$; $\lambda$ log-uniform in
   $[10^{-4}, 10]$; $\theta \in \lbrace 0.2, 0.5, 0.8\rbrace$; $\kappa = 10^U \cdot \tfrac{\Delta x}{2}(\max\lvert u\rvert+1)$,
   $U \sim U(-1,1)$; $\Delta t$ halved from 1 until (CFLθ) holds.
9. Errors: discrete norms $\Delta x\sum_i\lvert e_i\rvert$ ($L^1$), $(\Delta x\sum_i e_i^2)^{1/2}$ ($L^2$),
   $\max_i\lvert e_i\rvert$ ($L^\infty$) at cell centres; fine-grid references are restricted by averaging blocks of
   $N_{ref}/N$ cells.
10. §5.6 caps the Remark 3.22 step by $\Delta t \le 0.9\Delta x/(\max_i\lvert u^n_i\rvert + 1)$. **Not used** (author's
    decision): both schemes take $\Delta t$ from their $X$ bound only. Remark 3.23 needs $\max\lvert u\rvert\Delta t/\Delta x \le 1$
    for $\kappa^n$ to satisfy (53); it is $\le 0.45$ in all runs of §7.

## 10. Second-order scheme — §4, `SecondOrder()`

Unknowns $(\rho, q)$, $q = \rho u$; $v = q/\rho$. Code: `SecondOrder(; limiter = :minmod, c_D = 1.0, cfl = 0.4)`;
in a parameter file `kind = "second"` under `[scheme]`. Every test case runs with it.

* Reconstruction (54)–(55) on $(\rho, v)$: $w^-_{i+1/2} = w_i + \tfrac12\Delta x\sigma_i$,
  $w^+_{i+1/2} = w_{i+1} - \tfrac12\Delta x\sigma_{i+1}$, $\sigma$ minmod (`limiter = :minmod`) or central (`:none`).
* Fluxes at a state $(\rho, q)$, with $\phi$ (56)–(63):
```math
\bar\rho_e = \mathrm{LM}(\rho^-_e, \rho^+_e), \quad E_e = \bar\rho_e \lbrace v\rbrace_e, \quad
D_e = c_D \Delta x^2 (\partial_E \rho)_e \ (72), \quad
G_e = v^-_e F_e^+ + v^+_e F_e^- \ (61),
```
  $S_i = -\tfrac12\big[(\bar\rho\partial_E\phi)_{i+1/2} + (\bar\rho\partial_E\phi)_{i-1/2}\big]$ (15), $C = \mathrm{div}_M G$.
  The momentum flux is (61) only (author's decision); (64) is not implemented.
* ARS(2,2,2) (73)–(74): $\gamma = 1 - 1/\sqrt2$, $\delta = 1 - 1/(2\gamma)$. A stage mass flux is
  $F = E - \gamma\Delta t\bar\rho\partial_E\phi - D$, with $E$, $\bar\rho$ and $D$ from the stage's explicit state.
* Stage sequence (76)–(87), from $(\rho^n, q^n, \phi^n)$:
  1. $F^{(1)} = E - D$, $C^{(1)}$, $S^{(1)}$ at $(\rho^n, q^n, \phi^n)$.
  2. $\rho^{(2)}_E = \rho^n - \gamma\Delta t \mathrm{div}_M F^{(1)}$, $\tilde q^{(2)} = q^n - \gamma\Delta t C^{(1)}$,
     $q^{(2)}_E = \tilde q^{(2)} + \gamma\Delta t S^{(1)}$. $E^{(2)}$, $D^{(2)}$ at $(\rho^{(2)}_E, \tilde q^{(2)})$;
     solve (88) with $\hat\rho^{(2)} = \rho^n - \gamma\Delta t \mathrm{div}_M(E^{(2)} - D^{(2)})$ for $\phi^{(2)}$;
     $F^{(2)}$, $S^{(2)}$ with $\phi^{(2)}$.
  3. $C^{(2)}$ from $F^{(2)}$ and the reconstruction of $(\rho^{(2)}_E, q^{(2)}_E)$ (81).
     $\rho^{(3)}_E = \rho^n - \Delta t \mathrm{div}_M(\delta F^{(1)} + (1-\delta)F^{(2)})$,
     $\hat q^{(3)} = q^n - \Delta t(\delta C^{(1)} + (1-\delta)C^{(2)}) + (1-\gamma)\Delta t S^{(2)}$.
     $E^{(3)}$, $D^{(3)}$ at $(\rho^{(3)}_E, \hat q^{(3)})$; solve (88) with
     $\hat\rho^{(3)} = \rho^n - \Delta t \mathrm{div}_M\big((1-\gamma)F^{(2)} + \gamma(E^{(3)} - D^{(3)})\big)$ for $\phi^{n+1}$;
     $\rho^{n+1} = \rho^n - \Delta t \mathrm{div}_M\big((1-\gamma)F^{(2)} + \gamma F^{(3)}\big)$,
     $q^{n+1} = \hat q^{(3)} + \gamma\Delta t S^{(3)}$.
* (88) is the solver of §6 with $\mu_e = \gamma^2\Delta t \bar\rho_e$ and $\kappa = 0$, so that
  $e^{\phi} - \lambda^2\Delta_M\phi = \rho^{(k)}$ holds for the stage density.
* Time step (§4.11): $\Delta t = \min\big((44),\ (91)\big)$, both with $\eta$ of (43); in (91)
  $\Delta t \le \mathrm{cfl} \Delta x / \max_i\big(\lvert v_i\rvert + \eta_i \max(\lvert\phi_i - \phi_{i-1}\rvert, \lvert\phi_{i+1} - \phi_i\rvert)/\Delta x\big)$,
  $\eta_i$ the larger of the two edge values, cfl = 0.4. A step that makes any stage or reconstructed density
  non-positive is repeated with $\Delta t/2$; $\Delta t < 10^{-4}\Delta x$ is a collapse. Order studies use
  `TimeStep(mode = :fixed, cfl = 0.4)` (§5.4).
* Walls: $\sigma = 0$ in cells 1 and $N$; $F = G = 0$ and $\partial_E\phi = 0$ on the walls, as in §7.
* Tests: the rate (93) equals (68) on random states (with the $\beta_e$ of Remark 4.2, both limiters, $c_D \in \lbrace 0, 1, 10\rbrace$);
  $-\mathrm{div}_M F \to -\partial_x(\rho u)$ at order 2 (Proposition 4.9); mass to round-off and $\rho > 0$ on every case;
  $\lVert e^\phi - \rho\rVert_\infty \propto \lambda^2$; fully discrete $L^1$ order ≈ 2 for $\rho$, $m$, $\phi$ (smooth,
  simple wave, soliton; $\Delta t = 0.4\Delta x$). No energy test: §4.13 claims no fully discrete inequality.

## 11. Second order: unclear items and the reading used

1. (76)–(87) subtract $\gamma\Delta t S$ from $q$, while (63) and (66) give $\dot q = -C + S$. **Used:** $q \mathrel{+}= \Delta t S$
   (the force $-\rho\partial_x\phi$ enters with the sign of (66)); Remark 4.11 is consistent with this.
2. Which $\bar\rho$ in (59), (88) and $S^{(k)}$: **the log mean of the reconstruction of the stage's explicit state**
   (stage 1: $\rho^n$; stage 2: $\rho^{(2)}_E$; stage 3: $\rho^{(3)}_E$), the same at a face in $E$, $D^{si}$, (88) and $S$,
   as Remark 4.12 and Condition 2 require.
3. $D$ is not placed in the stage sequence. **Used:** $D^{(k)}$ at the stage's explicit density, inside $F^{(k)}$ and $\hat\rho^{(k)}$,
   so that the Poisson equation holds for $\rho^{(k)}$.
4. (81) lists $C^{(2)} = C(\rho^{(2)}_E, q^{(2)}_E, \phi^{(2)})$. **Used:** $G^{(2)}$ from the mass flux $F^{(2)}$ and
   the velocities reconstructed from $(\rho^{(2)}_E, q^{(2)}_E)$.
5. (91) names neither cfl nor $\eta$. **Used:** cfl = 0.4 and $\eta$ of (43) (author's decision). The positivity retry
   is not in the manuscript (author's decision).
6. (93): the kinetic part is $u_i\dot q_i - \tfrac12 u_i^2\dot\rho_i$.
7. Expansion at $\lambda = 1$: the flow converges on $x = \tfrac12$ and the near-pressureless solution concentrates
   there from $t \approx 1/(2\pi)$ (first order reaches $\max\rho \approx 130$ at $T = 0.6$). Second order has no
   positivity guarantee and collapses at $t \approx 0.17$; at $\lambda = 0.1$ it completes. The manuscript runs
   this case with first order only (§5.6).
