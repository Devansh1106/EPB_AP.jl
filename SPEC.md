# SPEC — EPB_AP.jl (first-order schemes)

Source: Arun & Ghorai, *A collocated finite volume scheme for the Euler–Poisson–Boltzmann system*
(referee version 1). Equation numbers refer to that manuscript; "[1]" is Arun & Ghorai, CAMWA 185 (2025).

## 1. Equations (1)–(3)

    ∂t ρ + ∂x(ρu) = 0,   ∂t(ρu) + ∂x(ρu²) = −ρ ∂x φ,   −ε² ∂xx φ + e^φ = ρ.

All problems: periodic interval, N uniform cells, centres x_i = a + (i−½)Δx, all unknowns cell centred (6).

## 2. Discrete framework (§2)

* Edge e = i+½, indices mod N. (∂E q)_e = (q_{i+1}−q_i)/Δx, (divM v)_i = (v_{i+½}−v_{i−½})/Δx, ΔM = divM∘∂E (8).
* Log mean LM(a,b) = (b−a)/(log b − log a), LM(a,a)=a (9). Harmonic edge mean 1/ρʰ_e = ½(1/ρ_i + 1/ρ_{i+1}) (11).
* {u}_e = ½(u_i+u_{i+1}), Δu_e = u_{i+1}−u_i, a^± = ½(a±|a|) (12).
* Energy (14): E_h = Σ Δx ½ρ_i u_i² + Σ Δx e^{φ_i}(φ_i−1) + (ε²/2) Σ_e Δx (∂E φ)_e².
* Source (15): S_i = −½[ ρ̄_{i+½}(∂E φ)_{i+½} + ρ̄_{i−½}(∂E φ)_{i−½} ].

## 3. Fluxes (17)–(19)

ρ̄_e = LM(ρ_i, ρ_{i+1}); F_e = ρ̄_e{u}_e − Q_e; G_e = u^up_e F_e with u^up_e = u_i if F_e ≥ 0, else u_{i+1}.

## 4. Scheme A — §3.3 (Theorem 3.5), `scheme = :shift`

Step (27)–(29): ρ^{n+1} = ρ^n − Δt divM F; (ρu)^{n+1} = (ρu)^n − Δt divM G + Δt S; −ε²ΔM φ^{n+1} + e^{φ^{n+1}} = ρ^{n+1};
F, G with ρ̄^n, u^n and Q_e = η_e Δt (∂E φ^{n+1})_e; S with ρ̄^n, φ^{n+1}.
* η_e = (5/4)(ρ̄^n_e)²/ρ^{n,h}_e (43). Option `eta = :global`: max_e of (43) at every edge (§5.2).
* Time step (44): X = Δt/Δx ≤ (2/5)ρ_i / (A_i + √(A_i² + (4/5)C_i ρ_i)), A_i = Σ_{e∈N(i)} ρ̄_e|{u}_e|,
  **C_i = Δx Σ_{e∈N(i)} η_e |(∂E φ^n)_e|** (see §8 item 1). Root used as is; min over i.
* A posteriori (CFL): (Δt/Δx)(|F_{i+½}|+|F_{i−½}|) ≤ ρ^n_i/5 with the computed F; if violated, repeat the step
  with Δt multiplied by the smallest ratio (ρ_i/5)/(…) over i.

## 5. Scheme B — §3.4, linearised variant (52), `scheme = :dissipation`

Same as A except Q_e = μ_e(∂E φ^{n+1})_e − κε²(∂E ΔM φ^{n+1})_e (46), μ_e = η_e Δt + κ ρ̂^n_e,
ρ̂^n_e = LM(e^{φ^n_i}, e^{φ^n_{i+1}}), ν_e = 0.
* θ = ½; η_e = (ρ̄^n_e)²/(2(1−θ)ρ^{n,h}_e) (Cor. 3.17(i)); κ^n = (Δx/2) max_i |u^n_i| (Remark 3.23).
  For §5.6 "choice of κ": κ^n = c_κ Δx (max|u^n| + s_κ), (c_κ, s_κ) ∈ {(½,1),(½,0),(¼,0),(1/10,0)}.
* Time step (Remark 3.22): c_e = ρ̄^n_e{u^n}_e − κ(∂E ρ^n)_e; outward c_{i,e} = +c_e (e=i+½), −c_e (e=i−½);
  A_i = Σ c⁺_{i,e}, C_i = Δx Σ η_e|(∂E φ^n)_e|, X = 2θρ_i/(A_i + √(A_i² + 4θC_iρ_i)); Δt = 0.9 min_i X_i Δx,
  also Δt ≤ 0.9Δx/(max|u^n|+1). A posteriori (CFLθ): b_i = (Δt/Δx)Σ_e F⁺_{i,e} ≤ θρ^n_i; retry as in §4.

## 6. Implicit solve (Remark 3.22), both schemes

ρ* = ρ^n − Δt divM(ρ̄^n{u^n}). Solve for φ = φ^{n+1}:
e^φ − ε²ΔM φ − Δt divM(μ ∂E φ) + κε²Δt ΔM²φ = ρ* (52)  (scheme A: μ = ηΔt, κ = 0).
Newton with Jacobian diag(e^φ) + ε²L + ΔtL_μ + κε²ΔtL², L = −ΔM, periodic corners included; sparse SPD solve.
Then F from (18), ρ^{n+1} from (27), so (29) holds by construction; u^{n+1} = (ρu)^{n+1}/ρ^{n+1}.
φ⁰ from ρ⁰ by (29) (same solver, Δt = 0).

## 7. Test cases (all periodic; φ⁰ from (29))

| case | domain | ρ⁰, u⁰ | ε | T | N | Δt | reference |
|---|---|---|---|---|---|---|---|
| smooth (§5.2) | [0,1] | 1+0.3 sin2πx, 0.2 cos2πx | 1, 1e-2, 1e-4 | 0.5 (AP: 0.3; order: 0.2) | 200; order 100…800 | (44) | fine grid N=3200, cell-averaged |
| steep (§5.2) | [0,1] | 1+0.8 sin³2πx, 0.25 cos2πx | 1, 1e-3 | 0.5 | 200 | (44) | none |
| near vacuum (§5.2) | [0,1] | 1e-4+(1−1e-4)e^{−60(x−½)²}, 0 | 1, 1e-3 (η local/global: 1, 1e-2) | 0.5 (η study 0.3) | 200 | (44) | none |
| bump (§5.6) | [0,1] | 1+½e^{−(x−½)²/0.005}, 1 | 1, 0.1 | 2 | 200 | scheme rule | none |
| step (§5.6) | [0,1] | 1 on (¼,¾), ½ else; 0 | 0.1 | 1 | 200, 800, 3200 | scheme rule | min ρ stays ½ |
| expansion (§5.6) | [0,1] | 1, sin2πx | 1, 0.1 | 0.6 | 200 (400, 1600 at ε=0.1) | scheme rule | ρ(0,T)→1/(1+2πT) |
| soliton (§5.6) | [0,40] | travelling wave c=1.3, peak at x=20 | 1 | order 2; robustness L/c | 400; order 100…1600 | rule, ≤0.4Δx (order) | exact (translation) |
| ap (§5.6) | [0,1] | 1+0.2 sin2πx, 0.2 cos2πx | 1 … 1e-6 | 0.2 | 200 | fixed 0.25Δx | ‖e^φ−ρ‖∞ ∝ ε² |
| simple wave (§5.6) | [0,1] | e^{u⁰}, 0.2 sin2πx | 1e-4 | 0.3 | 100…1600 | rule, ≤0.4Δx | exact (characteristics) |
| riemann (§5.5, [1, §7.4]) | [−80,100] | 1 / n_r=0.5 at x=0, u=0 | 1e-4 | 8 (snapshot 4) | 2000, 4000, 8000 | (44) | ICE solution [1, (7.8)–(7.10)] |

* Data files: `data/<case>_<scheme>_N<N>_eps<ε>_T<T>.csv` (x, ρ, u, φ and exact columns if any), `…_t<t>.csv`
  snapshots, `…_history.csv` (t, Δt, E_h, min ρ, mass, ‖e^φ−ρ‖∞, CFL ratio); header lines `# key = value`.
* Soliton: u = c − √(c²−2φ), ρ = c/√(c²−2φ), ε²(φ′)² = 2V(φ), V = e^φ − 1 + c(√(c²−2φ) − c);
  φ_max from V(φ_max)=0. Computed by integrating ε²φ″ = e^φ − ρ(φ) from (φ_max, 0) with RK4 (h = 10⁻³) and
  Hermite interpolation; exact solution at t is the profile at ξ = x − 20 − ct (periodic).
* Simple wave: u(x,t) = u⁰(ξ), x = ξ + (u⁰(ξ)+1)t (Newton in ξ), ρ = e^u, φ = log ρ.
* Riemann: u_m solves (1 − n_r e^{u_m})(u_m² − 2u_m − 2 log n_r) − 2u_m² = 0, u_s = u_m/(1 − n_r e^{u_m}).
  Shock speed (x_s(8) − x_s(4))/4, x_s by equal areas over [x_a, 60], x_a = ½(u_m−1+u_s)t.
* Unit-test states (no time integration, §5.3 order-1 column): smooth, steep (u⁰ = 0.5cos), near vacuum,
  shock-like (1, 0.4 | 0.3, −0.4 at x=½), random rough (ρ = 0.5+r, u = 2s−1, `Random.seed!(20260923)`).

## 8. Properties checked by tests

Mass conservation to round-off; positivity ρ^{n+1} ≥ ⅘ρ^n (A) / (1−θ)ρ^n (B); E_h^{n+1} ≤ E_h^n per step
(tolerance 10³ϵ_mach max(1,|E⁰|)); semidiscrete rate (93) equals (23) with D ≡ 0; one-step random states
(Theorem 3.16); ‖e^φ − ρ‖∞ falls ≈100× per decade of ε; observed order ≈ 1.

## 9. Unclear items and the reading used

1. (44) prints C_i without Δx; the derivation, Remark 3.22 and (90) need C_i = Δx Σ η|∂φ|. **Used: with Δx** (author's
   decision). Note: the printed form reproduces the §5.2 tables to all digits, so those step counts and energy
   changes differ here; all properties (positivity, energy decrease, AP, order) hold with either form.
2. Soliton robustness T is missing: **T = L/c = 40/1.3** (one domain traversal, Degond et al. t_L, as instructed).
3. Riemann data from [1, §7.4]: domain [−80,100], ε = 1e-4, n_r = 0.5 only, log mean only; posed **periodic** as all
   problems here (manuscript §5), not with the BCs of [1]; the exact solution ignores the waves from the
   periodic wrap point, which stay outside [−60, 80] for t ≤ 8. Run with scheme A.
4. §5.6 uses the linearised variant (52); the literal variant (45) is not implemented (same within 3–4 digits per the
   manuscript). ν_e = 0.
5. Soliton: the text writes ε²(φ′)² = V but the quadrature ξ = ε∫dφ/√(2V) and a direct derivation give ε²(φ′)² = 2V; used 2V.
6. Newton stop: ‖δφ‖∞ ≤ 10⁻¹²(1+‖φ‖∞), or ‖R‖∞ ≤ 10⁻¹⁴ s (1+‖φ‖∞) with s = max diag J (the stencil scale of
   Remark 4.13); max 50 iterations; each update damped to ‖δφ‖∞ ≤ 1.
7. Collapse: Δt < 10⁻⁴Δx stops the run (§5.6). The last step is shortened to hit T (and snapshot times).
8. Random one-step states (§5.1) are loosely specified: 48 cells, log ρ ~ N(0, σ²), σ ∈ (0,3], optionally smoothed by
   a 3-point average; u = A·N(0,1), A ∈ (0,3]; ε log-uniform in [1e-4,10]; θ ∈ {0.2,0.5,0.8};
   κ = 10^U · (Δx/2)(max|u|+1), U ~ U(−1,1); Δt halved from 1 until (CFLθ) holds.
9. Errors: discrete norms Δx Σ|e_i| (L¹), (Δx Σ e_i²)^{½} (L²), max|e_i| (L∞) at cell centres; fine-grid
   references are restricted by averaging blocks of N_ref/N cells.
