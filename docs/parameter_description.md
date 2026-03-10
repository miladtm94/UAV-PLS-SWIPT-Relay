# System Parameter Reference

This document provides a complete reference for all system parameters used in the
simulation code for the paper on UAV-assisted SWIPT cooperative relaying with
physical layer security.

---

## Network Topology

| Symbol | Variable | Description | Units |
|--------|----------|-------------|-------|
| D | `D` | Total Alice-to-Bob horizontal distance | m (or normalized) |
| H | `H` | UAV flight altitude | m |
| t_A | `t_A` | Alice 3D position vector `[x, y, z]` | m |
| t_U | `t_U` | UAV 3D position vector `[x, y, z]` | m |
| t_B | `t_B` | Bob 3D position vector `[x, y, z]` | m |
| t_E | `t_E` | Eve 3D position vector `[x, y, z]` | m |
| d_au | `d_au` | 3D distance Alice → UAV = `norm(t_A - t_U)` | m |
| d_ub | `d_ub` | 3D distance UAV → Bob = `norm(t_U - t_B)` | m |
| d_ue | `d_ue` | 3D distance UAV → Eve = `norm(t_U - t_E)` | m |
| d_ae | `d_ae` | Distance Alice → Eve (ground link) | m |
| d_be | `d_be` | Distance Bob → Eve (ground link) | m |

---

## Channel Model Parameters

### Path Loss (Air-to-Ground, A2G)

| Symbol | Variable | Description | Typical Value |
|--------|----------|-------------|---------------|
| α_L | `alpha_L` | LoS link path-loss exponent | 2.0 |
| α_N | `alpha_N` | NLoS link path-loss exponent | 3.5 |
| a | `a` | LoS probability sigmoid parameter (urban) | 0.28 |
| b | `b` | LoS probability sigmoid parameter (urban) | 9.61 |
| θ | `theta_xy` | Elevation angle from ground to UAV (radians) | `asin(H/d_3D)` |
| P_LoS | `P_los_xy` | Probability of LoS link | `1/(1+a*exp(-b*(θ-a)))` |
| α_eff | `alpha_xy` | Effective path-loss exponent (weighted) | `(αL−αN)*PLoS + αN` |
| L_xy | `L_xy` | Large-scale path loss (linear) | `d_3D^(−α_eff)` |
| K_xy | `K_xy` | Rician K-factor (linear, elevation-dependent) | `db2pow(km+(kM−km)*2θ/π)` |
| k_m | `km` | Minimum K-factor in dB (near-horizontal) | 1 dB |
| k_M | `kM` | Maximum K-factor in dB (near-vertical) | 10 dB |

### Path Loss (Ground-to-Ground, G2G — Rayleigh)

| Symbol | Variable | Description | Typical Value |
|--------|----------|-------------|---------------|
| α_G | `alpha_N` | G2G path-loss exponent | 3.5 |
| L_ae | `L_ae` | Alice→Eve path loss = `d_ae^(-αG)` | computed |
| L_be | `L_be` | Bob→Eve path loss = `d_be^(-αG)` | computed |

### Small-Scale Fading

| Link | Distribution | Parameter | Variable |
|------|-------------|-----------|----------|
| A→U (Alice→UAV) | Rician (ncx2) | K_au | `K_au`, `lambda_au = sqrt(2*K_au)` |
| U→B (UAV→Bob) | Rician (ncx2) | K_ub | `K_ub`, `lambda_ub = sqrt(2*K_ub)` |
| U→E (UAV→Eve) | Rician (ncx2) | K_ue | `K_ue`, `lambda_ue = sqrt(2*K_ue)` |
| A→E (Alice→Eve) | Rayleigh (exp) | L_ae | `V ~ L_ae * exprnd(1)` |
| B→E (Bob→Eve) | Rayleigh (exp) | L_be | `W ~ L_be * exprnd(1)` |

**Channel normalization:** Rician channels are generated as:
```
X = ncx2rnd(2, sqrt(2*K_au), [1,N])   % raw
X = X / mean(X) * L_au                % normalized so E[X] = L_au
```

---

## SWIPT Parameters

| Symbol | Variable | Description | Range |
|--------|----------|-------------|-------|
| β | `beta` | Power Splitting Ratio at UAV relay | (0, 1) |
| η | `eff` | RF-to-DC energy harvesting efficiency | (0, 1), typical: 0.7 |
| ζ | `zeta` | Pilot-to-noise power ratio N_p/N_0 | ≥ 0, typical: 2 |
| ε | `epsilon` | Effective noise due to channel estimation error | computed |

The effective noise:
```
epsilon = zeta * N0^2 ./ (Pa*X + Pb*Y)
```

---

## Power Allocation Parameters

| Symbol | Variable | Description | Range |
|--------|----------|-------------|-------|
| P | `P` | Total network transmit power (linear) | W |
| λ | `lambda` | Power Allocation Factor (Alice's share) | [0, 1] |
| P_a | `Pa` | Alice transmit power = λ·P | W |
| P_b | `Pb` | Bob jamming power = (1–λ)·P | W |
| λ* | `lambda_opt` | Optimal PAF (closed-form) | see below |

**Optimal PAF:**
```
lambda_opt = (-W*Y + sqrt(W*Y*(V*Y + W*X))) / ((X-Y)*W + V*Y)
```
(instantaneous; averaged over channel realizations)

---

## SNR Expressions

**Desired link SNR at Bob (after SWIPT DF relay):**
```
Gamma_AB = (eff*beta*(1-beta)*Pa.*X.*Y) / 
           (eff*beta*(1-beta+zeta)*Y*N0 + (1-beta)*N0 + epsilon)
```

**Eavesdropper SNR — direct path (Alice→Eve):**
```
Gamma_e1 = (Pa/N0)*V / ((Pb/N0)*W + 1)
```

**Eavesdropper SNR — relay path (UAV→Eve):**
```
Gamma_e2 = eff*beta*(1-beta)*Pa*X*Z / 
           ((eff*beta*(1-beta)*Pb*Y*Z + Z*eff*beta*(1-beta+zeta)*N0 + (1-beta)*N0 + epsilon))
```

**Total eavesdropper SNR:**
- **Worst case (max combiner):** `Gamma_E = max(Gamma_e1, Gamma_e2)`
- **MRC bound:** `Gamma_E = Gamma_e1 + Gamma_e2`

---

## Performance Metrics

| Metric | Symbol | Definition |
|--------|--------|-----------|
| Connection Probability | CP | `Pr{0.5*log2(1+Gamma_AB) > Rt}` |
| Intercept Probability | IP | `Pr{0.5*log2(1+Gamma_E) > Re}` |
| Ergodic Secrecy Rate | ESR | `E[max(Ct - Ce, 0)]` where Ct = 0.5*log2(1+Gamma_AB), Ce = 0.5*log2(1+Gamma_E) |
| Secrecy Outage Probability | SOP | `Pr{Rs_inst < Rs_target}` |

---

## Rate Parameters

| Symbol | Variable | Description | Typical |
|--------|----------|-------------|---------|
| R_t | `Rt` | Target transmission rate | 0.5 – 2 bits/s/Hz |
| R_s | `Rs` | Target secrecy rate | 0.2 – 1 bits/s/Hz |
| R_e | `Re = Rt - Rs` | Leakage rate threshold | computed |
| δ_t | `delta_t = 2^(2*Rt)-1` | SNR threshold for CP | computed |
| δ_e | `delta_e = 2^(2*Re)-1` | SNR threshold for IP | computed |

---

## Monte Carlo Settings

| Variable | Description | Typical Value |
|----------|-------------|---------------|
| `N` | Number of channel realizations (inner loop) | 1e4 – 1e6 |
| `n` | Number of parameter sweep points | 20 – 50 |
| `Ne` | Number of PPP realizations (SOP script only) | 1e3 |

---

## Series Truncation Orders

The closed-form expressions use truncated series. Accuracy improves with order:

| Variable | Parameter | Effect | Recommended |
|----------|-----------|--------|-------------|
| `D_ord` | ncx2 (Alice channel) series order | Higher → better accuracy near origin | 10–20 |
| `R_ord` | ncx2 (Bob channel) series order | Higher → better accuracy at high SNR | 10–20 |
| `R_terms` (in g1, g2) | ESR lower bound series order | Higher → tighter bound | 5–30 |
