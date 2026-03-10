# UAV-Assisted SWIPT Cooperative Relaying with Physical Layer Security

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![IEEE Access](https://img.shields.io/badge/Paper-IEEE%20Access-orange.svg)](https://doi.org/10.1109/ACCESS.2019.2948384)

---

## 📖 Overview

This repository contains the full MATLAB simulation code accompanying the paper:

> **"Physical Layer Security of UAV-Assisted Cooperative Communication with SWIPT and Cooperative Jamming"**  
> *Published in IEEE Access*, 2019  
> DOI: [10.1109/ACCESS.2019.2948384](https://doi.org/10.1109/ACCESS.2019.2948384)

The paper investigates the **secrecy performance** of a UAV (Unmanned Aerial Vehicle) relay network that simultaneously supports **energy harvesting** and **information decoding** via a Power Splitting (PS) receiver, while using Bob as a **cooperative jammer** to confuse eavesdroppers.

---

## 🔬 Research Summary

### System Model


- **Alice (A):** Source node (ground), transmits information with power `λP`.
- **UAV Relay (U):** Flying relay, uses a **Power Splitting Receiver (PSR)** with ratio `β`:
  - Fraction `β` of received power → energy harvesting (supplies relay transmit power).
  - Fraction `(1–β)` → information decoding (decode-and-forward relay).
- **Bob (B):** Destination (ground), also transmits cooperative jamming power `(1–λ)P`.
- **Eve (E):** Passive eavesdropper (ground), wiretaps both the Alice→U and U→B links.

### Key Features

| Feature | Description |
|---------|-------------|
| **Channel Model** | Elevation-angle-dependent Rician fading for A2G links; Rayleigh for G2G |
| **SWIPT Protocol** | Power Splitting Receiver (PSR) at UAV |
| **Security Mechanism** | Bob-aided cooperative jamming |
| **Power Allocation** | Joint optimization of `λ` (PAF) and `β` (PSR) |
| **Benchmarks** | Equal Power Allocation (EPA), Optimal PSR, Static ground relay |

### Key Contributions

1. **Closed-form connection probability** expression for the Rician-fading SWIPT relay channel using a series expansion of the Marcum-Q function.
2. **Analytical lower bound** on the Ergodic Secrecy Rate (ESR) via Jensen's inequality approximation for non-central chi-squared random variables.
3. **Closed-form optimal power allocation** `λ*` that maximizes instantaneous secrecy rate.
4. **Asymptotic analysis** of ESR at high SNR, yielding the secrecy pre-log and offset.
5. **Flying UAV vs static relay comparison** demonstrating superior secrecy performance with trajectory optimization.

---

## 📁 Repository Structure

```
UAV-PLS-SWIPT-Relay/
│
├── README.md                   ← This file
├── LICENSE                     ← MIT License
├── .gitignore
│
├── src/                        ← Core reusable functions
│   ├── analytics/              ← Mathematical helper functions
│   │   ├── g1.m                ← E[log X] approximation (ncx2 RV)
│   │   ├── g2.m                ← E[log(X+b)] approximation
│   │   ├── g3.m                ← E[1/X] approximation
│   │   ├── PHI.m               ← PHI(r,b) auxiliary function
│   │   ├── PSI.m               ← PSI(m,r,β) auxiliary function
│   │   ├── Gamma_series.m      ← Gamma series helper
│   │   ├── I1.m                ← I1(t,β) integral helper
│   │   ├── Ei.m                ← Exponential integral Ei(x)
│   │   └── coeff.m             ← Series coefficient helper
│   ├── channel/
│   │   └── channel_utils.m     ← A2G/G2G channel generation utilities
│   └── optimization/
│       ├── obj.m               ← Secrecy rate objective (for λ optimization)
│       └── objfunc.m           ← Polynomial ratio objective (for β optimization)
│
├── simulations/                ← Paper figure reproduction scripts
│   ├── connection_probability/
│   │   ├── run_CP_vs_Power.m   ← Fig: CP vs transmit power (EPA)
│   │   ├── run_CP_vs_PSR.m     ← Fig: CP vs power splitting ratio
│   │   └── run_CP_OPA_vs_EPA.m ← Fig: CP — OPA vs EPA comparison
│   ├── intercept_probability/
│   │   └── run_IP_vs_Power.m   ← Fig: Intercept probability vs power
│   ├── secrecy_rate/
│   │   ├── run_ESR_vs_SNR.m    ← Fig: ESR vs SNR (analytical LB vs simulation)
│   │   └── run_ESR_vs_PAF.m    ← Fig: ESR vs power allocation factor
│   ├── secrecy_outage/
│   │   └── run_SOP_vs_Power.m  ← Fig: SOP with PPP eavesdroppers
│   ├── uav_trajectory/
│   │   └── run_FlyingUAV_vs_Static.m ← Fig: Flying UAV vs static relay
│   └── comprehensive/
│       ├── run_All_Metrics_vs_PSR.m  ← Fig: All metrics vs PSR (subplot)
│       └── run_Comprehensive_Simulation.m ← All main paper figures
│
├── tests/                      ← Debug and validation scripts
│   ├── test_g1_g2_accuracy.m
│   ├── test_marcumq_approximation.m
│   └── test_channel_normalization.m
│
├── examples/
│   └── quick_start.m           ← Minimal working example (~30 sec)
│
├── data/
│   └── processed/              ← Generated .mat result files (see below)
│
├── results/
│   └── figures/                ← Generated figures (.fig and .png)
│
└── docs/
    └── parameter_description.md ← Full parameter reference
```

---

## 🗂️ Missing `.mat` Files

Some numerical results were pre-computed and saved as `.mat` files, but **are not included** in this repository (they are generated by running the scripts). The table below documents each expected file:

| File | Generated By | Contents |
|------|-------------|----------|
| `data/processed/CP_vs_Power.mat` | `run_CP_vs_Power.m` | `P_dBW`, `Pc_Simulation`, `Pc_Exact`, `Pc_Approximate` |
| `data/processed/ESR_vs_SNR.mat` | `run_ESR_vs_SNR.m` | `P_dBW`, `ESR_sim`, `ESR_analytical` |
| `data/processed/ESR_vs_PSR.mat` | `run_All_Metrics_vs_PSR.m` | All metrics vs β |
| `data/processed/Comprehensive_Results.mat` | `run_Comprehensive_Simulation.m` | All main figures data |

**To regenerate:** Simply run the corresponding script. Each simulation script saves its output automatically to `data/processed/`.

---

## ⚙️ Dependencies & Requirements

| Requirement | Version | Notes |
|-------------|---------|-------|
| MATLAB | R2019b or later | Core engine |
| Statistics & Machine Learning Toolbox | Any | `ncx2rnd`, `marcumq`, `chi2rnd` |
| Signal Processing Toolbox | Any | `db2pow`, `pow2db` |
| Symbolic Math Toolbox | Any | Used only in `Ultimate_Simulation.m` for PSR optimization via `optimproblem`; optional for other scripts |

**Tested on:** MATLAB R2021a, R2022b, R2023a.

---

## 🚀 Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/UAV-PLS-SWIPT-Relay.git
cd UAV-PLS-SWIPT-Relay
```

### 2. Run the quick-start example

```matlab
cd examples
quick_start
```

This generates an ESR vs. transmit power plot in under 30 seconds.

### 3. Reproduce all paper figures

```matlab
cd simulations/comprehensive
run_Comprehensive_Simulation   % ~20-40 minutes
```

Or run individual figures:

```matlab
cd simulations/connection_probability
run_CP_vs_Power               % Connection Probability vs Power (~3 min)

cd ../secrecy_rate
run_ESR_vs_SNR                % ESR analytical vs simulation (~5 min)

cd ../uav_trajectory
run_FlyingUAV_vs_Static       % Flying UAV vs static relay (~15 min)
```

---

## 📊 How to Reproduce Each Figure

| Paper Figure | Script | Output File |
|---|---|---|
| CP vs Transmit Power | `simulations/connection_probability/run_CP_vs_Power.m` | `results/figures/CP_vs_Power.png` |
| CP vs PSR (β) | `simulations/connection_probability/run_CP_vs_PSR.m` | `results/figures/CP_vs_PSR.png` |
| CP: OPA vs EPA | `simulations/connection_probability/run_CP_OPA_vs_EPA.m` | `results/figures/CP_OPA_vs_EPA.png` |
| IP vs Transmit Power | `simulations/intercept_probability/run_IP_vs_Power.m` | `results/figures/IP_vs_Power.png` |
| ESR vs SNR | `simulations/secrecy_rate/run_ESR_vs_SNR.m` | `results/figures/ESR_vs_SNR.png` |
| ESR vs PAF (λ) | `simulations/secrecy_rate/run_ESR_vs_PAF.m` | `results/figures/ESR_vs_PAF.png` |
| SOP vs Power | `simulations/secrecy_outage/run_SOP_vs_Power.m` | `results/figures/SOP_vs_Power.png` |
| All Metrics vs PSR | `simulations/comprehensive/run_All_Metrics_vs_PSR.m` | `results/figures/All_Metrics_vs_PSR.png` |
| Flying UAV vs Static | `simulations/uav_trajectory/run_FlyingUAV_vs_Static.m` | `results/figures/FlyingUAV_vs_StaticRelay.png` |

---

## 🔑 Key Parameters

| Symbol | Variable | Description | Typical Value |
|--------|----------|-------------|---------------|
| `β` | `beta` | Power Splitting Ratio at UAV | 0.5 |
| `λ` | `lambda` | Power Allocation Factor (Alice share) | 0.5 – 0.8 |
| `η` | `eff` | Energy Harvesting Efficiency | 0.7 |
| `ζ` | `zeta` | Pilot-to-Noise Power Ratio Np/N0 | 2 |
| `αL` | `alpha_L` | LoS Path-Loss Exponent | 2.0 |
| `αN` | `alpha_N` | NLoS Path-Loss Exponent | 3.5 |
| `km, kM` | `km, kM` | Rician K-factor range [dB] | 1, 10 |
| `Rt` | `Rt` | Target Transmission Rate [bits/s/Hz] | 0.5 – 2 |
| `Rs` | `Rs` | Target Secrecy Rate [bits/s/Hz] | 0.2 – 1 |
| `N0` | `N0` | Noise Power (linear) | 1e-2 |

---

## 📐 Mathematical Background

The paper derives:

**1. Connection Probability (closed-form):**
$$P_c = (1+K_B)\exp(-K_A-K_B-(1+K_A)A) \cdot \mathcal{F}(K_A, K_B, A, B)$$

where $A$ and $B$ are SNR-threshold-dependent constants and $\mathcal{F}$ is a double series involving modified Bessel functions.

**2. ESR Lower Bound:**
$$\bar{R}_s \geq \frac{1}{2\log 2}\left[\log(1+e^{T_1}) - \log(1+T_2)\right]^+$$

where $T_1 \approx \mathbb{E}[\log \Gamma_{AB}]$ uses the `g1()` and `g2()` functions, and $T_2 \approx \mathbb{E}[\Gamma_E]$ uses the mean-based approximation.

**3. Optimal Power Allocation:**
$$\lambda^* = \frac{-WY + \sqrt{WY(VY + WX)}}{(X-Y)W + VY}$$

---

## 📝 Citation

If you use this code in your research, please cite:

```bibtex
@article{mamaghani2019performance,
  title={On the performance of low-altitude UAV-enabled secure AF relaying with cooperative jamming and SWIPT},
  author={Tatar Mamaghani, Milad and Hong, Yi},
  journal={IEEE access},
  volume={7},
  pages={153060--153073},
  year={2019},
}
```


---

## 📜 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

