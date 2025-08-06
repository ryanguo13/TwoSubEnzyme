# Two-Substrate Two-Product Reversible Enzyme Reaction Derivation

This project uses Julia Catalyst package to derive the kinetic equations for the two-substrate two-product reversible enzyme reaction `S + T + E <=> EST <=> E + P + Q`, referencing the derivation method for single-substrate single-product reactions shown in the figure.

## Project Structure

- `main.jl`: Main derivation code, including complete mathematical derivation and visualization
- `derivation.jl`: Detailed mathematical derivation process
- `simulation.jl`: Numerical simulation to validate derivation results
- `Project.toml`: Project dependency file

## Reaction Mechanism

```
S + T + E <=> EST <=> E + P + Q
```

Where:

- S, T: Substrates
- E: Enzyme
- EST: Enzyme-substrate complex
- P, Q: Products

## Main Derivation Results

### 1. Steady-State Enzyme-Substrate Complex Concentration

Under steady-state conditions, `d[EST]/dt = 0`, we obtain:

```
EST = (k1f*S*T + k2r*P*Q)*Etot / (k1r + k2f + k1f*S*T + k2r*P*Q)
```

### 2. Net Reaction Rate Equation

```
v = (Vmf*S*T/(Ks*Kt) - Vmr*P*Q/(Kp*Kq)) / (1 + S*T/(Ks*Kt) + P*Q/(Kp*Kq))
```

Where:

- `Vmf = k2f * Etot`: Maximum forward velocity
- `Vmr = k1r * Etot`: Maximum reverse velocity
- `Ks = (k1r + k2f) / k1f`: Michaelis constant for substrate S
- `Kt = (k1r + k2f) / k1f`: Michaelis constant for substrate T
- `Kp = (k1r + k2f) / k2r`: Michaelis constant for product P
- `Kq = (k1r + k2f) / k2r`: Michaelis constant for product Q

### 3. Thermodynamic Relationships

- Reaction quotient: `Q = [P][Q]/([S][T])`
- Equilibrium constant: `Keq = k1f*k2f/(k1r*k2r)`
- Standard Gibbs free energy change: `ΔG° = -RT*ln(Keq)`
- Actual Gibbs free energy change: `ΔGr = ΔG° + RT*ln(Q)`

### 4. Flux Ratio

```
J+/J- = e^(-ΔGr/RT) = Keq/Q
```

## Comparison with Single-Substrate Reaction

| Item              | Single-Substrate Reaction                             | Two-Substrate Reaction                                                                    |
| ----------------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Reaction          | S + E <=> ES <=> E + P                                | S + T + E <=> EST <=> E + P + Q                                                           |
| Rate Equation     | v = (Vmf*[S]/Ks - Vmr*[P]/Kp) / (1 + [S]/Ks + [P]/Kp) | v = (Vmf*[S][T]/(Ks*Kt) - Vmr*[P][Q]/(Kp*Kq)) / (1 + [S][T]/(Ks*Kt) + [P][Q]/(Kp*Kq)) |
| Reaction Quotient | Q = [P]/[S]                                           | Q = [P][Q]/([S][T])                                                                       |

## Running the Project

1. Install dependencies:

```julia
using Pkg
Pkg.instantiate()
```

2. Run main derivation:

```julia
include("main.jl")
```

3. Run numerical simulation:

```julia
include("simulation.jl")
```
