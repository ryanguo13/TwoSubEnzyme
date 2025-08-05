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

| Item | Single-Substrate Reaction | Two-Substrate Reaction |
|------|---------------------------|------------------------|
| Reaction | S + E <=> ES <=> E + P | S + T + E <=> EST <=> E + P + Q |
| Rate Equation | v = (Vmf*[S]/Ks - Vmr*[P]/Kp) / (1 + [S]/Ks + [P]/Kp) | v = (Vmf*[S][T]/(Ks*Kt) - Vmr*[P][Q]/(Kp*Kq)) / (1 + [S][T]/(Ks*Kt) + [P][Q]/(Kp*Kq)) |
| Reaction Quotient | Q = [P]/[S] | Q = [P][Q]/([S][T]) |

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

3. Run detailed derivation:
```julia
include("derivation.jl")
```

4. Run numerical simulation:
```julia
include("simulation.jl")
```

## Output Files

- `flux_ratio_plot.png`: Plot of flux ratio vs Gibbs free energy change
- `simulation_results.png`: Numerical simulation results plot

## Key Findings

1. **Two-Substrate Effect**: The substrate term changes from `[S]/Ks` to `[S][T]/(Ks*Kt)`, reflecting the cooperative effect of two substrates.

2. **Product Inhibition**: The product term changes from `[P]/Kp` to `[P][Q]/(Kp*Kq)`, indicating the inhibitory effect of two products on the reaction.

3. **Thermodynamic Consistency**: The relationship between flux ratio and Gibbs free energy change remains unchanged, consistent with thermodynamic principles.

4. **Steady-State Assumption**: When enzyme concentration is much smaller than substrate concentration, the steady-state assumption holds, and the derivation results are highly consistent with numerical simulation.

## Applications

This derivation result can be applied to:
- Two-substrate enzyme reaction kinetics analysis
- Metabolic network modeling
- Bioreactor design
- Enzyme inhibitor research

## References

This work references the derivation method for single-substrate reversible enzyme reactions shown in the figure, extending it to the two-substrate two-product case. This method is based on:
- Michaelis-Menten kinetics
- Steady-state assumption
- Basic thermodynamic principles
- Flux balance analysis 