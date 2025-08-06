# Two-Substrate Two-Product Reversible Enzyme Reaction Derivation

This project uses Julia Catalyst package to derive the kinetic equations for the two-substrate two-product reversible enzyme reaction `S + T + E <=> EST <=> E + P + Q`, referencing the derivation method for single-substrate single-product reactions shown in the figure.

## Project Structure

- `main.jl`: Main derivation code, including complete mathematical derivation and visualization
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

Under steady-state conditions, $\frac{d[EST]}{dt} = 0$, we obtain:

$$[EST] = \frac{(k_{1f}[S][T] + k_{2r}[P][Q])[E]_{tot}}{k_{1r} + k_{2f} + k_{1f}[S][T] + k_{2r}[P][Q]}$$

### 2. Net Reaction Rate Equation

$$v = \frac{V_{mf}[S][T]/(K_s K_t) - V_{mr}[P][Q]/(K_p K_q)}{1 + [S][T]/(K_s K_t) + [P][Q]/(K_p K_q)}$$

Where:

- $V_{mf} = k_{2f} \cdot [E]_{tot}$: Maximum forward velocity
- $V_{mr} = k_{1r} \cdot [E]_{tot}$: Maximum reverse velocity
- $K_s = \frac{k_{1r} + k_{2f}}{k_{1f}}$: Michaelis constant for substrate S
- $K_t = \frac{k_{1r} + k_{2f}}{k_{1f}}$: Michaelis constant for substrate T
- $K_p = \frac{k_{1r} + k_{2f}}{k_{2r}}$: Michaelis constant for product P
- $K_q = \frac{k_{1r} + k_{2f}}{k_{2r}}$: Michaelis constant for product Q

### 3. Thermodynamic Relationships

- Reaction quotient: $Q = \frac{[P][Q]}{[S][T]}$
- Equilibrium constant: $K_{eq} = \frac{k_{1f}k_{2f}}{k_{1r}k_{2r}}$
- Standard Gibbs free energy change: $\Delta G^\circ = -RT\ln(K_{eq})$
- Actual Gibbs free energy change: $\Delta G_r = \Delta G^\circ + RT\ln(Q)$

### 4. Flux Ratio

$$\frac{J_+}{J_-} = e^{-\Delta G_r/RT} = \frac{K_{eq}}{Q}$$

## Comparison with Single-Substrate Reaction

| Item              | Single-Substrate Reaction                             | Two-Substrate Reaction                                                                    |
| ----------------- | ----------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Reaction          | $S + E \rightleftharpoons ES \rightleftharpoons E + P$ | $S + T + E \rightleftharpoons EST \rightleftharpoons E + P + Q$                           |
| Rate Equation     | $v = \frac{V_{mf}[S]/K_s - V_{mr}[P]/K_p}{1 + [S]/K_s + [P]/K_p}$ | $v = \frac{V_{mf}[S][T]/(K_s K_t) - V_{mr}[P][Q]/(K_p K_q)}{1 + [S][T]/(K_s K_t) + [P][Q]/(K_p K_q)}$ |
| Reaction Quotient | $Q = \frac{[P]}{[S]}$                                 | $Q = \frac{[P][Q]}{[S][T]}$                                                               |

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
