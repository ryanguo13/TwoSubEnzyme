# Two-Substrate Two-Product Reversible Enzyme Reaction Derivation

This project uses Julia Catalyst package to derive the kinetic equations for the two-substrate two-product reversible enzyme reaction `S + T + E <=> EST <=> E + P + Q`, referencing the derivation method for single-substrate single-product reactions shown in the figure.

## Project Structure

- `main.jl`: Main derivation code, including complete mathematical derivation and visualization
- `simulation.jl`: Numerical simulation to validate derivation results
- `convenience.jl`: General convenience kinetics implementation (Liebermeister & Klipp 2006)
- `FAS.jl`: Fatty acid synthesis reaction modeling using convenience kinetics
- `test_fas.jl`: Test script for FAS.jl implementation
- `FAS_usage_example.jl`: Usage examples for FAS.jl
- `Project.toml`: Project dependency file

## Convenience Kinetics Implementation

### `convenience.jl`

A general implementation of convenience kinetics for arbitrary stoichiometric reversible enzyme reactions:

```julia
v = convenience_rate(species_left, species_right, stoich_left, stoich_right,
                     Km_left, Km_right, kcat_plus, kcat_minus, E_tot)
```

**Features:**

- Supports arbitrary stoichiometries (e.g., 2A + B ⇌ 3C)
- Includes activation/inhibition factors
- Provides LaTeX output for mathematical expressions
- Includes 3D visualization examples

### `FAS.jl`

Specific implementation for the fatty acid synthesis reaction:

```
Acetoacetyl-CoA + NADPH + H+ ⇌ (S)-3-Hydroxybutyryl-CoA + NADP+
```

**Key Features:**

- Thermodynamically consistent with ΔG°' = -16.8 kJ/mol
- Enforces Haldane relationship: Keq = (kcat_plus/kcat_minus) × Π(Km_right^β)/Π(Km_left^α)
- Supports custom temperature and equilibrium constants
- Returns symbolic rate expression ready for ODE/DAE integration

**Usage:**

```julia
# Basic usage with default parameters
result = build_fas_convenience()

# Use specific Keq value from literature
result = build_fas_convenience(Keq_override = 2.7e2)

# Custom temperature (physiological)
result = build_fas_convenience(T = 310.0)

# Access components
v = result.v                    # Symbolic rate expression
vars = result.vars              # Species variables
params = result.params          # Thermodynamic parameters
symbols = result.symbols        # Kinetic parameters
```

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

$$
[EST] = \frac{(k_{1f}[S][T] + k_{2r}[P][Q])[E]_{tot}}{k_{1r} + k_{2f} + k_{1f}[S][T] + k_{2r}[P][Q]}
$$

### 2. Net Reaction Rate Equation

$$
v = \frac{V_{mf}[S][T]/(K_s K_t) - V_{mr}[P][Q]/(K_p K_q)}{1 + [S][T]/(K_s K_t) + [P][Q]/(K_p K_q)}
$$

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

$$
\frac{J_+}{J_-} = e^{-\Delta G_r/RT} = \frac{K_{eq}}{Q}
$$

## Comparison with Single-Substrate Reaction

| Item              | Single-Substrate Reaction                                           | Two-Substrate Reaction                                                                                  |
| ----------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------- |
| Reaction          | $S + E \rightleftharpoons ES \rightleftharpoons E + P$            | $S + T + E \rightleftharpoons EST \rightleftharpoons E + P + Q$                                       |
| Rate Equation     | $v = \frac{V_{mf}[S]/K_s - V_{mr}[P]/K_p}{1 + [S]/K_s + [P]/K_p}$ | $v = \frac{V_{mf}[S][T]/(K_s K_t) - V_{mr}[P][Q]/(K_p K_q)}{1 + [S][T]/(K_s K_t) + [P][Q]/(K_p K_q)}$ |
| Reaction Quotient | $Q = \frac{[P]}{[S]}$                                             | $Q = \frac{[P][Q]}{[S][T]}$                                                                           |

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

4. Test FAS convenience kinetics:

```julia
include("FAS.jl")
```

## FAS Reaction Parameters

The FAS reaction is modeled with the following thermodynamic parameters:

- **Standard Gibbs Free Energy**: ΔG°' = -16.8 kJ/mol
- **Equilibrium Constant**: Keq = 8.7 × 10² (from literature)
- **Temperature**: T = 298.15 K (default)
- **Gas Constant**: R = 8.314462618 × 10⁻³ kJ/(mol·K)

The implementation ensures thermodynamic consistency by enforcing the Haldane relationship between forward and reverse catalytic constants.

## Results

```

=== FAS Convenience Kinetics Test ===
Reaction: Acetoacetyl-CoA + NADPH + H+ ⇌ (S)-3-Hydroxybutyryl-CoA + NADP+
ΔG°' = -16.8 kJ/mol
Keq = 877.4707765445186
T = 298.15 K

Rate expression:
\begin{equation}
\frac{\mathtt{Etot} \left( \frac{ - \mathtt{kcat\_plus} \mathtt{NADP}\left( t \right) \mathtt{HB3CoA}\left( t \right)}{877.47 \mathtt{K\_AcacCoA} \mathtt{K\_H} \mathtt{K\_NADPH}} + \frac{\mathtt{kcat\_plus} H\left( t \right) \mathtt{NADPH}\left( t \right) \mathtt{AcacCoA}\left( t \right)}{\mathtt{K\_AcacCoA} \mathtt{K\_H} \mathtt{K\_NADPH}} \right)}{-1 + \left( 1 + \frac{\mathtt{HB3CoA}\left( t \right)}{\mathtt{K\_HB3CoA}} \right) \left( 1 + \frac{H\left( t \right)}{\mathtt{K\_H}} \right) \left( 1 + \frac{\mathtt{NADPH}\left( t \right)}{\mathtt{K\_NADPH}} \right) \left( 1 + \frac{\mathtt{AcacCoA}\left( t \right)}{\mathtt{K\_AcacCoA}} \right) \left( 1 + \frac{\mathtt{NADP}\left( t \right)}{\mathtt{K\_NADP}} \right)}
\end{equation}


=== Numerical Test ===
Reaction rate under given conditions: 4.4525719281013987e-5 mmol/L/s
Reaction quotient Q = 0.04000000000000001 
Actual ΔG = -24.779454853327145 kJ/mol
ΔG°' = -16.8 kJ/mol
v = 4.4525719281013987e-5 mmol/L/s
Reaction proceeds forward (consuming AcacCoA + NADPH)

LaTeX format rate expression:
\begin{equation}
\frac{\mathtt{Etot} \left( \frac{ - \mathtt{kcat\_plus} \mathtt{NADP}\left( t \right) \mathtt{HB3CoA}\left( t \right)}{877.47 \mathtt{K\_AcacCoA} \mathtt{K\_H} \mathtt{K\_NADPH}} + \frac{\mathtt{kcat\_plus} H\left( t \right) \mathtt{NADPH}\left( t \right) \mathtt{AcacCoA}\left( t \right)}{\mathtt{K\_AcacCoA} \mathtt{K\_H} \mathtt{K\_NADPH}} \right)}{-1 + \left( 1 + \frac{\mathtt{HB3CoA}\left( t \right)}{\mathtt{K\_HB3CoA}} \right) \left( 1 + \frac{H\left( t \right)}{\mathtt{K\_H}} \right) \left( 1 + \frac{\mathtt{NADPH}\left( t \right)}{\mathtt{K\_NADPH}} \right) \left( 1 + \frac{\mathtt{AcacCoA}\left( t \right)}{\mathtt{K\_AcacCoA}} \right) \left( 1 + \frac{\mathtt{NADP}\left( t \right)}{\mathtt{K\_NADP}} \right)}
\end{equation}


=== Example completed ===
```
