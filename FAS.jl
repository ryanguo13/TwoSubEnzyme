# FAS.jl — Acetoacetyl‑CoA + NADPH + H+ ⇌ (S)-3‑Hydroxybutyryl‑CoA + NADP+
# Modeling using convenience kinetics with thermodynamic constraints via ΔG°'/Keq

using ModelingToolkit, Symbolics, Latexify
import Symbolics: value

include("convenience.jl")

export build_fas_convenience

"""
    build_fas_convenience(; T = 298.15, ΔG0p_kJ = -16.8, Keq_override = nothing)

Builds the convenience kinetics rate expression for fatty acid synthesis step (reductase)
acetoacetyl‑CoA + NADPH + H+ ⇌ (S)-3‑hydroxybutyryl‑CoA + NADP+

- Calculates equilibrium constant from ΔG°' (in kJ/mol): `Keq = exp(-ΔG°'/(R*T))`,
  where `R = 8.314462618e-3 kJ/(mol·K)`; can also override with `Keq_override` (e.g., 2.7e2).
- Uses Haldane relationship: `Keq = (kcat_plus/kcat_minus) * Π(Km_right^β)/Π(Km_left^α)`,
  to derive `kcat_minus` from `kcat_plus` and Km values.

Returns a named tuple containing:
- `v` : symbolic rate expression (unit: mmol/L/s)
- `vars` : species variables tuple `(AcacCoA, NADPH, H, HB3CoA, NADP)` (unit: mmol/L)
- `params` : parameter dictionary (containing Keq, ΔG0p_kJ, etc.)
- `symbols` : related parameter symbols (Km, kcat, Etot, etc.)

Unit specifications:
- All concentrations: mmol/L
- Rate: mmol/L/s
- kcat: 1/s
- Km: mmol/L
- ΔG°': kJ/mol
- T: K
"""
function build_fas_convenience(; T = 298.15, ΔG0p_kJ = -16.8, Keq_override = nothing)
    # Constants and equilibrium constant
    R_kJ = 8.314462618e-3  # kJ/(mol·K)
    Keq_val = isnothing(Keq_override) ? exp(-ΔG0p_kJ / (R_kJ * T)) : Keq_override

    # Species (time-dependent symbols for ODE/DAE system integration)
    @variables t
    @variables AcacCoA(t) NADPH(t) H(t) HB3CoA(t) NADP(t)

    # Stoichiometric coefficients: AcacCoA + NADPH + H+ ⇌ HB3CoA + NADP+
    stoich_left  = (1, 1, 1)
    stoich_right = (1, 1)

    # Km, kcat, total enzyme (symbols)
    @parameters K_AcacCoA K_NADPH K_H K_HB3CoA K_NADP
    @parameters kcat_plus kcat_minus Etot

    # Calculate kcat_minus using Haldane relationship
    K_prod_left  = K_AcacCoA^stoich_left[1] * K_NADPH^stoich_left[2] * K_H^stoich_left[3]
    K_prod_right = K_HB3CoA^stoich_right[1] * K_NADP^stoich_right[2]
    kcat_minus_expr = kcat_plus * K_prod_right / (Keq_val * K_prod_left)

    # Use convenience_rate function directly from convenience.jl
    v = convenience_rate((AcacCoA, NADPH, H), (HB3CoA, NADP),
                         stoich_left, stoich_right,
                         (K_AcacCoA, K_NADPH, K_H), (K_HB3CoA, K_NADP),
                         kcat_plus, kcat_minus_expr, Etot)

    vars = (AcacCoA, NADPH, H, HB3CoA, NADP)
    symbols = (K_AcacCoA=K_AcacCoA, K_NADPH=K_NADPH, K_H=K_H,
               K_HB3CoA=K_HB3CoA, K_NADP=K_NADP,
               kcat_plus=kcat_plus, kcat_minus=kcat_minus, Etot=Etot)

    params = Dict(
        :R_kJ => R_kJ,
        :T => T,
        :ΔG0p_kJ => ΔG0p_kJ,
        :Keq => Keq_val
    )

    return (v=v, vars=vars, params=params, symbols=symbols)
end

println("=== FAS Convenience Kinetics Test ===")
    
# Build model
result = build_fas_convenience()

println("Reaction: Acetoacetyl-CoA + NADPH + H+ ⇌ (S)-3-Hydroxybutyryl-CoA + NADP+")
println("ΔG°' = $(result.params[:ΔG0p_kJ]) kJ/mol")
println("Keq = $(result.params[:Keq])")
println("T = $(result.params[:T]) K")

println("\nRate expression:")
println(latexify(result.v))

# Numerical test
println("\n=== Numerical Test ===")

# Set kinetic parameter values (Km, kcat, enzyme concentration)
param_vals = Dict(
    result.symbols.K_AcacCoA => 0.01,  # mmol/L - Km for AcacCoA
    result.symbols.K_NADPH => 0.05,   # mmol/L - Km for NADPH
    result.symbols.K_H => 1e-4,       # mmol/L - Km for H+ (pH 7)
    result.symbols.K_HB3CoA => 0.01,   # mmol/L - Km for HB3CoA
    result.symbols.K_NADP => 0.05,     # mmol/L - Km for NADP+
    result.symbols.kcat_plus => 10,  # 1/s - forward catalytic constant
    result.symbols.Etot => 1e-3       # mmol/L - total enzyme concentration
)

# Set reaction species concentrations (substrates and products)
conc_vals = Dict(
    result.vars[1] => 1.0,    # mmol/L - AcacCoA concentration
    result.vars[2] => 0.5,    # mmol/L - NADPH concentration
    result.vars[3] => 1e-4,   # mmol/L - H+ concentration (pH 7.0: 1e-4 mmol/L = 1e-7 mol/L)
    result.vars[4] => 0.1,    # mmol/L - HB3CoA concentration
    result.vars[5] => 0.2     # mmol/L - NADP+ concentration
)

# Calculate rate
v_num = substitute(result.v, merge(param_vals, conc_vals))
v_val = Symbolics.value(v_num)

println("Reaction rate under given conditions: $(v_val) mmol/L/s")

# Verify thermodynamic consistency - corrected reaction quotient calculation, excluding H+
Q = (conc_vals[result.vars[4]] * conc_vals[result.vars[5]]) / 
    (conc_vals[result.vars[1]] * conc_vals[result.vars[2]])
ΔG_actual = result.params[:ΔG0p_kJ] + result.params[:R_kJ] * result.params[:T] * log(Q)

println("Reaction quotient Q = $Q ")
println("Actual ΔG = $(ΔG_actual) kJ/mol")
println("ΔG°' = $(result.params[:ΔG0p_kJ]) kJ/mol")
println("v = $(v_val) mmol/L/s")
if v_val > 0
    println("Reaction proceeds forward (consuming AcacCoA + NADPH)")
elseif v_val < 0
    println("Reaction proceeds backward (consuming HB3CoA + NADP+)")
else
    println("Reaction is at equilibrium")
end

println("\nLaTeX format rate expression:")
println(latexify(result.v))

println("\n=== Example completed ===") 