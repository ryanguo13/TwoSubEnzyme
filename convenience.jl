# ========================= convenience.jl =========================
# General convenience kinetics (Liebermeister & Klipp 2006)
# Applicable to: arbitrary stoichiometric reversible enzyme reactions
# ==================================================================

using ModelingToolkit, Symbolics, Latexify, Plots
import Symbolics: value


export convenience_rate, h_act, h_inh

"""
    convenience_rate(species_left, species_right, stoich_left, stoich_right,
                     Km_left, Km_right, kcat_plus, kcat_minus, E_tot)

Returns the symbolic rate expression for convenience kinetics.

Parameters
----------
- species_left / right : tuple/vector, species symbols Unit is mmol/L
- stoich_left / right  : corresponding stoichiometric coefficients (integers ≥ 0, dimensionless)
- Km_left / right      : corresponding half-saturation constants Unit is mmol/L
- kcat_plus / minus    : macroscopic forward/reverse catalytic constants Unit is 1/s
- E_tot                : total enzyme concentration Unit is mmol/L

Example
-------  
v = convenience_rate((S,T), (P,Q), (1,1), (1,1),
                     (Ks,Kt), (Kp,Kq), k⁺, k⁻, Etot) # Unit is mmol/L/s
"""
function convenience_rate(species_left, species_right,
                          stoich_left, stoich_right,
                          Km_left, Km_right,
                          kcat_plus, kcat_minus, E_tot)

    # Normalized concentrations
    x_left  = [s / k for (s, k) in zip(species_left,  Km_left)]
    x_right = [s / k for (s, k) in zip(species_right, Km_right)]

    # Numerator
    numerator = kcat_plus * prod(x_left  .^ stoich_left) -
                kcat_minus * prod(x_right .^ stoich_right)

    # Denominator: polynomial (1 + x + ... + x^α)
    poly_left  = [sum([xi^k for k in 0:α]) for (xi, α) in zip(x_left,  stoich_left)]
    poly_right = [sum([xi^l for l in 0:β]) for (xi, β) in zip(x_right, stoich_right)]

    denominator = prod(poly_left) * prod(poly_right) - 1

    return E_tot * numerator / denominator
end

# ---------- Optional: activation / inhibition factors ----------
h_act(d, Ka) = 1 + d / Ka
h_inh(d, Ki) = Ki / (Ki + d)

# ==================================================================
# Examples and visualization (only executed when running this file directly)
# ==================================================================
if abspath(PROGRAM_FILE) == @__FILE__
    # Example 1: single substrate ↔ single product  A ⇌ B
    @variables t A(t) B(t) Ka Kb k⁺ k⁻ Etot
    v1 = convenience_rate((A,), (B,), (1,), (1,), (Ka,), (Kb,), k⁺, k⁻, Etot)

    # Example 2: two substrates ↔ two products  S + T ⇌ P + Q
    @variables S(t) T(t) P(t) Q(t) Ks Kt Kp Kq
    v2 = convenience_rate((S,T), (P,Q), (1,1), (1,1),
                          (Ks,Kt), (Kp,Kq), k⁺, k⁻, Etot)

    # Example 3: 2A + B ⇌ 3C
    @variables A(t) B(t) C(t) Ka Kb Kc
    v3 = convenience_rate((A,B), (C,), (2,1), (3,),
                          (Ka,Kb), (Kc,), k⁺, k⁻, Etot)

    # One-click LaTeX output
    println("Example 2 rate expression:")
    println(latexify(v2))

    # Example: two-substrate reaction with activation/inhibition
    @variables Act(t) Inh(t) Ka_act Ki_inh
    v2_mod = v2 * h_act(Act, Ka_act) * h_inh(Inh, Ki_inh)
    println("\nWith activation/inhibition:")
    println(latexify(v2_mod))

    # Set numerical ranges for variables
    S_vals = range(0, 10, length=100)
    T_vals = range(0, 10, length=100)

    # Fixed parameters
    P_val = 0.0
    Q_val = 0.0
    Ks_val = 1.0
    Kt_val = 1.0
    Kp_val = 1.0
    Kq_val = 1.0
    k⁺_val = 1.0
    k⁻_val = 1.0
    Etot_val = 1.0

    # Build grid
    S_grid = repeat(S_vals', length(T_vals), 1)
    T_grid = repeat(T_vals, 1, length(S_vals))

    # Calculate v2 values on the grid
    v2_vals = [Symbolics.value(substitute(v2, Dict(
        S=>Sg, T=>Tg, P=>P_val, Q=>Q_val,
        Ks=>Ks_val, Kt=>Kt_val, Kp=>Kp_val, Kq=>Kq_val,
        k⁺=>k⁺_val, k⁻=>k⁻_val, Etot=>Etot_val
    ))) for (Sg, Tg) in zip(S_grid[:], T_grid[:])]

    v2_mat = reshape(v2_vals, length(T_vals), length(S_vals))

    # Plot 3D surface
    surface(S_vals, T_vals, v2_mat,
        xlabel="S", ylabel="T", zlabel="v₂",
        title="Convenience rate v₂ as a function of S and T",
        color=:viridis,
        size=(800, 600)
    )
    savefig("convenience_rate_3d.svg")
end