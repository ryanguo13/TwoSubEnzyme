using Catalyst
using ModelingToolkit
using Symbolics
using Latexify
using Plots

# Define reaction system: S + T + E <=> EST <=> E + P + Q
@reaction_network TwoSubEnzyme begin
    # Step 1: S + T + E <=> EST
    (k1f, k1r), S + T + E <--> EST
    
    # Step 2: EST <=> E + P + Q  
    (k2f, k2r), EST <--> E + P + Q
end

# Temperature constants
T_val = 298  # K
R_val = 8.314  # J/(mol·K)

# Get reaction system parameters and variables
@parameters k1f k1r k2f k2r
@variables t S(t) T(t) E(t) EST(t) P(t) Q(t)

# Define total enzyme concentration
Etot = E + EST

# Steady-state condition: d[EST]/dt = 0
# From the reaction network we get:
# d[EST]/dt = k1f*S*T*E - k1r*EST - k2f*EST + k2r*P*Q*E = 0

# Solve for steady-state [EST]
# k1f*S*T*E - k1r*EST - k2f*EST + k2r*P*Q*E = 0
# k1f*S*T*E + k2r*P*Q*E = (k1r + k2f)*EST
# EST = (k1f*S*T*E + k2r*P*Q*E) / (k1r + k2f)

# Since E = Etot - EST, substitute to get:
# EST = (k1f*S*T*(Etot - EST) + k2r*P*Q*(Etot - EST)) / (k1r + k2f)
# EST = (k1f*S*T*Etot + k2r*P*Q*Etot - k1f*S*T*EST - k2r*P*Q*EST) / (k1r + k2f)
# EST*(k1r + k2f + k1f*S*T + k2r*P*Q) = (k1f*S*T + k2r*P*Q)*Etot
# EST = (k1f*S*T + k2r*P*Q)*Etot / (k1r + k2f + k1f*S*T + k2r*P*Q)

# Define Michaelis constants
Ks = (k1r + k2f) / k1f  # Michaelis constant for substrate S
Kt = (k1r + k2f) / k1f  # Michaelis constant for substrate T  
Kp = (k1r + k2f) / k2r  # Michaelis constant for product P
Kq = (k1r + k2f) / k2r  # Michaelis constant for product Q

# Define maximum velocities
Vmf = k2f * Etot  # Maximum forward velocity
Vmr = k1r * Etot  # Maximum reverse velocity

# Steady-state enzyme-substrate complex concentration
EST_ss = (k1f*S*T + k2r*P*Q)*Etot / (k1r + k2f + k1f*S*T + k2r*P*Q)

# Net reaction rate v = k2f*EST - k2r*P*Q*E
# Substitute E = Etot - EST and EST_ss
v = k2f*EST_ss - k2r*P*Q*(Etot - EST_ss)

# Simplified rate equation
# v = (Vmf*S*T/(Ks*Kt) - Vmr*P*Q/(Kp*Kq)) / (1 + S*T/(Ks*Kt) + P*Q/(Kp*Kq))

# Define reaction quotient Q = [P][Q]/([S][T])
Q = P*Q/(S*T)

# Define equilibrium constant Keq = k1f*k2f/(k1r*k2r)
Keq = k1f*k2f/(k1r*k2r)

# Standard Gibbs free energy change
ΔG° = -R_val*T_val*log(Keq)  # where R is gas constant, T is temperature

# Actual Gibbs free energy change
ΔGr = ΔG° + R_val*T_val*log(Q)

# Forward and reverse fluxes
J_plus = k2f*EST_ss  # Forward flux
J_minus = k2r*P*Q*(Etot - EST_ss)  # Reverse flux

# Flux ratio
flux_ratio = J_plus/J_minus

println("=== Two-Substrate Two-Product Reversible Enzyme Reaction Derivation ===")
println("Reaction: S + T + E <=> EST <=> E + P + Q")
println()
println("Steady-state enzyme-substrate complex concentration:")
println("EST = (k1f*S*T + k2r*P*Q)*Etot / (k1r + k2f + k1f*S*T + k2r*P*Q)")
println()
println("Net reaction rate:")
println("v = (Vmf*S*T/(Ks*Kt) - Vmr*P*Q/(Kp*Kq)) / (1 + S*T/(Ks*Kt) + P*Q/(Kp*Kq))")
println()
println("Gibbs free energy change:")
println("ΔGr = ΔG° + RT*ln(Q)")
println("where Q = [P][Q]/([S][T])")
println()
println("Flux ratio:")
println("J+/J- = e^(-ΔGr/RT) = Keq/Q")

# Create visualization function
function plot_flux_ratio_vs_gibbs()
    # Create Gibbs free energy change range
    ΔGr_range = -40:0.1:0  # kJ/mol
    
    # Calculate flux ratios for different S/T ratios
    S_T_ratios = [0.2, 1.0, 5.0]
    colors = [:blue, :orange, :gray]
    
    p = plot(xlabel="ΔGr (kJ/mol)", ylabel="J+/J-", 
             yscale=:log10, title="Two-Substrate Enzyme Reaction Flux Ratio",
             legend=:topleft)
    
    for (i, ratio) in enumerate(S_T_ratios)
        # For given S/T ratio, calculate flux ratio
        flux_ratios = exp.(-ΔGr_range ./ (R_val*T_val/1000)) .* ratio
        plot!(p, ΔGr_range, flux_ratios, 
              label="S/T = $ratio", color=colors[i], linewidth=2)
    end
    
    return p
end

# Generate plot
p = plot_flux_ratio_vs_gibbs()
display(p)

# Save plot
savefig(p, "flux_ratio_plot.png")

println("\nPlot saved as flux_ratio_plot.png")
