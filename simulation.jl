using DifferentialEquations
using Plots
using Statistics

# Define the ODE system manually
function two_substrate_system!(du, u, p, t)
    S, T, E, EST, P, Q = u
    k1f, k1r, k2f, k2r = p
    
    # Reaction rates
    r1f = k1f * S * T * E
    r1r = k1r * EST
    r2f = k2f * EST
    r2r = k2r * P * Q * E
    
    # Differential equations
    du[1] = -r1f + r1r  # dS/dt
    du[2] = -r1f + r1r  # dT/dt
    du[3] = -r1f + r1r + r2f - r2r  # dE/dt
    du[4] = r1f - r1r - r2f + r2r  # dEST/dt
    du[5] = r2f - r2r  # dP/dt
    du[6] = r2f - r2r  # dQ/dt
end

# Set parameter values
p = [1.0, 0.1, 2.0, 0.05]  # [k1f, k1r, k2f, k2r]

# Initial conditions [S, T, E, EST, P, Q]
u0 = [10.0, 5.0, 1.0, 0.0, 0.0, 0.0]

# Time range
tspan = (0.0, 10.0)

# Create ODE problem
prob = ODEProblem(two_substrate_system!, u0, tspan, p)

# Solve
sol = solve(prob, Tsit5(), saveat=0.1)

println("=== Numerical Simulation Results ===")
println("Simulation completed successfully!")
println("Final concentrations:")
println("  S = $(sol[1,end])")
println("  T = $(sol[2,end])")
println("  E = $(sol[3,end])")
println("  EST = $(sol[4,end])")
println("  P = $(sol[5,end])")
println("  Q = $(sol[6,end])")

# Create plot
p = plot(sol.t, sol[1,:], label="S", linewidth=2)
plot!(p, sol.t, sol[2,:], label="T", linewidth=2)
plot!(p, sol.t, sol[3,:], label="E", linewidth=2)
plot!(p, sol.t, sol[4,:], label="EST", linewidth=2)
plot!(p, sol.t, sol[5,:], label="P", linewidth=2)
plot!(p, sol.t, sol[6,:], label="Q", linewidth=2)
xlabel!("Time")
ylabel!("Concentration")
title!("Two-Substrate Enzyme Reaction Simulation")

savefig(p, "simulation_results.png")
println("\nPlot saved as simulation_results.png") 