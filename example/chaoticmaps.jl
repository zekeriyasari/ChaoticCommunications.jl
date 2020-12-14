# Simple file to simulate discrete time chaotic systems.

using DifferentialEquations 
using Plots 
using Statistics

# -------------------------------- Chaotic Maps ------------------------------------ # 

function logistic(dx, x, u, t) 
    dx[1] = 1 - 2 * x[1]^2
end

function cubic(dx, x, u, t) 
    dx[1] = 4 * x[1]^3 - 3x[1]
end

function tent(dx, x, u, t) 
    dx[1] = 0 ≤ x[1] ≤ 0.6 ? x[1] / 0.6 : (1 - x[1]) / 0.4
end

function henon(dx, x, u, t)
    dx[1] = 1 + x[2] - 1.4 * x[1]^2
    dx[2] = 0.3x[1]
end

function bernoulli(dx, x, u, t) 
    dx[1] = x[1] < 0 ? 1.2x[1] + 1 : 1.2x[1] - 1
end

# -------------------------------- Normalization -------------------------------------- # 
"Returns normalized (with zero mean and unity variance) `x`."
normalize(x) = (x .- mean(x)) / std(x)

"Displays some of the statistical properties of `x`"
function statprops(x) 
    @show mean(x) 
    @show mean(x.^2) 
    @show mean(x.^3) 
    @show var(x) 
    @show var(x.^2) 
    nothing 
end

# -------------------------------- Data generation ------------------------------------ # 

f = tent
tspan = (0., 1e4) # Note: tspan should be (0., 1e5) to get the table in the benchmarks. 
if f == henon  
    x0 = rand(2) 
else 
    x0 = rand(1)
end
prob = DiscreteProblem(f, x0, tspan) 
sol = solve(prob) 

# ----------------------------------- Plots ------------------------------------------- #
if length(sol.u[1]) == 1 
    plt = plot(sol.t, getindex.(sol.u, 1))
else length(sol.u[1]) == 2
    ms = 1
    l = @layout [a c; b]
    plt = plot(layout=l)
    plot!(sol.t, getindex.(sol.u, 1), subplot=1, ms=ms) 
    plot!(sol.t, getindex.(sol.u, 2), subplot=2, ms=ms) 
    scatter!(getindex.(sol.u, 1), getindex.(sol.u, 2), ms=ms, subplot=3) 
end 
display(plt)

statprops(normalize(getindex.(sol.u, 1)))
