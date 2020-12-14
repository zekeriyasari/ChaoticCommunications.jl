# Simple file to simulation discrete time chaotic systems.

using DifferentialEquations 
using Plots 

function f(dx, x, u, t) 
    dx[1] = -x[1] 
end

x0 = ones(1) 
tspan = (0., 10.) 
prob = ODEProblem(f, x0, tspan) 
sol = solve(prob) 
