# This file includes chaotic maps 

export Logistic, Cubic, Tent, Henon, Bernoulli, trajectory, normalize 

abstract type AbstractMap end

for cmap in [:Logistic, :Cubic, :Tent, :Henon, :Bernoulli]
    @eval begin 
        """
            $TYPEDEF 
        # Fields 
            $TYPEDFIELDS
        """
        mutable struct $cmap <: AbstractMap 
            "Initial condition"
            x::Vector{Float64} 
            "Initial time"
            t::Float64 
        end
    end
    if cmap == :Henon 
        @eval $cmap() = $cmap(rand(2), 0.)
    else 
        @eval $cmap() = $cmap(rand(1), 0.)
    end
end

(cmap::Logistic)(dx, x, u, t)  = ( dx[1] = 1 - 2 * x[1]^2 )
(cmap::Cubic)(dx, x, u, t)     = ( dx[1] = 4 * x[1]^3 - 3x[1] )
(cmap::Tent)(dx, x, u, t)      = ( dx[1] = 0 ≤ x[1] ≤ 0.6 ? x[1] / 0.6 : (1 - x[1]) / 0.4 )
(cmap::Henon)(dx, x, u, t)     = ( dx[1] = 1 + x[2] - 1.4 * x[1]^2; dx[2] = 0.3x[1] )
(cmap::Bernoulli)(dx, x, u, t) = ( dx[1] = x[1] < 0 ? 1.2x[1] + 1 : 1.2x[1] - 1 )

"""
   $SIGNATURES

Normalizes (zero mean and unity variance)  `x`.  
"""
normalize(x) = (x .- mean(x)) / std(x)

"""
    $SIGNATURES

Returns the dimension of the state space of `cmap`. 
"""
statedim(cmap::AbstractMap) = length(cmap.x)

"""
    $SIGNATURES

Returns a trajectory of `camp` for a time span of `tspan`. `idx` is the indices of trajectory to be returned. 
"""
function trajectory(cmap::AbstractMap, dt, idx=1)
    sol = solve(DiscreteProblem(cmap, cmap.x, (cmap.t, cmap.t + dt)))
    cmap.t = sol.t[end] 
    cmap.x = sol.u[end]
    map(i -> normalize(getindex.(sol.u, i)), idx) 
end

