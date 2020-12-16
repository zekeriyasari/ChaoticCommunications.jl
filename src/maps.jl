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
        struct $cmap <: AbstractMap 
            "Initial condition"
            x0::Vector{Float64} 
        end
    end
    if cmap == :Henon 
        @eval $cmap() = $cmap(rand(2))
    else 
        @eval $cmap() = $cmap(rand(1))
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

Returns a trajectory of `camp` for a time span of `tspan`.

!!! note 
    If state space of `cmap` is multidimensional, the first dimension is returned. 
"""
function trajectory(cmap::AbstractMap, tspan)
    sol = solve(DiscreteProblem(cmap, cmap.x0, tspan))
    normalize(getindex.(sol.u, 1)) 
end

