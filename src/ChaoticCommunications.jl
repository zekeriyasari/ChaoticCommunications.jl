#= 
    A package for chaotic communications 
=# 
module ChaoticCommunications

using DifferentialEquations 
using Statistics 
using LinearAlgebra 
using SpecialFunctions
using DocStringExtensions

include("maps.jl") 
include("theoretical.jl") 
include("blocks.jl") 

end # module
