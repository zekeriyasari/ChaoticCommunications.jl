# This file includes chaotic maps 

export Logistic, Cubic, Tent, Henon, Bernoulli,
    Lorenz, Chua, Rossler, Chen,
    trajectory!, normalize

abstract type AbstractOscillator end
abstract type AbstractDiscreteOscillator <: AbstractOscillator end
abstract type AbstractContinuousOscillator <: AbstractOscillator end

# ----------------------------- Discrete Time Maps --------------------------------- # 

mutable struct Logistic <: AbstractDiscreteOscillator
    a::Float64
    x::Vector{Float64}
    t::Float64
end
Logistic() = Logistic(2, rand(1), 0.0)
(cmap::Logistic)(dx, x, u, t) = (dx[1] = 1 - cmap.a * x[1]^2)

mutable struct Cubic <: AbstractDiscreteOscillator
    a::Float64
    b::Float64
    x::Vector{Float64}
    t::Float64
end
Cubic() = Cubic(4, 3, rand(1), 0.0)
(cmap::Cubic)(dx, x, u, t) = (dx[1] = cmap.a * x[1]^3 - cmap.b * x[1])

mutable struct Tent <: AbstractDiscreteOscillator
    a::Float64
    b::Float64
    c::Float64
    x::Vector{Float64}
    t::Float64
end
Tent() = Tent(0.6, 0.6, 0.4, rand(1), 0.0)
(cmap::Tent)(dx, x, u, t) = (dx[1] = 0 ≤ x[1] ≤ cmap.a ? x[1] / cmap.b : (1 - x[1]) / cmap.c)

mutable struct Henon <: AbstractDiscreteOscillator
    a::Float64
    b::Float64
    x::Vector{Float64}
    t::Float64
end
Henon() = Henon(1.4, 0.3, rand(2), 0.0)

function (cmap::Henon)(dx, x, u, t)
    dx[1] = 1 + x[2] - cmap.a * x[1]^2
    dx[2] = cmap.b * x[1]
end

mutable struct Bernoulli <: AbstractDiscreteOscillator
    a::Float64
    x::Vector{Float64}
    t::Float64
end
Bernoulli() = Bernoulli(1.2, rand(1), 0.0)
(cmap::Bernoulli)(dx, x, u, t) = (dx[1] = x[1] < 0 ? cmap.a * x[1] + 1 : cmap.a * x[1] - 1)

# ----------------------------- Continuous Time Maps --------------------------------- # 

mutable struct Lorenz <: AbstractContinuousOscillator
    σ::Float64
    β::Float64
    ρ::Float64
    γ::Float64
    x::Vector{Float64}
    t::Float64
end
Lorenz() = Lorenz(10, 8 / 3, 28, 1, rand(3), 0)

function (cmap::Lorenz)(dx, x, u, t)
    dx[1] = cmap.σ * (x[2] - x[1])
    dx[2] = x[1] * (cmap.ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - cmap.β * x[3]
end

Base.@kwdef mutable struct Chen <: AbstractContinuousOscillator
    a::Float64 = 35.0
    c::Float64 = 28.0
    β::Float64 = 8 / 3
    γ::Float64 = 1.0
    x::Vector{Float64} = rand(3)
    t::Float64 = 0.0
end

function (cmap::Chen)(dx, x, u, t)
    dx[1] = cmap.a * (x[2] - x[1])
    dx[2] = (cmap.c - cmap.a - x[3]) * x[1] + cmap.c * x[2]
    dx[3] = x[1] * x[2] - cmap.β * x[3]
    dx .*= cmap.γ
end


struct Diode
    a::Float64
    b::Float64
    bp::Float64
end
Diode() = Diode(-1.143, -0.714, 1.0)
function (diode::Diode)(x)
    if x < -cmap.bp
        cmap.b * x + (cmap.b - cmap.a) * cmap.bp
    elseif -cmap.bp ≤ x ≤ cmap.bp
        cmap.a * x
    elseif x > cmap.bp
        cmap.b * x + (cmap.a - cmap.b) * bp1
    end
end

mutable struct Chua <: AbstractContinuousOscillator
    α::Float64
    β::Float64
    h::Diode
    x::Vector{Float64}
    t::Float64
end
Chua() = Chua(15, 28, Diode(), rand(3), 0.0)

function (cmap::Chua)(dx, x, u, t)
    dx[1] = cmap.α * (x[2] - x[1] - cmap.h(x[1]))
    dx[2] = x[1] - x[2] + x[3]
    dx[3] = -β * x[2]
end


mutable struct Rossler <: AbstractContinuousOscillator
    a::Float64
    b::Float64
    c::Float64
    x::Vector{Float64}
    t::Float64
end
Rossler() = Rossler(0.38, 0.3, 4.82, rand(3), 0.0)

function (cmap::Rossler)(dx, x, u, t)
    dx[1] = -x[2] - x[3]
    dx[2] = x[1] + cmap.a * x[2]
    dx[3] = cmap.b + x[3] * (x[1] - cmap.c)
end

# --------------------------------------- Methods  -------------------------------------- # 

"""
   $SIGNATURES

Normalizes (zero mean and unity variance)  `x`.  
"""
normalize(x) = (x .- mean(x)) / std(x)

"""
    $SIGNATURES

Returns the dimension of the state space of `cmap`. 
"""
statedim(cmap::AbstractOscillator) = length(cmap.x)

"""
    $SIGNATURES

Returns a trajectory! of `camp` for a time span of `trange`. `idx` is the indices of trajectory! to be returned. 
"""
function trajectory! end

function trajectory!(cmap::AbstractDiscreteOscillator, trange::Real, idx::Union{<:AbstractVector,<:Int}=1; normalized::Bool=true)
    probf = (dx, x, u, t) -> cmap(dx, x, u, t)
    sol = solve(DiscreteProblem(probf, cmap.x, (cmap.t, cmap.t + trange)))
    cmap.t = sol.t[end]
    cmap.x = sol.u[end]
    normalized ? map(i -> normalize(getindex.(sol.u, i)), idx) : map(i -> getindex.(sol.u, i), idx)
end

function trajectory!(cmap::AbstractContinuousOscillator, trange::Real, tsample::Real, idx::Union{<:AbstractVector,<:Int}=1; normalized::Bool=true)
    sol = solve(ODEProblem(cmap, cmap.x, (cmap.t, cmap.t + trange)), saveat=tsample)
    ns = floor(Int, trange / tsample) + 1
    if length(sol) > ns
        sol = sol[1:ns]
    end
    cmap.t = sol.t[end]
    cmap.x = sol.u[end]
    normalized ? map(i -> normalize(getindex.(sol.u, i)), idx) : map(i -> getindex.(sol.u, i), idx)
end


