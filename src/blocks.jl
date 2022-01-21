# This file includes blocks of the communication systems 

export SymbolGenerator, Modulator, AWGNChannel, Detector, CSK, DCSK, iscontinuous, dbtoval, valtodb

# ------------------------------ Symbol Generator --------------------------------- # 

"""
    $TYPEDEF

Bit generator 

# FIELDS 

    $TYPEDFIELDS
"""
struct SymbolGenerator
    "Generated symbols"
    symbols::Vector{Int}
    SymbolGenerator(ns::Int, M::Int) = new(rand(1:M, ns)) 
end

# ------------------------------ CSK Modulator --------------------------------- # 

abstract type AbstractScheme end

"""
    $TYPEDEF

Chaos Shift Keying scheme 
"""
struct CSK <: AbstractScheme end 

"""
    $TYPEDEF

Differential Chaos Shift Keying
""" 
struct DCSK <: AbstractScheme end 

"""
    $TYPEDEF

Chaos-Shif-Keying Modulator 

# FIELDS

    $TYPEDFIELDS
"""
struct Modulator{T<:AbstractScheme, S<:AbstractOscillator} 
    "Modulation scheme"
    scheme::T 
    "Chaotic generator"
    gens::Vector{S}
    "Symbol duration"
    tsymbol::Float64 
    "Sampling period"
    tsample::Float64 
end 

(mdltr::Modulator)(symbols) = iscontinuous(mdltr) ? contmodulate!(mdltr, symbols) : discmodulate!(mdltr, symbols) 

function discmodulate!(modulator::Modulator, symbols)
    K = numsamples(modulator) 
    refs = map(gen -> trajectory!(gen, K - 1), modulator.gens)
    refs, map(enumerate(symbols)) do (i, symbol) 
        refs[symbols][i]
    end 
end

function contmodulate!(modulator::Modulator, symbols)
    numsymbols = length(symbols)
    tf = modulator.tsymbol 
    ts = modulator.tsample
    refs = map(gen -> [trajectory!(gen, tf - ts, ts) for i in 1 : numsymbols], modulator.gens)
    refs, map(enumerate(symbols)) do (i, symbol) 
        refs[symbol][i] 
    end 
end

"""
    $SIGNATURES 

Returns true if `modulator` operates in continuous time 
"""
iscontinuous(modulator::Modulator) = eltype(modulator.gens) <: AbstractContinuousOscillator

"""
    $SIGNATURES

Returns the number of samples per symbol of `modulator` 
"""
numsamples(modulator::Modulator) = floor(Int, modulator.tsymbol / modulator.tsample)


# ------------------------------ Channel  --------------------------------- # 

"""
    $SIGNATURES 

Additive white Gaussian noise channel 

# Fields 
    $TYPEDFIELDS
"""
mutable struct AWGNChannel 
    "Energy per symbol to noise power spectral density ratio (in dB)"
    esno::Float64 
    "Symbol duration"
    tsymbol::Float64 
    "Sampling period"
    tsample::Float64 
end

function (channel::AWGNChannel)(tx) 
    numsymbols = length(tx) 
    numsamples = length(tx[1]) 
    # fs = 1 / channel.tsample
    Es = mean(energy.(tx))
    # σ = sqrt(Es * fs / 2 / dbtoval(channel.esno))
    σ = sqrt(Es / 2 / dbtoval(channel.esno))
    n = collect(eachrow(σ * randn(numsymbols, numsamples)))
    tx + n
end 

"""
    $SIGNATURES 

Converts `dB` to numerical value.
"""
dbtoval(db) = 10^(db / 10) 

"""
    $SIGNATURES 

Converts `dB` to numerical value.
"""
valtodb(db) = 10 * log10(db) 

"""
    $SIGNATURES 

Returns the energy of `x` sampled with a sampling frequency `ts`.
"""
energy(x, ts=1) = sum(abs.(x).^2) * ts

# ---------------------------------- Detector --------------------------------- # 
"""
    $TYPEDEF

Correlation detector 
"""
struct Detector end

function (detector::Detector)(refs, rx, getstats::Bool=false)
    y = map(rxi -> dot(rxi, (refs[1] - refs[2])), rx)
    symbols = map(yi -> yi ≥ 0 ? 1 : 2,  y)
    getstats ? (symbols, mean(y), var(y)) : symbols
    # map(enumerate(rx)) do (i, rxi) 
    #     argmax(map(ref -> rxi ⋅ ref[i], refs)) 
    # end
end
