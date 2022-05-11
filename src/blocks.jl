# This file includes blocks of the communication systems 

export SymbolGenerator, Modulator, AWGNChannel, Detector, CSK, DCSK, iscontinuous, dbtoval, valtodb,
    CommunicationSystem, ber_numerical, ber_theoretical, getrefs

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
struct Modulator{T<:AbstractScheme,S<:AbstractOscillator}
    "Modulation scheme"
    scheme::T
    "Chaotic generator"
    gens::Vector{S}
    "Symbol duration"
    tsymbol::Float64
    "Sampling period"
    tsample::Float64
end

# When called, a modulator returns the references and transmitted signals.
(mdltr::Modulator)(symbols::AbstractVector, refs=getrefs(mdltr, length(symbols))) = refs, modulate(mdltr, symbols, refs)

"""
    $TYPEDSIGNATURES

Returns the reference signals for `numsymbols` number of symbols. 
To reduces execution time, this function runs in parallel threads.
"""
function getrefs(mdltr::Modulator, numsymbols::Int)
    map(mdltr.gens) do gen
        refgen = Vector{Vector{Float64}}(undef, numsymbols)
        if typeof(gen) <: AbstractContinuousOscillator
            Threads.@threads for i in 1:numsymbols
                refgen[i] = trajectory!(gen, mdltr.tsymbol, mdltr.tsample)
            end
        else
            Threads.@threads for i in 1:numsymbols
                refgen[i] = trajectory!(gen, numsamples(mdltr))
            end
        end
        refgen
    end
end

"""
    $TYPEDSIGNATURES

Modulates the `symbols` to transmitted signals.
"""
function modulate(mdltr::Modulator, symbols::AbstractVector, refs::AbstractVector=getrefs(mdltr, length(symbols)))
    map(enumerate(symbols)) do (i, symbol)
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
    Es = mean(energy.(tx))
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
energy(x, ts=1) = sum(abs.(x) .^ 2) * ts

# ---------------------------------- Detector --------------------------------- # 
"""
    $TYPEDEF

Correlation detector 
"""
struct Detector end

function (detector::Detector)(refs, rx, getstats::Bool=false)
    y = dot.(rx, refs[1] - refs[2])
    symbols = map(yi -> yi ≥ 0 ? 1 : 2, y)
    getstats ? (symbols, mean(y), var(y)) : symbols
end

# ---------------------------------- Communication System --------------------------------- # 
"""
    $TYPEDEF

Communication system ready to be simulated for performance computation. 

# Fields 

    $TYPEDFIELDS
"""
struct CommunicationSystem{T1<:SymbolGenerator,T2<:Modulator,T3<:AWGNChannel,T4<:Detector}
    "Symbol generator"
    symbolgen::T1
    "Modulator"
    modulator::T2
    "Additive white Gaussian channel"
    channel::T3
    "MAP detector"
    detector::T4
end

"""
    $TYPEDSIGNATURES

Simulates `comsys` to compute the symbol error versus `esno` numerically. 
"""
function ber_numerical(comsys::CommunicationSystem, esno::AbstractVector)
    # Transmitter
    message = comsys.symbolgen.symbols
    refs, tx = comsys.modulator(message)

    # Performance computation
    symerr = zeros(length(esno))
    for i = 1:length(symerr)
        @info i
        comsys.channel.esno = esno[i]                       # Update channel noise level
        rx = comsys.channel(tx)                             # Received symbols
        extractedmessage = comsys.detector(refs, rx)        # Extracted message
        symerr[i] = sum(message .!= extractedmessage) / length(message)  # Symbol error
    end

    return symerr
end

"""
    $TYPEDSIGNATURES

Simulates `comsys` to compute the symbol error versus `esno` theoretically. The referecne signals `refs` is constructed if empty.
"""
function ber_theoretical(comsys::CommunicationSystem, esno::AbstractVector; refs::AbstractVector=[])
    # Transmission
    if isempty(refs)
        message = comsys.symbolgen.symbols
        refs, _ = comsys.modulator(message)
    end

    # Performance computation
    mas = zeros(length(esno))
    for i in 1:length(mas)
        rxallones = comsys.channel(refs[1])
        _, μ, σ² = comsys.detector(refs, rxallones, true)
        mas[i] = 1 / 2 * erfc(μ / sqrt(2 * σ²))
    end

    return mas
end

