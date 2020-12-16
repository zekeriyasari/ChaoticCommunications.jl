# This file includes blocks of the communication systems 

export Generator, CSKModulator

"""
    $TYPEDEF

Bit generator 

# FIELDS 

    $TYPEDFIELDS
"""
struct Generator
    "Generated bits"
    bits::Vector{Bool}
    Generator(nbits::Int) = new(rand(Bool, nbits)) 
end

"""
    $TYPEDEF

Chaos-Shif-Keying Modulator 

# FIELDS

    $TYPEDFIELDS
"""
struct CSKModulator{T} 
    "Chaotic generator"
    gens::Vector{T}
    "Spreading factor"
    β::Int
end 

choosegen(modulator::CSKModulator, idx::Bool) = idx ? modulator.gens[1] : modulator.gens[2]

function (modulator::CSKModulator)(bits)
    β = modulator.β
    vcat(
        map(enumerate(bits)) do (l, bit)
            trajectory(choosegen(modulator, bit), ((l - 1) * β, l * β - 1))
        end...
    )
end
