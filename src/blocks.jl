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

# NOTE: 
# The difference between the classical and chaotic modulation schemes is that 
# in the classical modulation schemes, for the same message symbol the same signal 
# is transmitted. However, in the chaotic modulation schemes, for the same message symbol 
# different signal is transmitted. In other words, the chaotic modulation schmes do not have 
# a fixed alphabet from which transmission signals can be chosen. Instead, for each message 
# symbol transmission, we have construct different chaotic signals by solving the differential 
# equation of the chaotic systems. 

function (modulator::CSKModulator)(bits)
    vcat(
        map(bits) do bit
            trajectory(choosegen(modulator, bit), modulator.β)
        end...
    )
end
