# This file is used for case studies.

using ChaoticCommunications
using Plots 

# Simulation settinngs 
nbits = 10
β = 100 

# Construct blocks 
generator = Generator(nbits) 
modulator = CSKModulator([Logistic(), Logistic()], β)

# Runs simulation 
# tx = trues(nbits) |> modulator 
tx = generator.bits |> modulator 

# Plots 
plot(tx, marker=(:circle, 1))
