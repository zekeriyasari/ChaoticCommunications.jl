using ChaoticCommunications 
using SpecialFunctions
using Plots 

# Settings 
k = 1 
M = 2^k
nusers = 1
cmap = Bernoulli
ns = Int(1e4)
ebno = collect(0 : 1 : 10)
tsymbol = 250. 
tsample = 1.
nsamples = floor(Int, tsymbol / tsample)

# Blocks 
gen = SymbolGenerator(ns, M) 
modulator = Modulator(CSK(), [cmap(), cmap()], tsymbol, tsample)
channel = AWGNChannel(1., tsymbol, tsample) 
detector = Detector()

# Run communication system 
message = gen.symbols
refs, tx = modulator(message) 
symerr = zeros(length(ebno))
for i = 1 : length(symerr) 
    # Monte carlo simulation 
    channel.esno = ebno[i]
    rx = channel(tx)
    extractedmessage = detector(refs, rx)
    symerr[i] = sum(message .!= extractedmessage) / ns 
end 

# Compute theoretical bervals 
bervals = map(val -> bermacsk(val, modulator.gens[1], nsamples, nusers), ebno)

# Plots 
plt = plot(legend=:bottomleft) 
plot!(ebno, symerr, yscale=:log10, lw=0.5, 
    markershape=:auto, color=:black, gridalpha=0.9, 
    minorgrid=true, minorgridalpha=0.5, 
    label="CSK-$(typeof(cmap))-montecarlo")
plot!(ebno, bervals, yscale=:log10, lw=0.5, 
    markershape=:auto, color=:black, gridalpha=0.9, 
    minorgrid=true, minorgridalpha=0.5, 
    label="CSK-$(typeof(cmap))-theoretical")
