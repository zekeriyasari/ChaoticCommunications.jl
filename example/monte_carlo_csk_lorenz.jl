using ChaoticCommunications 
using SpecialFunctions
using Plots 

# Settings 
k = 1 
M = 2^k
ns = Int(1e5)
esno = collect(0 : 2 : 8)
tsymbol = 10. 
tsample = 0.005 

# Blocks 
gen = SymbolGenerator(ns, M) 
modulator = Modulator(CSK(), [Lorenz(), Lorenz()], tsymbol, tsample)
channel = AWGNChannel(1., tsymbol, tsample) 
detector = Detector()

# Run communication system 
message = gen.symbols
refs, tx = modulator(message) 
symerr = zeros(length(esno))
mas = zeros(length(esno))
for i = 1 : length(symerr) 
    @info i 
    # Monte carlo simulation 
    channel.esno = esno[i]
    rx = channel(tx)
    extractedmessage = detector(refs, rx)
    symerr[i] = sum(message .!= extractedmessage) / ns 
    # MAS simulation 
    rxallones = channel(refs[1])
    _, μ, σ² = detector(refs, rxallones, true)
    mas[i] = 1 / 2 * erfc(μ / sqrt(2 * σ²))
end 

# Plots 
plt = plot(legend=:bottomleft) 
plot!(esno, symerr, marker=(:circle, 3), yscale=:log10, label="montecarlo")
plot!(esno, mas, marker=(:circle, 3), yscale=:log10, label="mas")
xlabel!("esno") 
ylabel!("Pe") 
title!("Single User CSK") 
