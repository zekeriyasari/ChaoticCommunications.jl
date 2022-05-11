using ChaoticCommunications
using SpecialFunctions
using PyPlot
pygui(true)
close("all")

# Construct communication system 
k = 1
M = 2^k
ns = Int(1e4)
tsymbol = 10.0
tsample = 0.01
comsys = CommunicationSystem(
    SymbolGenerator(ns, M),
    Modulator(CSK(), [Lorenz(), Lorenz()], tsymbol, tsample),
    AWGNChannel(1.0, tsymbol, tsample),
    Detector()
)

# Run simulation 
esno = 0:2:8
symerr = ber_numerical(comsys, esno)

# Plots 
plot(esno, symerr, marker="o")
yscale("log")
grid(which="both", axis="both")
xlabel("EsNo [dB]")
ylabel("Pe")
