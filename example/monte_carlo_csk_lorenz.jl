using ChaoticCommunications
using PyPlot
pygui(true)
close("all")

# Construct communication system 
k = 1
M = 2^k
ns = Int(5e4)
tsymbol = 10.0
tsample = 0.01
comsys = CommunicationSystem(
    SymbolGenerator(ns, M),
    Modulator(CSK(), [Chen(), Chen()], tsymbol, tsample),
    AWGNChannel(1.0, tsymbol, tsample),
    Detector()
)

# Run simulation 
esno = 0:2:10

# Transmisssion
msg = comsys.symbolgen.symbols
refs, tx = comsys.modulator(msg)

# Rescale symbol energy 
ϵ = 10
γ = √2ϵ

let tx = tx
    for γi ∈ [1, γ]
        tx *= γi
        # Corruption and detection 
        symerr_numerical = map(esno) do val
            @info "Running $val dB"
            comsys.channel.esno = val
            rx = comsys.channel(tx)
            msgext = comsys.detector(refs, rx)
            sum(msg .!= msgext) / length(msg)
        end
        nusers = 1
        symerr_theoretical = map(val -> bermacsk(val, comsys.modulator.gens[1], tsymbol, tsample, nusers), esno)

        # Plots 
        figure()
        plot(esno, symerr_numerical, marker="o", label="Numerical")
        plot(esno, symerr_theoretical, marker="*", label="Theoretial")
        legend()
        yscale("log")
        grid(which="both", axis="both")
        xlabel("EsNo [dB]")
        ylabel("Pe")
        savefig(joinpath(@__DIR__, "bersym_gamma_$(round(γi, digits=3)).pdf"))
    end
end
