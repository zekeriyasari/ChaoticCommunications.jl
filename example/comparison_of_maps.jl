using ChaoticCommunications
using SpecialFunctions 
using Plots 
using Statistics

# Settings 
for nusers in [1, 3, 5]
    nsamples = 100
    ebno = 0 : 1 : 20 
    theme(:default)
    plt = plot(legend=:bottomleft) 
    for cmap in [Logistic(), Cubic(), Tent(), Henon(), Bernoulli()]
        # Plots 
        bervals = map(val -> bermacsk(val, cmap, nsamples, nusers), ebno)
        plot!(ebno, bervals, yscale=:log10, lw=0.5, 
                        markershape=:auto, color=:black, gridalpha=0.9, 
                        minorgrid=true, minorgridalpha=0.5, 
                        label="CSK-$(typeof(cmap))-theoretical")
    end 
    display(plt) 
    savefig(joinpath(@__DIR__, "comparison_of_maps_N_$nusers.pdf"))
end 
