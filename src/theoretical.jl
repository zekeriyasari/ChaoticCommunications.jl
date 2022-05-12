
export bermacsk

"""
    $TYPEDSIGNATURES 
        
Return probability error for `ebno` for a given `cmap`. `numsamples` is the number of samples per symbol and `numusers` is
the number of users.
"""
function bermacsk(ebno::Real, cmap::AbstractDiscreteOscillator, numsamples::Int, numusers::Int)
    x2 = trajectory!(cmap, numsamples) .^ 2
    Ω = var(x2) / (mean(x2)^2)
    1 / 2 * erfc(
        (2 * Ω / numsamples + 2 * (2 * numusers - 1) / numsamples + 2 * (dbtoval(ebno)^(-1)))^(-1 / 2)
    )
end

function bermacsk(ebno::Real, cmap::AbstractContinuousOscillator, tsymbol::Real, tsample::Real, numusers::Int)
    nsamples = floor(Int, tsymbol / tsample)
    x2 = trajectory!(cmap, tsymbol, tsample) .^ 2
    Ω = var(x2) / (mean(x2)^2)
    1 / 2 * erfc(
        (2 * Ω / nsamples + 2 * (2 * numusers - 1) / nsamples + 2 * (dbtoval(ebno)^(-1)))^(-1 / 2)
    )
end

