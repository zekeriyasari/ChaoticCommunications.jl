
export bermacsk

"""
    $SIGNATURES 

Return probability error for `ebno` for a given `cmap`. `numsamples` is the number of samples per symbol and `numusers` is
the number of users.
"""
function bermacsk(ebno, cmap, numsamples, numusers)
    x2 = trajectory!(cmap, numsamples - 1).^2
    Ω = var(x2) / (mean(x2)^2)
    1 / 2 * erfc(
        (2 * Ω / numsamples + 2 * (2 * numusers - 1)/numsamples + 2 * (dbtoval(ebno)^(-1)))^(-1 / 2)
        )
end
