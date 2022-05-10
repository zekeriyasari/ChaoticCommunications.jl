
using ChaoticCommunications
using PyPlot
using Statistics
pygui(true)
close("all")

n = 100
Tc = 0.01
Tb = 20.0
ϵ = 10.0
cmaps = [Chen(), Chen()]
x1, x2 = [map(i -> trajectory!(cmap, Tb, Tc, normalized=false), 1:n) for cmap in cmaps]
err = ϵ * (x1 - x2)

Eb1 = map(item -> sum(item .^ 2), x1) |> mean
Eb2 = map(item -> sum(item .^ 2), x2) |> mean
Eberr = map(item -> sum(item .^ 2), err) |> mean

@show Eb1, Eb2, Eberr

figure()
plot(x1[1])
