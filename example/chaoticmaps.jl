# This file includes an example file to simulate chaotic maps 

using ChaoticCommunications
using Plots

# Generate data  
cmaps = [Logistic(), Tent(), Logistic(), Henon(), Bernoulli()]
plt = plot(layout=(length(cmaps), 1))
for (i, cmap) in enumerate(cmaps)
    tx = trajectory!(cmap, 100.)
    plot!(tx, marker=(:circle, 1), subplot=i)
end 
display(plt)

