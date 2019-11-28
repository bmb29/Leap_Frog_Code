using DrWatson
quickactivate(@__DIR__, "Watson")
using DynamicalBilliards, PyPlot, LinearAlgebra

include(srcdir("plot_perturbationgrowth.jl"))
include(srcdir("unitcells.jl"))

@tagsave(datadir("mushrooms, "Λ_N=$N.bson"), (@dict Λ Λσ ws hs description))