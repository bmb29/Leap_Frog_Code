
using DelimitedFiles
using DifferentialEquations
using ProgressMeter
using PyCall
include("PSS_function_BIG.jl")
pygui(:qt)
using PyPlot
pygui(true)


Energy=BigFloat("0.125789")
H=(2*Energy)^2
max_hit=BigInt(100000)
t_end=BigFloat("1000")


Q=BigFloat(".0000000001")
P=BigFloat(".0000000001")

Y,Q_PSS, P_PSS=PSS_function_BIG(Q,P,Energy,t_end, max_hit)
plot(P_PSS,Q_PSS, color="k",".",markersize=1, markeredgewidth=.1)
