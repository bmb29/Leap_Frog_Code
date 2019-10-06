using DelimitedFiles
using Distributed
using DifferentialEquations
using ProgressMeter
using TaylorIntegration
using Base.Threads

include("PSS_Taylor.jl")
include("Yfind.jl")
include("is_it.jl")
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

Energy=.3
t_end=100.
n_iter_P=101;
n_iter_Q=50;

ArrP=range(1e-8,stop=.99,length=n_iter_P)
ArrQ=range(1e-8,stop=1,length=n_iter_Q)
Q_0=zeros(0)
P_0=zeros(0)

Q_1=zeros(0)
P_1=zeros(0)

Q_2=zeros(0)
P_2=zeros(0)

Q_3=zeros(0)
P_3=zeros(0)

Threads.@threads for j=1:n_iter_P
    for k=1:n_iter_Q
    Q0,P0,Q1,P1,Q2,P2,Q3,P3=PSS_function_taylor(ArrQ[k],ArrP[j], Energy, t_end)
        if Q0!=0 && P0!=0
            push!(Q_0,Q0)
            push!(Q_1,Q1)
            push!(Q_2,Q2)
            push!(Q_3,Q3)

            push!(P_0,P0)
            push!(P_1,P1)
            push!(P_2,P2)
            push!(P_3,P3)
        end
    end
end
figure()
# axis([ -2.5, 2.5,-1.1, 1.1])
Q0=vcat(Q_0,-Q_0,Q_0,-Q_0)
P0=vcat(P_0,P_0,-P_0,-P_0)
T0=hcat(Q0,P0)
plot(P0,Q0 ,color="r",".",markersize=2, markeredgewidth=.1)

Q1=vcat(Q_1,-Q_1,Q_1,-Q_1)
P1=vcat(P_1,P_1,-P_1,-P_1)
T1=hcat(Q1,P1)
plot(P1,Q1 ,color="b",".",markersize=2, markeredgewidth=.1)

Q2=vcat(Q_2,-Q_2,Q_2,-Q_2)
P2=vcat(P_2,P_2,-P_2,-P_2)
T2=hcat(Q2,P2)
plot(P2,Q2 ,color="g",".",markersize=2, markeredgewidth=.1)

Q3=vcat(Q_3,-Q_3,Q_3,-Q_3)
P3=vcat(P_3,P_3,-P_3,-P_3)
T3=hcat(Q3,P3)
plot(P3,Q3 ,color="m",".",markersize=2, markeredgewidth=.1)


outfile_0 = "Taylor_outfile0_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_1 = "Taylor_outfile1_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_2 = "Taylor_outfile2_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_3 = "Taylor_outfile3_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"


a = open(outfile_0, "w")
b = open(outfile_1, "w")
c = open(outfile_2, "w")
d = open(outfile_3, "w")

writedlm(outfile_0, T0)
writedlm(outfile_1, T1)
writedlm(outfile_2, T2)
writedlm(outfile_3, T3)



close(a)
close(b)
close(c)
close(d)
