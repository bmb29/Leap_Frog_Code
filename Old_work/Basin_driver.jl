using DelimitedFiles
using DifferentialEquations
using ProgressMeter
include("PSS_vern9.jl")
 using PyCall
pygui(:qt)
using PyPlot
pygui(true)
max_hit=2
Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
return
end
condition(u,t,integrator)= u[5]>max_hit
function condition2(u,t,integrator) # Event when event_f(u,t) == 0
u[1]
end
affect!(integrator) = terminate!(integrator)
function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
end
cb2 = ContinuousCallback(condition2,affect2!,nothing)
cb1 = DiscreteCallback(condition,affect!)

cb=CallbackSet(cb2,cb1)
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))


Energy=.25
t_end=100.
n_iter_P=500
n_iter_Q=501

# ArrP=range(1.1,stop=2,length=n_iter_P)
# ArrQ=range(1e-10,stop=1,length=n_iter_Q)
# ArrP=range(0,stop=.1,length=n_iter_P)
# ArrQ=range(1e-10,stop=.1,length=n_iter_Q)
ArrP=zeros(n_iter_P)
ArrQ=zeros(n_iter_Q)
Q_bound=zeros(0)
P_bound=zeros(0)

Q_0=zeros(0)
P_0=zeros(0)

Q_1=zeros(0)
P_1=zeros(0)

Q_2=zeros(0)
P_2=zeros(0)

Q_3=zeros(0)
P_3=zeros(0)

Q_4=zeros(0)
P_4=zeros(0)

Q_5=zeros(0)
P_5=zeros(0)

@showprogress 1 "Computing..." for j=1:n_iter_P
    for k=1:n_iter_Q
        # ArrQ[k]=.1*rand()
        # ArrP[j]=.01*rand()
        K=PSS_vern9(ArrQ[k],ArrP[j], Energy, t_end)
        if K==0
            push!(Q_0,ArrQ[k])
            push!(P_0,ArrP[j])
        elseif K==1
            push!(Q_1,ArrQ[k])
            push!(P_1,ArrP[j])
        elseif K==2
            push!(Q_2,ArrQ[k])
            push!(P_2,ArrP[j])
        elseif K==3
            push!(Q_3,ArrQ[k])
            push!(P_3,ArrP[j])
        elseif K==4
            push!(Q_4,ArrQ[k])
            push!(P_4,ArrP[j])

        elseif K==-1
            push!(Q_bound,ArrQ[k])
            push!(P_bound,ArrP[j])
        end
    end
end

figure()
# axis([ -2.5, 2.5,-1.1, 1.1])
Qbound=vcat(Q_bound,-Q_bound,Q_bound,-Q_bound)
Pbound=vcat(P_bound,P_bound,-P_bound,-P_bound)
Tbound=hcat(Qbound,Pbound)
plot(Pbound,Qbound ,color="k",".",markersize=1, markeredgewidth=.1)

Q0=vcat(Q_0,-Q_0,Q_0,-Q_0)
P0=vcat(P_0,P_0,-P_0,-P_0)
T0=hcat(Q0,P0)
plot(P0,Q0 ,color="r",".",markersize=1, markeredgewidth=.1)

Q1=vcat(Q_1,-Q_1,Q_1,-Q_1)
P1=vcat(P_1,P_1,-P_1,-P_1)
T1=hcat(Q1,P1)
plot(P1,Q1 ,color="r",".",markersize=1, markeredgewidth=.1)

Q2=vcat(Q_2,-Q_2,Q_2,-Q_2)
P2=vcat(P_2,P_2,-P_2,-P_2)
T2=hcat(Q2,P2)
plot(P2,Q2 ,color="g",".",markersize=1, markeredgewidth=.1)

Q3=vcat(Q_3,-Q_3,Q_3,-Q_3)
P3=vcat(P_3,P_3,-P_3,-P_3)
T3=hcat(Q3,P3)
plot(P3,Q3 ,color="m",".",markersize=1, markeredgewidth=.1)

Q4=vcat(Q_4,-Q_4,Q_4,-Q_4)
P4=vcat(P_4,P_4,-P_4,-P_4)
T4=hcat(Q4,P4)
plot(P4,Q4 ,color="c",".",markersize=1, markeredgewidth=.1)

Q5=vcat(Q_5,-Q_5,Q_5,-Q_5)
P5=vcat(P_5,P_5,-P_5,-P_5)
T5=hcat(Q5,P5)
plot(P5,Q5 ,color="y",".",markersize=1, markeredgewidth=.1)

outfile_0 = "Vern9_outfile0_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_1 = "Vern9_outfile1_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_2 = "Vern9_outfile2_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_3 = "Vern9_outfile3_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_4 = "Vern9_outfile4_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_5 = "Vern9_outfile5_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
outfile_bound = "Vern9_outfile1_bound"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"


a = open(outfile_0, "w")
b = open(outfile_1, "w")
c = open(outfile_2, "w")
d = open(outfile_3, "w")
e = open(outfile_4, "w")
f = open(outfile_5, "w")
g = open(outfile_bound, "w")

writedlm(outfile_0, T0)
writedlm(outfile_1, T1)
writedlm(outfile_2, T2)
writedlm(outfile_3, T3)
writedlm(outfile_4, T4)
writedlm(outfile_5, T5)

writedlm(outfile_bound, Tbound)


close(a)
close(b)
close(c)
close(d)
close(e)
close(f)
close(g)
