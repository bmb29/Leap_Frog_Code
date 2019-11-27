 using Distributed
@everywhere using DifferentialEquations
@everywhere using ProgressMeter

include("PSS_function_parallel.jl")
include("Yfind.jl")
# include("is_it.jl")

@everywhere using PyCall
pygui(:qt)
@everywhere using PyPlot
pygui(true)
#PSS_function(Q,P, Energy, t_end, Max_hit)
@everywhere max_hit=10
@everywhere Energy=.25
@everywhere H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
@everywhere Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
@everywhere ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
@everywhere ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
@everywhere function Eq_of_M(du,u,p,t)
     du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
     du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
     du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
     du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
     du[5]=0
 return
end
@everywhere condition(u,t,integrator)= u[5]>max_hit
@everywhere function condition2(u,t,integrator) # Event when event_f(u,t) == 0
 u[1]
end
@everywhere affect!(integrator) = terminate!(integrator)
@everywhere function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
 end
@everywhere cb2 = ContinuousCallback(condition2,nothing,affect2!,  save_positions = (false,true))
@everywhere cb1 = DiscreteCallback(condition,affect!)

@everywhere cb=CallbackSet(cb2,cb1)
@everywhere H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))


n_iter_P=200
n_iter_Q=200

N=n_iter_P*n_iter_Q

T=@showprogress pmap(PSS_function_parallel,1:N)
COLOR=["b", "g", "c", "m", "y", "k", "r"]
#
# outfileQ = "PSSoutfileQ_"*string(Energy)*"nP"*string(n_iter_P)*"nQ"*string(n_iter_Q)*".dat"
# outfileP = "PSSoutfileP_"*string(Energy)*"nP"*string(n_iter_P)*"nQ"*string(n_iter_Q)*".dat"
# outfilePAIR = "PSSoutfileDuple_"*string(Energy)*"nP"*string(n_iter_P)*"nQ"*string(n_iter_Q)*".dat"
#
# a = open(outfileQ, "w")
# b = open(outfileP, "w")
# c = open(outfilePAIR, "w")
# Q_tabs=zeros(0)
# P_tabs=zeros(0)
# QQ=zeros(max_hit+2,N)
# PP=zeros(max_hit+2,N)
figure()
axis([ -2.5, 2.5,-1.5, 1.5])
# axis([ -.02,.02,-.2, .2])

@showprogress 1 "Computing..." for k=1:N
    if T[k]!=nothing
        current_color=COLOR[mod(k,length(COLOR))+1]
        plot(T[k][2],T[k][1], color=current_color,".",markersize=1, markeredgewidth=.1)
        plot(T[k][2],-T[k][1], color=current_color,".",markersize=1, markeredgewidth=.1)
        plot(-T[k][2],T[k][1], color=current_color,".",markersize=1, markeredgewidth=.1)
        plot(-T[k][2],-T[k][1], color=current_color,".",markersize=1, markeredgewidth=.1)
        # for j=1:length(T[k][1])
        #     # push!(Q_tabs,T[k][1][j])
        #     # push!(P_tabs,T[k][2][j])
        #     QQ[j,k]=T[k][1][j]
        #     PP[j,k]=T[k][2][j]
        # end
    end
end
# Q=vcat(Q_tabs,-Q_tabs,Q_tabs,-Q_tabs)
# P=vcat(P_tabs,P_tabs,-P_tabs,-P_tabs)
# # PAIR=hcat(Q,P)
# plot(P,Q ,color="k",".",markersize=1, markeredgewidth=.1)
# writedlm(outfileQ, QQ)
# writedlm(outfileP, PP)
#
# # writedlm(outfilePAIR, PAIR)
#
# close(a)
# close(b)
