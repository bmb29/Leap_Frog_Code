using Distributed
@everywhere using DelimitedFiles
@everywhere using DifferentialEquations
@everywhere using ProgressMeter
include("where_did_it_come_from_parallel.jl")
@everywhere using PyCall
pygui(:qt)
@everywhere using PyPlot
pygui(true)

@everywhere max_hit=3
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
@everywhere condition(u,t,integrator)= u[5]==max_hit
@everywhere function condition2(u,t,integrator) # Event when event_f(u,t) == 0
u[1]
end
@everywhere condition3b(u,t,integrator)=  u[5]==2
@everywhere condition3f(u,t,integrator)=  u[5]==2

@everywhere affect!(integrator) = terminate!(integrator)
@everywhere function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
end
# @everywhere cb2 = ContinuousCallback(condition2,affect2!, save_positions = (false,true) )
@everywhere cb2 = ContinuousCallback(condition2, affect2!, rootfind=true, save_positions = (true,false))
@everywhere cb2f = ContinuousCallback(condition2, affect2!, rootfind=true,  save_positions = (true,false))
@everywhere cb2b = ContinuousCallback(condition2, affect2!,rootfind=true, save_positions = (true,false))

@everywhere cb1 = DiscreteCallback(condition,affect!)
@everywhere cb3b = DiscreteCallback(condition3b,affect!)
@everywhere cb3f = DiscreteCallback(condition3f,affect!)

@everywhere cb=CallbackSet(cb2,cb1)

@everywhere cb_backward=CallbackSet(cb2b,cb3b)
@everywhere cb_forward=CallbackSet(cb2f,cb3f)

@everywhere X_check=zeros(0)
@everywhere Y_check=zeros(0)
@everywhere Q_check=zeros(0)
@everywhere P_check=zeros(0)
@everywhere hit_check1=zeros(0)
@everywhere hit_check2=zeros(0)

@everywhere Q_init=zeros(0)
@everywhere P_init=zeros(0)
@everywhere Y_initS=zeros(0)
@everywhere Y_initB=zeros(0)
@everywhere test=zeros(0)


@everywhere Y_from_nowhere=zeros(0)
@everywhere Q_from_nowhere=zeros(0)
@everywhere P_from_nowhere=zeros(0)

@everywhere X_1before=zeros(0)
@everywhere Y_1before=zeros(0)
@everywhere Q_1before=zeros(0)
@everywhere P_1before=zeros(0)
# ArrQ[k]=.1*rand()
# ArrP[j]=.01*rand()
@everywhere Q_end=1.5
@everywhere Q_start=-Q_end
@everywhere P_end=3
@everywhere P_start=-P_end
@everywhere n_iter_P=8001
@everywhere n_iter_Q=8001
@everywhere N=n_iter_P*n_iter_Q

@everywhere delta_P=(P_end-P_start)/(n_iter_P)
@everywhere delta_Q=(Q_end-Q_start)/(n_iter_Q)

@everywhere ArrP=P_start: delta_P:P_end
@everywhere ArrQ=Q_start: delta_Q: Q_end

T=@showprogress pmap(where_did_it_come_from_parallel,1:N)


for i=1:N
    Q=ArrQ[Int(ceil(i/n_iter_Q))]
    P=ArrP[mod(i,n_iter_P)+1]
    if T[i]!=nothing
        push!(test,0)
        # if T[i][5]<2
        #     push!(Q_from_nowhere,Q)
        #     push!(P_from_nowhere,P)
        #     push!(Y_from_nowhere, T[i][6])
        # elseif T[i][5]==2
            push!(Q_init,Q)
            push!(P_init,P)
            push!(Y_initS,T[i][6])
            # push!(Y_initB,T[i][7])

            q=T[i][1]
            p=T[i][2]
            y=Yfind(q,p,.25)
            # if ~isempty(y) && ~isempty(yB)
            # y_min=min(abs(yS),abs(yB))
            if ~isempty(y) && T[i][5]==2 && abs(T[i][4]-y[1])<1e-2
            push!(Q_1before,T[i][1])
            push!(P_1before,T[i][2])
            push!(X_1before,T[i][3])
            push!(Y_1before,T[i][4])
            # push!(Q_check,T[i][6])
            # push!(P_check,T[i][7])
            # push!(X_check,T[i][8])
            # push!(Y_check,T[i][9])
            push!(hit_check1,T[i][5])
             # push!(hit_check2,T[i][6])
         end
    end
end

figure()
 axis([ -3, 3,-1.5, 1.5])
Q_exit=vcat(Q_init,Q_init,-Q_init,-Q_init)
P_exit=vcat(P_init,-P_init,P_init,-P_init)
# plot(P_exit,Q_exit ,color="r",".",markersize=2, markeredgewidth=.1)
plot(P_init,Q_init ,color="r",".",markersize=2, markeredgewidth=.1)

figure()
axis([ -3, 3,-1.5, 1.5])
Q_previous=vcat(Q_1before,Q_1before,-Q_1before,-Q_1before)
P_previous=vcat(P_1before,-P_1before,P_1before,-P_1before)
# plot(P_previous,Q_previous ,color="b",".",markersize=2, markeredgewidth=.1)
plot(P_1before,Q_1before ,color="b",".",markersize=2, markeredgewidth=.1)


# figure()
# axis([ -3, 3,-1.5, 1.5])
# Q_nowhere=vcat(Q_from_nowhere,Q_from_nowhere,-Q_from_nowhere,-Q_from_nowhere)
# P_nowhere=vcat(P_from_nowhere,-P_from_nowhere,P_from_nowhere,-P_from_nowhere)
# plot(P_nowhere,Q_nowhere ,color="r",".",markersize=2, markeredgewidth=.1)
#
# figure()
# axis([ -3, 3,-1.5, 1.5])
# Qcheck=vcat(Q_check,Q_check,-Q_check,-Q_check)
# Pcheck=vcat(P_check,-P_check,P_check,-P_check)
# plot(Pcheck,Qcheck ,color="k",".",markersize=2, markeredgewidth=.1)

figure()
axis([ -3, 3,-1.5, 1.5])
# plot(P_exit,Q_exit ,color="r",".",markersize=2, markeredgewidth=.1)
# plot(P_previous,Q_previous ,color="b",".",markersize=2, markeredgewidth=.1)
#plot(P_nowhere,Q_nowhere ,color="r",".",markersize=2, markeredgewidth=.1)
# plot(Pcheck,Qcheck ,color="k",".",markersize=2, markeredgewidth=.1)
plot(P_init,Q_init ,color="r",".",markersize=2, markeredgewidth=.1)
plot(P_1before,Q_1before ,color="b",".",markersize=2, markeredgewidth=.1)

#
# outfile_0 = "ParVern9_outfile0_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_1 = "ParVern9_outfile1_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_2 = "ParVern9_outfile2_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_3 = "ParVern9_outfile3_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_4 = "ParVern9_outfile4_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_5 = "ParVern9_outfile5_"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
# outfile_bound = "ParVern9_outfile1_bound"*string(Energy)*"_nP_"*string(n_iter_P)*"_nQ_"*string(n_iter_Q)*".dat"
#
#
# a = open(outfile_0, "w")
# b = open(outfile_1, "w")
# c = open(outfile_2, "w")
# d = open(outfile_3, "w")
# e = open(outfile_4, "w")
# f = open(outfile_5, "w")
# g = open(outfile_bound, "w")
#
# writedlm(outfile_0, T0)
# writedlm(outfile_1, T1)
# writedlm(outfile_2, T2)
# writedlm(outfile_3, T3)
# writedlm(outfile_4, T4)
# writedlm(outfile_5, T5)
# writedlm(outfile_bound, Tbound)
#
#
# close(a)
# close(b)
# close(c)
# close(d)
# close(e)
# close(f)
# close(g)
