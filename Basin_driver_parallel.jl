using Distributed
# @everywhere using DelimitedFiles
@everywhere using DifferentialEquations
@everywhere using ProgressMeter
@everywhere using MATLAB
@everywhere using Printf

include("BASIN_PSS_function_parallel.jl")

# @everywhere using PyCall
# pygui(:qt)
# @everywhere using PyPlot
# pygui(true)

@everywhere max_hit=5
@everywhere Energy=.25
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
@everywhere cb2 = ContinuousCallback(condition2,affect2!,rootfind=false)
@everywhere cb1 = DiscreteCallback(condition,affect!)

@everywhere cb=CallbackSet(cb2,cb1)
@everywhere H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
@everywhere ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
@everywhere ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
@everywhere function Eq_of_M(du,u,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
return
end

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

# ArrQ[k]=.1*rand()
# ArrP[j]=.01*rand()
@everywhere Q_start=-1.5
@everywhere Q_end=1.5
@everywhere P_start=-2.5
@everywhere P_end=2.5
@everywhere n_iter_P=2000
@everywhere n_iter_Q=2001
@everywhere N=n_iter_P*n_iter_Q

@everywhere delta_P=(P_end-P_start)/n_iter_P
@everywhere delta_Q=(Q_end-Q_start)/n_iter_Q

@everywhere ArrP=P_start: delta_P:P_end
@everywhere ArrQ=Q_start: delta_Q: Q_end

Energy=.25
location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
h=replace(@sprintf("%.13f",Energy),"."=>"_")
file_name=location*"Escape_"*h*".fig"


T=@showprogress pmap(BASIN_PSS_function_parallel,1:N)
for i=1:N
    Q=ArrQ[Int(ceil(i/n_iter_Q))]
    P=ArrP[mod(i,n_iter_P)+1]
    if floor(T[i])==0
        push!(Q_0,Q)
        push!(P_0,P)
    elseif floor(T[i])==1
        push!(Q_1,Q)
        push!(P_1,P)
    elseif floor(T[i])==2
        push!(Q_2,Q)
        push!(P_2,P)
    elseif floor(T[i])==3
        push!(Q_3,Q)
        push!(P_3,P)
    elseif floor(T[i])==4
        push!(Q_4,Q)
        push!(P_4,P)
    elseif floor(T[i])==5
        push!(Q_5,Q)
        push!(P_5,P)
    elseif T[i]==-1

        push!(Q_bound,Q)
        push!(P_bound,P)
    end
end

mat"figure(); hold on;"
mat"axis([ -$width,$width,-$height,$height ])"
# mat"plot($P_bound,$Q_bound ,'k.','MarkerSize',10)"
mat"plot($P_0,$Q_0 ,'b.','MarkerSize',$size)"
mat"plot($P_1,$Q_1 ,'r.','MarkerSize',$size)"
mat"plot($P_2,$Q_2 ,'g.','MarkerSize',$size)"
mat"plot($P_4,$Q_4 ,'c.','MarkerSize',$size)"
mat"plot($P_5,$Q_5 ,'y.','MarkerSize',$size)"
mat"savefig($file_name)"
mat"close"
