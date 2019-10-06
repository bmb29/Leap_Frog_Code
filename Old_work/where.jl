using DelimitedFiles
using DifferentialEquations
using ProgressMeter
include("YfindNP.jl")
include("is_it_exit.jl")
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
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
condition(u,t,integrator)= u[5]==max_hit
function condition2(u,t,integrator) # Event when event_f(u,t) == 0
u[1]
end
condition3b(u,t,integrator)=  u[5]==2
condition3f(u,t,integrator)=  u[5]==3

affect!(integrator) = terminate!(integrator)
function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
end
# @everywhere cb2 = ContinuousCallback(condition2,affect2!, save_positions = (false,true) )
cb2 = ContinuousCallback(condition2, affect2!, rootfind=true, save_positions = (true,true))
cb2f = ContinuousCallback(condition2,  affect2!, rootfind=true, save_positions = (true,true))
cb2b = ContinuousCallback(condition2, affect2!,rootfind=true)

cb1 = DiscreteCallback(condition,affect!)
cb3b = DiscreteCallback(condition3b,affect!)
cb3f = DiscreteCallback(condition3f,affect!)

cb=CallbackSet(cb2,cb1)

cb_backward=CallbackSet(cb2b,cb3b)
cb_forward=CallbackSet(cb2f,cb3f)
max_hit=4
Energy=.25
t_end=1e5
tol_dist=1e-3
H=(2*Energy)^2

Q_start=1e-10
Q_end=1.5
P_start=0
P_end=3

n_iter_P=501
n_iter_Q=501

delta_P=(P_end-P_start)/n_iter_P
delta_Q=(Q_end-Q_start)/n_iter_Q

ArrP= P_start: delta_P: P_end
ArrQ= Q_start: delta_Q: Q_end
#up
Q= 0.4275000000715
P= 1.53
#down
Q_check=zeros(0)
P_check=zeros(0)

figure()
@showprogress 1 "Computing..." for i=1:length(ArrP)
    for j=1:length(ArrQ)
        P=ArrP[i]
        Q=ArrQ[j]
        Y=YfindNP(Q,P,H)
if ~isempty(Y) && Y>0
    u0=zeros(5)
    u0[1]=0 #X=0
    u0[2]=P #P
    u0[3]=Q #Q
    u0[4]=Y #Y
    u0[5]=0 # counter variable
    #are the initial conditions right?
    # prob_test= ODEProblem(Eq_of_M,u0,(0., 100.))
    # sol_test=solve(prob_test, RK4(),maxiters=1e10 )
    if true #sol_test.u[2][1]>0 #check if upcrossing
        #THESE GUYS HIT ONCE AND LEAVE

        prob_forward = ODEProblem(Eq_of_M,u0,(0., t_end))
        sol_forward=solve(prob_forward, Vern9(),maxiters=1e20, reltol=1e-10,abstol=1e-12,callback=cb, save_start=true ,save_end=true,save_everystep=true)
        # plot(sol_forward.t,sol_forward[1,:])
        B=hcat(sol_forward.u[end-10],sol_forward.u[end])'
        if is_it_exit(B,H,tol_dist) && sol_forward.u[end][5]==2
         #THESE GUYS HIT ONCE AND LEAVE
            prob_backward = ODEProblem(Eq_of_M,u0,(0., -t_end))
            sol_backward=solve(prob_backward, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_backward,save_start=true,save_end=true,save_everystep=true)
            uf=zeros(5)
            uf[1]=0 #X
            uf[2]=sol_backward[2,end] #P
            uf[3]=sol_backward[3,end] #Q
            uf[4]=sol_backward[4,end] #Y
            uf[5]=0
            dH=abs(H_test(uf)-H)
            prob_check= ODEProblem(Eq_of_M,uf,(0., t_end))
            sol_check=solve(prob_check, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_forward,save_start=true,save_end=true,save_everystep=true)
            if sol_backward.u[end][5]==2
                push!(Q_check, Q)
                push!(P_check, P)
                plot(sol_check.t,sol_check[1,:])
                # plot(sol_backward.t,sol_backward[1,:])
            end
        end
  end
end
end
end
