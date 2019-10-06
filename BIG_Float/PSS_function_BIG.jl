
using DifferentialEquations
using ProgressMeter
include("YfindBIG.jl")

using PyCall
pygui(:qt)
using PyPlot
pygui(true)

condition(u,t,integrator)= u[5]>max_hit
function condition2(u,t,integrator) # Event when event_f(u,t) == 0
 u[1]
end
affect!(integrator) = terminate!(integrator)
function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
 end
cb2 = ContinuousCallback(condition2,affect2!,nothing,rootfind = true,   save_positions = (false,true))
cb1 = DiscreteCallback(condition,affect!)
cb=CallbackSet(cb2,cb1)
Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
ODE1(z,w)=conj(  im * w.*( BigFloat("1") ./(z.^2-w.^2)+BigFloat("1") ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( BigFloat("1") ./(w.^2-z.^2)+BigFloat("1") ./(1+z.^2) ))
function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
    return
end


function PSS_function_BIG(Q,P,Energy,t_end, max_hit)
    H=(2*Energy)^2
    Y=YfindBIG(Q,P,H)

    if ~isempty(Y) && Y>0

        u0=zeros(BigFloat ,5)
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y #Y
        #u0[1]=u0[5]=0

        prob = ODEProblem(Eq_of_M,u0,(BigFloat("0."), t_end))
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-30,abstol=1e-28,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-10,abstol=1e-12,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol=solve(prob, RK4(),reltol=1e-8,abstol=1e-10,callback=cb,save_start=true,save_end=true,save_everystep=false)

        uf=zeros(BigFloat ,5)
        uf[1]=sol[1,end] #X
        uf[2]=sol[2,end] #P
        uf[3]=sol[3,end] #Q
        uf[4]=sol[4,end] #Y
        uf[5]=0
        dH=abs(H_test(uf)-H)
        #
        if dH<1e-10
            return Y, sol[:,2:end-1][3,:],sol[:,2:end-1][2,:]
        end

    end

end
