

using DifferentialEquations
using Roots
using PyCall
include("YfindBIG.jl")
pygui(:qt)
using PyPlot
pygui(true)




Energy=BigFloat("0.1249")
H=(2*Energy)^2
t_end=100

barrier=20
Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
function g(resid,u,p,t)
  resid[1] = ( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))-H
  resid[2] = 0
  resid[3] = 0
  resid[4] = 0
end
condition(u,t,integrator)=abs(u[1])>barrier || abs(u[2])>barrier || abs(u[3])>barrier||abs(u[4])>barrier
affect!(integrator) = terminate!(integrator)
cb1 = DiscreteCallback(condition,affect!)

# Restrict to energy surface
cb2 = ManifoldProjection(g)
cb=CallbackSet(cb1)

ODE1(z,w)=conj(  im * w.*( BigFloat("1") ./(z.^2-w.^2)+BigFloat("1") ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( BigFloat("1") ./(w.^2-z.^2)+BigFloat("1") ./(1+z.^2) ))
function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    return
end
# refernce
# z=X+iP W=Q+iY . Plane Q=0, P=0
# z=u[1]+i u[2] ; w=u[3]+i u[4]

#Q=0, have P find Y


Q=BigFloat(".1")
P=BigFloat(".1")
Y=YfindBIG(Q,P,H)
u0=zeros(4)
u0[1]=0 #X
u0[2]=P #P
u0[3]=Q #Q
u0[4]=Y #Y
TEST1_0=H_test(u0)/H
prob = ODEProblem(Eq_of_M,u0,(0., t_end))
# t,A=solve(prob, RK4(),reltol=1e-15,abstol=1e-16,maxiters=1e15)

sol=solve(prob, Vern9(),maxiters=BigInt(1e10), reltol=BigFloat("1e-15"),abstol=BigFloat("1e-17"),callback=cb)
# figure()

A=sol[:,:]
A=convert(Array{Float64},A)
figure()
 plot3D(A[3,:],A[2,:],A[1,:],linewidth=.5)




N=length(A[1,:])
H_v=zeros(BigFloat,N)
umax=zeros(BigFloat ,length(A[1,:]))
dH=zeros(BigFloat ,length(A[1,:]))
for i=1:length(A[1,:])
    uf=zeros(BigFloat ,4)
    uf[1]=A[1,i] #X
    uf[2]=A[2,i] #P
    uf[3]=A[3,i]#Q
    uf[4]=A[4,i] #Y
    dH[i]=abs(H_test(uf)/H)
    umax[i]=maximum(uf)
end
MAX=maximum(umax)
worse=log(maximum(dH))
dH=convert(Array{Float64}, dH)

figure()
plot(1:N, log.(dH))
