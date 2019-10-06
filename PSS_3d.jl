

using DifferentialEquations
using Roots
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

include("YfindNP.jl")

Energy=.25
# Energy=0.1251;
# Energy=0.1249;
H=(2*Energy)^2
h=sqrt(H)
t_end=10000
barrier=20
max_hit=500


Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))

condition(u,t,integrator)= u[5]>p[1] || maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])])>p[2]
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

ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
    return
end
# refernce
# z=X+iP W=Q+iY . Plane Q=0, P=0
# z=u[1]+i u[2] ; w=u[3]+i u[4]

#Q=0, have P find Y


Q=.1
P=.1
Y=YfindNP(Q,P,H)
u0=zeros(5)
u0[1]=0 #X
u0[2]=P #P
u0[3]=Q #Q
u0[4]=Y[1] #Y
u0[5]=0
TEST1_0=H_test(u0)/H
p=zeros(2)
p[1]=max_hit
p[2]=barrier
prob = ODEProblem(Eq_of_M,u0,(0., t_end),p)
# t,A=solve(prob, RK4(),reltol=1e-15,abstol=1e-16,maxiters=1e15)

sol=solve(prob, Vern9(),reltol=1e-13,abstol=1e-15,maxiters=1e20,callback=cb)
figure()
# =1e-16,
# figure()
A=sol[:,:]
plot3D(A[3,:],A[2,:],A[1,:],linewidth=.5)

#testing if energy is conserved at end
uf=zeros(4)
uf[1]=A[1,end]
uf[2]=A[2,end]
uf[3]=A[3,end]
uf[4]=A[4,end]

#finding maxes for exit criteria
u1max=maximum(abs.(A[1,:]))
u2max=maximum(abs.(A[2,:]))
u3max=maximum(abs.(A[3,:]))
u4max=maximum(abs.(A[4,:]))
N=length(A[1,:])
H_v=zeros(N)
for k=1:N
    um=zeros(4)
    um[1]=A[1,k]
    um[2]=A[2,k]
    um[3]=A[3,k]
    um[4]=A[4,k]
    H_v[k]=H_test(um)/H
end


TEST1_f=H_test(uf)/H
# figure()
# plot(1:N, log.(H_v))
