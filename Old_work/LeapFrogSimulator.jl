using DifferentialEquations
# using Plots
include("YfindNP.jl")
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

H=.25
h=1/(2*H)
max_time=50
a=(2+h-2*sqrt(h+1))/h; # paramter a (ratio) in terms of h

#exit
# Q=0.0015000000999000002
# P=0.0411
# Y=0.5780004334287135
#before
# Q=0.7040082229574488
# P=-1.798419378432707
# Y=1.523090142608121
#nowhere
Q=-1.0010302468906507
P=1.7385260494565795
Y=YfindNP(Q,P,H)

u0=[0;P;Q;Y[1];0;0;0]
# u0= [0;0;pert;-Y_0;0;0;0]
tspan = (0.0,max_time)
#10.531140269496431
#4.75752093999673
ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
Centroid_ODE(z,w)=conj( 1 ./(1+z.^2)+1 ./(1+w.^2) )
# refernce
# z=X+iP W=Q+iY . Plane Q=0, P=0
# z=u[1]+i u[2] ; w=u[3]+i u[4]

function Eq_of_M(du,u,p,t)
du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
du[5]=real(Centroid_ODE(u[1]+im*u[2],u[3]+im*u[4]));
du[6]=imag(Centroid_ODE(u[1]+im*u[2],u[3]+im*u[4]))
du[7]=0
return
end
condition(u,t,integrator)= u[7]==20
affect!(integrator) = terminate!(integrator)
condition2(u,t,integrator)=  u[1]
function affect2!(integrator)
integrator.u[7]=integrator.u[7]+1
end
cb1 = DiscreteCallback(condition,affect!)
cb2 =ContinuousCallback(condition2,affect2!,rootfind=true)
cb=CallbackSet(cb2,cb1)

prob = ODEProblem(Eq_of_M,u0,tspan)
t,A=solve(prob, Vern9(),reltol=1e-13,abstol=1e-20,maxiters=1e15,callback=cb)
N=length(t.t)
Z=A[1,:]+im *A[2,:];
W=A[3,:]+im *A[4,:];
Energy=zeros(N)
for k=1:N
    Energy[k]=abs(1/(1+Z[k]^2)-1/(1+W[k]^2))
end


Centeroid=zeros(N);
Centeroid=A[5,:]+im *A[6,:];
lin_impulse=im*(1+a )*ones(N);
z_1pos=.5*(Centeroid+Z+lin_impulse+W);
z_1neg=.5*(Centeroid+Z-lin_impulse-W);
z_2pos=.5*(Centeroid-Z+lin_impulse-W);
z_2neg=.5*(Centeroid-Z-lin_impulse+W);
DDD=(z_2pos+z_1pos-z_2neg-z_1neg)/2;

z_1posX=real(z_1pos)
z_1posY=imag(z_1pos)

z_2posX=real(z_2pos)
z_2posY=imag(z_2pos)

z_1negX=real(z_1neg)
z_1negY=imag(z_1neg)

z_2negX=real(z_2neg)
z_2negY=imag(z_2neg)

figure()
plot(t.t,A[1,:])
uf=zeros(7)
uf[1]=A[1,end]#X
uf[2]=A[2,end] #P
uf[3]=A[3,end] #Q
uf[4]=A[4,end] #Y
uf[5]=A[5,end] #Y
uf[6]=A[6,end]
#
# figure()
# tspan = (0.0,-max_time)
# prob = ODEProblem(Eq_of_M,uf,tspan)
# conditionM(u,t,integrator)= u[7]==20
# cb1 = DiscreteCallback(conditionM,affect!)
# cb2 =ContinuousCallback(condition2, affect2!, rootfind=true)
# cb=CallbackSet(cb2,cb1)
# t,B=solve(prob, Vern9(),reltol=1e-13,abstol=1e-20,maxiters=1e15,callback=cb)
# plot(t.t,B[1,:])
#
#
#
# figure()
# tspan = (0.0,-max_time)
# prob = ODEProblem(Eq_of_M,u0,tspan)
# conditionM(u,t,integrator)= u[7]==20
#
# cb1 = DiscreteCallback(conditionM,affect!)
# cb2 =ContinuousCallback(condition2, affect2!, rootfind=true)
# cb=CallbackSet(cb2,cb1)
# t,B=solve(prob, Vern9(),reltol=1e-13,abstol=1e-20,maxiters=1e15,callback=cb)
# plot(t.t,B[1,:])

# # col=rand(1)
# plot!(z_1pos,linewidth=1,seriescolor=:red)
# plot!(z_1neg, linewidth=1,seriescolor=:blue)
# plot!(z_2pos,linewidth=1)
# plot!(z_2neg,linewidth=1)
