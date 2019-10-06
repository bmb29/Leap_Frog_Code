using DifferentialEquations
using Plots
include("is_it.jl")
barrier=300
H=.19
h=1/(2*H)
max_time=exp(10)
a=(2+h-2*sqrt(h+1))/h; # paramter a (ratio) in terms of h

pert=1e-3
Y_0=1/sqrt(h+1)
u0= [0;0;pert;Y_0]
tspan = (0.0,max_time)

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
return
end
condition(u,t,integrator)=abs(u[1])>barrier || abs(u[2])>barrier || abs(u[3])>barrier||abs(u[4])>barrier
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition,affect!)

prob = ODEProblem(Eq_of_M,u0,tspan)
sol=solve(prob, Vern9(),reltol=1e-13,abstol=1e-15,maxiters=1e15,callback=cb)
A=sol[:,:]
N=length(sol.t)
Z=A[1,:]+im *A[2,:];
W=A[3,:]+im *A[4,:];
Energy=zeros(N)
for k=1:N
    Energy[k]=abs(1/(1+Z[k]^2)-1/(1+W[k]^2))
end


Centeroid=zeros(N);

lin_impulse=im*(1+a )*ones(N);
z_1pos=.5*(Z+lin_impulse+W);
z_1neg=.5*(Z-lin_impulse-W);
z_2pos=.5*(-Z+lin_impulse-W);
z_2neg=.5*(-Z-lin_impulse+W);

#1pos and 1neg bound
d_11=abs.(z_1pos-z_1neg)
#1pos and 2neg bound
d_12=abs.(z_1pos-z_2neg)

z_1posX=real(z_1pos)
z_1posY=imag(z_1pos)

z_2posX=real(z_2pos)
z_2posY=imag(z_2pos)

z_1negX=real(z_1neg)
z_1negY=imag(z_1neg)

z_2negX=real(z_2neg)
z_2negY=imag(z_2neg)
gr()

p=plot()

plot!(z_1pos,linewidth=1,seriescolor=:red)
plot!(z_1neg, linewidth=1,seriescolor=:blue)
plot!(z_2pos,linewidth=1)
plot!(z_2neg,linewidth=1)
# legend('$z_1^+$', '$z_1^-$', '$z_2^+$', '$z_2^-$','Interpreter','latex')
uf=zeros(4)
uf[1]=A[1,end] #X
uf[2]=A[2,end] #P
uf[3]=A[3,end]#Q
uf[4]=A[4,end] #Y
tol=1e-1
is_it_exit(uf,H,tol)
