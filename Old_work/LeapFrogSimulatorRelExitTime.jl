using DifferentialEquations
using Plots
n_iter=50
barrier=4
HH=range(.125001,stop=.499,length=n_iter)
exit_time=zeros(length(HH))
dE=zeros(length(HH))
max_time=exp(8)
tspan = (0.0,max_time)
pert=1e-6

ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
condition(u,t,integrator)=abs(u[1])>barrier || abs(u[2])>barrier || abs(u[3])>barrier||abs(u[4])>barrier
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition,affect!)
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

for j=1:length(HH)
    H=HH[j]
    h=1/(2*H)
    Y_0=1/sqrt(h+1)
    u0= [0;0;pert;Y_0;0;0]
    prob = ODEProblem(Eq_of_M,u0,tspan)
    t,A=solve(prob, RK4(),reltol=1e-8,abstol=1e-10,maxiters=1e15,callback=cb)
    exit_time[j]=t.t[end]
    Z=A[1,:]+im *A[2,:];
    W=A[3,:]+im *A[4,:];
    N=length(t.t)
    Energy=zeros(N)
    for k=1:N
        Energy[k]=abs(1/(1+Z[k]^2)-1/(1+W[k]^2))
    end
    dE[j]=Energy[end]-Energy[1]

end

p=plot(HH,log.(exit_time),seriestype=:scatter,marker=:circle,msw = 0,ms=2,seriescolor = :black,legend=false)
savefig(p,"ExitLogLog1000b.pdf")
 plot(-HH,-log.(abs.(dE)))
