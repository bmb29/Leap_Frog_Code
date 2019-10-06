using IntervalArithmetic
using DifferentialEquations
using Plots

function SHO(du,u,p,t)
du[1]=u[2]+@interval(pi)
du[2]=-u[1]
end

T_end=2*pi
u0=[1.,0.]
tspan=(0.,T_end)
prob = ODEProblem(SHO,u0,tspan)
t,u=solve(prob, RK4(),reltol=1e-8,abstol=1e-10,maxiters=1e15)
# plot(u[1,:],u[2,:])
diff=u[1,end]-u[1,1]
println(diff)
