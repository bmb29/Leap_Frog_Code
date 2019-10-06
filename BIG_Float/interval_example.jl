# using IntervalArithmetic
using DifferentialEquations
using LinearAlgebra

function Flo(du,u,p,t)
    h=BigFloat(".125")
    A11=(-1) *(1+4. *h^2+4. *h*cos(2. *t)).^(-1/2).*sin(2. *t)
    A12=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+ 4. *h*cos(2. *t)).^(1/2))).^(-1) *((-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))+h*(1+(-4).*h+16. *h^2+2. *(1+ 4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-8).*h*(1+4. *h^2+4. *h*cos(2. * t)).^(1/2)+(-1) *cos(4. *t)))
    A21=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t) ).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))).^(-1) *(h+4. *h^2. *(1+4. *h)+2. *h*(1+4. *h).*(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos( 2. *t)).^(1/2))+(-1) *h*cos(4. *t))
    A22=(1+4. *h^2+4. *h*cos(2. *t)).^( .-1/2).*sin(2. *t)

    du[1]=A11*u[1]+A12*u[2]
    du[2]=A21*u[1]+A22*u[2]
end
#

T_end=BigFloat(pi)
tspan=(BigFloat("0"),T_end)
relT=BigFloat("1e-72")
absT=relT*1e-2

u0=[BigFloat("1"),BigFloat("0")]
prob = ODEProblem(Flo,u0,tspan)
t,u1=solve(prob, Feagin14(),reltol=relT,abstol=absT,maxiters=1e20)
diff=sqrt((u1[1,end]-u1[1,1])^2+(u1[2,end]-u1[2,1])^2)
println(diff)
#
# u0=[BigFloat("0"),BigFloat("1")]
# prob = ODEProblem(Flo,u0,tspan)
# t,u2=solve(prob, RK4(),reltol=relT,abstol=absT,maxiters=1e15)
#
# monodromy_matrix=[u1[1,end] u2[1,end]; u1[2,end] u2[2,end]]
# discriminant=abs(tr(monodromy_matrix))-2
# println(discriminant)
# Flo1=eigvals(monodromy_matrix)[1]
# Flo2=eigvals(monodromy_matrix)[2]
