using IntervalArithmetic
using DifferentialEquations
using LinearAlgebra
using TaylorIntegration
using DiffEqBase
# setprecision(100)
function Flo(du,u,p,t)
    h=@interval(.125)
    A11=(-1) *(1+4. *h^2+4. *h*cos(2. *t)).^(-1/2).*sin(2. *t)
    A12=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+ 4. *h*cos(2. *t)).^(1/2))).^(-1) *((-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))+h*(1+(-4).*h+16. *h^2+2. *(1+ 4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-8).*h*(1+4. *h^2+4. *h*cos(2. * t)).^(1/2)+(-1) *cos(4. *t)))
    A21=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t) ).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))).^(-1) *(h+4. *h^2. *(1+4. *h)+2. *h*(1+4. *h).*(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos( 2. *t)).^(1/2))+(-1) *h*cos(4. *t))
    A22=(1+4. *h^2+4. *h*cos(2. *t)).^( .-1/2).*sin(2. *t)

    du[1]=A11*u[1]+A12*u[2]
    du[2]=A21*u[1]+A22*u[2]
end
#

T_end=pi
tspan=(0.,T_end)
tol=1e-16

u0=[1.,0.]
prob = ODEProblem(Flo,u0,tspan)

# t,u=taylorinteg(Flo,u0, tspan[1],tspan[2],25,1e-50, maxsteps=1_000_000_000)
t,u=solve(prob,TaylorMethod(10),abstol=tol,maxiters=1_000_000_000)
diff=sqrt((u[1,end]-u0[1])^2+(u[2,end]-u0[2])^2)
# diff=sqrt((u[1,end]-u1[1,1])^2+(u1[2,end]-u1[2,1])^2)
# println(diff)
# using IntervalArithmetic
# using DifferentialEquations
# using LinearAlgebra
# using TaylorIntegration
# using DiffEqBase
# setprecision(100)
# @taylorize function Flo(du,u,p,t)
#     h=@interval(.125)
#     A11=(-1) *(1+4. *h^2+4. *h*cos(2. *t)).^(-1/2).*sin(2. *t)
#     A12=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+ 4. *h*cos(2. *t)).^(1/2))).^(-1) *((-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))+h*(1+(-4).*h+16. *h^2+2. *(1+ 4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-8).*h*(1+4. *h^2+4. *h*cos(2. * t)).^(1/2)+(-1) *cos(4. *t)))
#     A21=(1+4. *h^2+(1+4. *h^2+4. *h*cos(2. *t) ).^(1/2)+2. *h*cos(2. *t).*(2+(1+4. *h^2+4. *h*cos(2. *t)).^(1/2))).^(-1) *(h+4. *h^2. *(1+4. *h)+2. *h*(1+4. *h).*(1+4. *h^2+4. *h*cos(2. *t)).^(1/2)+(-1) *cos(2. *t).*(1+(-8).*h^2+(1+4. *h^2+4. *h*cos( 2. *t)).^(1/2))+(-1) *h*cos(4. *t))
#     A22=(1+4. *h^2+4. *h*cos(2. *t)).^( .-1/2).*sin(2. *t)
#
#     du[1]=A11*u[1]+A12*u[2]
#     du[2]=A21*u[1]+A22*u[2]
# end
# #
#
# T_end=BigFloat(pi)
# tspan=(BigFloat("0"),T_end)
# tol=BigFloat("1e-25")
#
# u0=[BigFloat("1"),BigFloat("0")]
# prob = ODEProblem(Flo,u0,tspan)
#
# # t,u=taylorinteg(Flo,u0, tspan[1],tspan[2],25,1e-50, maxsteps=1_000_000_000)
# t,u=solve(prob,TaylorMethod(25),abstol=tol,maxiters=1_000_000_000_000)
# diff=sqrt((u[1,end]-u0[1])^2+(u[2,end]-u0[2])^2)
# # diff=sqrt((u[1,end]-u1[1,1])^2+(u1[2,end]-u1[2,1])^2)
# println(diff)
