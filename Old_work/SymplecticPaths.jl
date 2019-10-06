using DifferentialEquations
using Plots
using LinearAlgebra



function DeQ(du,u,p,t)
du[1]=-u[2]*(1+u[1]^2) / ( (u[1]^2+u[2]^2)*(1-u[2]^2));
du[2]=u[1]*(1-u[2]^2)/ ( (u[1]^2+u[2]^2)*(1+u[1]^2));
du[3]= 2*u[1]*u[2]/(u[1]^2+u[2]^2)^2 *u[3]-1/(u[2]^2+u[1]^2)^2* (1-u[2]^2)*(3*u[1]^4+u[2]^2*u[1]^2-u[2]^2+u[1]^2) /( (1+u[1]^2)^2 )*u[4];
du[4]=-(1/(u[1]^2+u[2]^2)^2* (1+u[1]^2)*(3*u[2]^4+u[1]^2*u[2]^2+u[1]^2-u[2]^2)/( (1-u[2]^2)^2) )*u[3]-(2*u[1]*u[2]/(u[1]^2+u[2]^2)^2)*u[4];
du[5]= 2*u[1]*u[2]/(u[1]^2+u[2]^2)^2 *u[5]-1/(u[2]^2+u[1]^2)^2* (1-u[2]^2)*(3*u[1]^4+u[2]^2*u[1]^2-u[2]^2+u[1]^2) /( (1+u[1]^2)^2 )*u[6];
du[6]=-(1/(u[1]^2+u[2]^2)^2* (1+u[1]^2)*(3*u[2]^4+u[1]^2*u[2]^2+u[1]^2-u[2]^2)/( (1-u[2]^2)^2) )*u[5]-(2*u[1]*u[2]/(u[1]^2+u[2]^2)^2)*u[6];
end
condition(u,t,integrator)= u[1]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)


#
H=range(.05,stop=.2,length=10)
# for k=1:length(H)
H=.125
    h=1/(2*H)
    Y0=1.0/sqrt(h+1)
    tspan0=(0.,1e4)
    u0=[0.,Y0,0.,0.,0.,0.]
    prob = ODEProblem(DeQ,u0,tspan0)
    t,A1=solve(prob, Vern9(),reltol=1e-12,abstol=1e-20,maxiters=1e15,callback=cb)
    FinalT=t.t[end]


    u0=[0.,Y0,1.,0.,0.,1.]
    tspan=(0.,FinalT)
    prob = ODEProblem(DeQ,u0,tspan)
    t,A=solve(prob, Vern9(),reltol=1e-12,abstol=1e-12,maxiters=1e15)
    N=length(t.t)
    r=zeros(N)
    z=zeros(N)
    theta=zeros(N)
    for j=1:N
        #Create Monodromy matrix
        a=A[3,j]
        b=A[4,j]
        c=A[5,j]
        d=A[6,j]
        M=[a c; b d]
        W,S,V = svd(M)
        P=W*Diagonal(S)*W'
        U=W*V'
        r[j]=P[1,1]
        z[j]=P[1,2]
        theta[j]=asin(U[2,1])
    end
    x=r.*cos.(theta)
    y=r.*sin.(theta)


plot3d(x,y,z)
