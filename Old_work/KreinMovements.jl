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
H=.1251
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
    t,A=solve(prob, Vern9(),reltol=1e-15,abstol=1e-15,maxiters=1e15)
    N=length(t.t)
    theta1=zeros(N)
    theta2=zeros(N)
    for j=1:N
        #Create Monodromy matrix
        a=A[3,j]
        b=A[4,j]
        c=A[5,j]
        d=A[6,j]
        M=[a c; b d]


        Flo1=eigvals(M)[1]
        Flo2=eigvals(M)[2]


        theta1[j]=angle(Flo1/abs(Flo1))
        theta2[j]=angle(Flo2/abs(Flo2))
    end
    plot(t.t, theta1)
    # plot(cos.(theta1), sin.(theta1))
    # plot!(cos.(theta2), sin.(theta2))
