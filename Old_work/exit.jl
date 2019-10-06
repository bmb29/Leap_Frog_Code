
using DifferentialEquations
using Plots
barrier=4

condition(u,t,integrator)=abs(u[1])>barrier || abs(u[2])>barrier || abs(u[3])>barrier||abs(u[4])>barrier
affect!(integrator) = terminate!(integrator)
cb=DiscreteCallback(condition,affect!)
Energy=0.18;
H=(2*Energy)^2
h=sqrt(H)


PofY_find(YY)=sqrt((h+YY.^2-h*YY.^2 )./(1+h-h*YY.^2))
H1(p,q,t)= ( (q[2]-q[1])^2+(p[2]-p[1])^2 )*( (q[2]+q[1])^2+(p[2]+p[1])^2 )
H2(p,q,t)= (p[2]^4+2*p[2]^2*(q[1]^2-1)+(1+q[1]^2)^2 )*(q[2]^4+2*q[2]^2*(p[1]^2+1)+(p[1]^2-1)^2 )
HAM(p,q,t)=H1(p,q,t)/H2(p,q,t)
PofY_find(YY)=sqrt((h+YY.^2-h*YY.^2 )./(1+h-h*YY.^2))

n_iter=100
Y_start=1e-6
Y_end=1e-2
t_end=1e10

Y=range(Y_start,stop=Y_end,length=n_iter)
P=zeros(n_iter)
T=zeros(n_iter)

for k=1:n_iter
    global
    q0,p0=[zeros(2) for i in 1:2]
    P[k]=PofY_find(Y[k])
    q0[1]=0
    p0[1]=P[k]
    q0[2]=0
    p0[2]=Y[k]
    prob= HamiltonianProblem{true}(HAM, p0, q0, (0., t_end))
    t,A=solve(prob, Vern9(),reltol=1e-13,abstol=1e-20,maxiters=1e15, callback=cb)
    T[k]=t.t[end]
end
println("Max exit time is  ", maximum(T))
