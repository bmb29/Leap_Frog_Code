

using DifferentialEquations
using Roots
using PyCall
pygui(:qt)
using PyPlot
pygui(true)


function condition(u,t,integrator) # Event when event_f(u,t) == 0
 u[1]
end

function affect!(integrator)
end
#cb1 = SavingCallback((u,t,integrator)->(u[4],u[2]), saved_values)

cb2= ContinuousCallback(condition, affect!)
Energy=0.112012;
H=(2*Energy)^2
h=sqrt(H)
n_iter=100
t_end=1000

# H1(p,q,t)= ( (q[2]-q[1])^2+(p[2]-p[1])^2 )*( (q[2]+q[1])^2+(p[2]+p[1])^2 )
# H2(p,q,t)= (p[2]^4+2*p[2]^2*(q[1]^2-1)+(1+q[1]^2)^2 )*(q[2]^4+2*q[2]^2*(p[1]^2+1)+(p[1]^2-1)^2 )
# HAM(p,q,t)=H1(p,q,t)/H2(p,q,t)
Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))


ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
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

#Q=0, have P find Y
P_start=8e-3
P_end=1e-1
P=range(P_start,stop=P_end,length=n_iter)
Y=zeros(n_iter)
TEST1_0=zeros(n_iter)
TEST1_f=zeros(n_iter)
u1max=zeros(n_iter)
u2max=zeros(n_iter)
u3max=zeros(n_iter)
u4max=zeros(n_iter)
YofP_find(PP)=sqrt((h+PP^2-h*PP^2 )./(1+h-h*PP^2))
axis([ -.2,.2,-.02, .02]);

# figure(dpi=3600)
for k=1:n_iter
    #  u0=zeros(4)
    # Y[k]=YofP_find(P[k])
    Pval=P[k]
    Y_find1(y)=Hamil(0,y,0,Pval)-H
    u0=zeros(4)
    Y[k]=find_zero(Y_find1,.01, maxeval=100, tol=1e-15)
    u0[1]=0 #X
    u0[2]=P[k] #P
    u0[3]=0 #Q
    u0[4]=Y[k] #Y
    TEST1_0[k]=H_test(u0)/H
    prob = ODEProblem(Eq_of_M,u0,(0., t_end))
    t,A=solve(prob, Vern9(),reltol=1e-15,abstol=1e-12,maxiters=1e15, callback=cb2,save_end=false,save_everystep=false)
    plot(A[3,:],A[2,:],"k,")
    plot(A[3,:],-A[2,:],"k,")
    plot(-A[3,:],A[2,:],"k,")
    plot(-A[3,:],-A[2,:],"k,")

    #testing if energy is conserved at end
    uf=zeros(4)
    uf[1]=A[1,end]
    uf[2]=A[2,end]
    uf[3]=A[3,end]
    uf[4]=A[4,end]
    TEST1_f[k]=H_test(uf)/H

    #finding maxes for exit criteria
    u1max[k]=



 end
#P=0, have Q find Y

Q_start=5e-3
Q_end=1e-1
Q=range(Q_start,stop=Q_end,length=n_iter)
Y=zeros(n_iter)
TEST2=zeros(n_iter)
for k=1:n_iter
    Y_find(YY)=Hamil(0,YY,Q[k],0)-H
    u0=zeros(4)
    Y[k]=find_zero(Y_find,.01, maxeval=100,tol =1e-15)
    u0[1]=0 #X
    u0[2]=0 #P
    u0[3]=Q[k] #Q
    u0[4]=Y[k] #Y
    TEST2[k]=H_test(u0)
    prob = ODEProblem(Eq_of_M,u0,(0., t_end))
    t,A=solve(prob, Vern9(),reltol=1e-15,abstol=1e-12,maxiters=1e15, callback=cb2,save_end=false,save_everystep=false)
    plot(A[3,:],A[2,:],"k,")
    plot(A[3,:],-A[2,:],"k,")
    plot(-A[3,:],A[2,:],"k,")
    plot(-A[3,:],-A[2,:],"k,")
end
# savefig("PSS11.png")
