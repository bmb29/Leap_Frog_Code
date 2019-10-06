
using DifferentialEquations
using Plots
using Roots
barrier=4
function condition(u,t,integrator) # Event when event_f(u,t) == 0
 u[3]
end

function affect!(integrator)
end
#cb1 = SavingCallback((u,t,integrator)->(u[4],u[2]), saved_values)
cb= ContinuousCallback(condition, affect!)


    Energy=.1;
    H=(2*Energy)^2
    h=sqrt(H)
    n_iter=3
    # Y_start=rand(1)*1e-6
    # Y_end=Y_start
    t_end=2.1e3
    # Y=range(Y_start,stop=Y_end,length=n_iter)
    # X=[rand(1)*1e-6]
    # P=zeros(n_iter)
    Final_H=zeros(n_iter,n_iter)
    H1(p,q,t)= ( (q[2]-q[1])^2+(p[2]-p[1])^2 )*( (q[2]+q[1])^2+(p[2]+p[1])^2 )
    H2(p,q,t)= (p[2]^4+2*p[2]^2*(q[1]^2-1)+(1+q[1]^2)^2 )*(q[2]^4+2*q[2]^2*(p[1]^2+1)+(p[1]^2-1)^2 )
    HAM(p,q,t)=H1(p,q,t)/H2(p,q,t)
    PofY_find(YY)=sqrt((h+YY.^2-h*YY.^2 )./(1+h-h*YY.^2))
    p=plot()
    DPI_VAL=1000

    Hsolve(XX,YY,PP)=((XX.^2+(PP-YY).^2)*(XX.^2+(PP+YY).^2)) ./((PP.^4+2*PP.^2*(XX.^2-1)+(1+XX.^2).^2)*(YY.^2-1).^2)
    A_start=1e-5
    A_end=.1
    t_end=2.1e3
    global
    B=range(A_start,stop=A_end,length=n_iter)


    for j=1:n_iter
    x=B[j]
      for k=1:n_iter
        y=B[k]
        q0,p0=[zeros(2) for i in 1:2]
        F(PP)=Hsolve(x,y,PP)-H
        pp=find_zero(F,.01)
        q0[1]=x
        p0[1]=pp
        q0[2]=0
        p0[2]=y
        Final_H[j,k]=HAM(p0,q0,0)
        prob= HamiltonianProblem{true}(HAM, p0, q0, (0., t_end))
        t,A=solve(prob, Vern9(),reltol=1e-8,abstol=1e-20,maxiters=1e15, callback=cb,save_end=false,save_everystep=false)
        plot!(A[2,:],A[4,:],seriestype=:scatter,marker=:circle,msw = 0,ms=1,seriescolor = :black,legend=false,dpi=DPI_VAL)
        plot!(A[2,:],-A[4,:],seriestype=:scatter,marker=:circle,msw = 0,ms=1,seriescolor = :black,legend=false)
        plot!(-A[2,:],A[4,:],dpi=DPI_VAL,seriestype=:scatter,marker=:circle,msw = 0,ms=1,seriescolor = :black,legend=false)
        plot!(-A[2,:],-A[4,:],dpi=DPI_VAL,seriestype=:scatter,marker=:circle,msw = 0,ms=1,seriescolor = :black,legend=false)
        plot!([A[2,end]],[A[4,end]],seriestype=:scatter,marker=:circle,msw = 0,ms= 1,seriescolor = :red,legend=false,dpi=DPI_VAL)
        plot!([≈],[-A[4,end]],dpi=DPI_VAL,seriestype=:scatter,marker=:circle,msw = 0,ms=1,seriescolor = :red,legend=false)
        plot!([A[2,end]],[A[4,end]],dpi=DPI_VAL,seriestype=:scatter,marker=:circle,msw = 0,ms=1,seriescolor = :red,legend=false)
        plot!([-A[2,end]],[-A[4,end]],dpi=DPI_VAL,seriestype=:scatter,marker=:circle,msw = 0,ms=1,seriescolor =:red,legend=false)
    # XX=A[2,end]
    # YY=A[4,end]
end
end
 # img = SVG("pss.svg", 14cm, 8cm)
 # draw(img, p)
