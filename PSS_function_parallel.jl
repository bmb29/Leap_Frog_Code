
using DifferentialEquations
using ProgressMeter
include("Yfind.jl")

using PyCall
pygui(:qt)
using PyPlot
pygui(true)


@everywhere function PSS_function_parallel(K)
    Q_start=1e-10
    Q_end=1
    P_start=0
    P_end=.3

    n_iter_P=200
    n_iter_Q=200

    delta_P=(P_end-P_start)/n_iter_P
    delta_Q=(Q_end-Q_start)/n_iter_Q

    ArrP=P_start: delta_P:P_end
    ArrQ=Q_start: delta_Q: Q_end

    Q=ArrQ[Int(ceil(K/n_iter_Q))]
    P=ArrP[mod(K,n_iter_P)+1]
    max_hit=20
    Energy=.25
    t_end=1e4
    H=(2*Energy)^2

    Y=Yfind(Q,P,H)

    if ~isempty(Y) && Y>0

        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y #Y
        u0[5]=0 # counter variable
        #are the initial conditions right?
        prob = ODEProblem(Eq_of_M,u0,(0., t_end))
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-10,abstol=1e-12,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol=solve(prob, RK4(),reltol=1e-8,abstol=1e-10,callback=cb,save_start=true,save_end=true,save_everystep=false)

        uf=zeros(5)
        uf[1]=sol[1,end] #X
        uf[2]=sol[2,end] #P
        uf[3]=sol[3,end] #Q
        uf[4]=sol[4,end] #Y
        uf[5]=0
        dH=abs(H_test(uf)-H)
        #
        # prob_test= ODEProblem(Eq_of_M,u0,(0., 10.))
        # sol_test=solve(prob_test, RK4(),maxiters=1e10 )
        # S=sign(sol_test.u[2][1])
        # if dH<1e-8 && S<0
        #     return sol[:,1:end-1][3,:],sol[:,1:end-1][2,:]
        if dH<1e-8 #&& S>0
            return sol[:,2:end-1][3,:],sol[:,2:end-1][2,:]
        end
    end
end
