#
include("Yfind.jl")
include("is_it_exit2.jl")
include("is_it.jl")

@everywhere function BASIN_PSS_function_parallel(K)


    Q_start=-1.5
    Q_end=1.5
    P_start=-2.5
    P_end=2.5
    n_iter_P=2000
    n_iter_Q=2001

    delta_P=(P_end-P_start)/n_iter_P
    delta_Q=(Q_end-Q_start)/n_iter_Q

    ArrP=P_start: delta_P:P_end
    ArrQ=Q_start: delta_Q: Q_end


    Q=ArrQ[Int(ceil(K/n_iter_Q))]
    P=ArrP[mod(K,n_iter_P)+1]

    max_hit=5
    Energy=.25
    t_end=1000
    tol_dist=1e-5
    H=(2*Energy)^2

    Q_0=0
    P_0=0

    Q_1=0
    P_1=0

    Q_2=0
    P_2=0

    Q_3=0
    P_3=0


    Y=Yfind(Q,P,H)
    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y[1] #Y
        u0[5]=0 #Y
        # counter variable
        #are the initial conditions right?
        # prob_test= ODEProblem(Eq_of_M,u0,(0., 10.))
        # sol_test=solve(prob_test, RK4(),maxiters=1e10 )
        # S=sign(sol_test.u[2][1])
        if true
            prob = ODEProblem(Eq_of_M,u0,(0., t_end))
            # sol=solve(prob,RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb)
            sol=solve(prob, Vern9(),maxiters=1e10, reltol=1e-8,abstol=1e-10,callback=cb)
            uf=zeros(5)
            uf[1]=sol[1,end] #X
            uf[2]=sol[2,end] #P
            uf[3]=sol[3,end] #Q
            uf[4]=sol[4,end] #Y
            uf[5]=0
            dH=abs(H_test(uf)-H)
            #
            if dH<1e-5
                B=hcat(sol.u[end-10],sol.u[end])'

                if is_it_exit2(B,H,tol_dist) && is_it_exit(sol.u[end],H,tol_dist)
                    return sol.u[end][5]/2.0
                else
                    return -1
                end
            end
        end
    end
    return -2
end
