include("Yfind.jl")
include("is_it_exit2.jl")
include("is_it.jl")
@everywhere function where_did_it_come_from_parallel(K)
    Q_end=1.5
    Q_start=-Q_end
    P_end=3
    P_start=-P_end

    n_iter_P=8001
    n_iter_Q=8001

    delta_P=(P_end-P_start)/(n_iter_P)
    delta_Q=(Q_end-Q_start)/(n_iter_Q)

    ArrP= P_start: delta_P: P_end
    ArrQ= Q_start: delta_Q: Q_end


    Q=ArrQ[Int(ceil(K/n_iter_Q))]
    P=ArrP[mod(K,n_iter_P)+1]

    max_hit=3
    Energy=.25
    t_end=1e8
    tol_dist=1e-7
    H=(2*Energy)^2
    #
    # Q_0=0
    # P_0=0
    #
    # Q_1=0
    # P_1=0
    #
    # Q_2=0
    # P_2=0
    #
    # Q_3=0
    # P_3=0
    #

    Y=Yfind(Q,P,H)
    # YB=YfindB(Q,P,H)

    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y[1] #Y
        u0[5]=0 #Y
        # counter variable
        #are the initial conditions right?


        prob_test= ODEProblem(Eq_of_M,u0,(0., 10.))
        sol_test=solve(prob_test, RK4(),maxiters=1e10 )
        S=sign(sol_test.u[2][1])
        if true#check if upcrossing
            prob_forward = ODEProblem(Eq_of_M,u0,(0., t_end))
            sol_forward=solve(prob_forward, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb, save_start=true,save_end=true,save_everystep=true)
            B=hcat(sol_forward.u[end-5],sol_forward.u[end])'
            if is_it_exit2(B,H,tol_dist) && is_it_exit(sol_forward.u[end],H,tol_dist) && sol_forward.u[end][5]==0
                prob_backward = ODEProblem(Eq_of_M,u0,(0.,-t_end))
                sol_backward=solve(prob_backward, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_backward,save_start=true,save_end=true,save_everystep=false)

                uf=zeros(5)
                uf[1]=sol_backward.u[end][1] #X
                uf[2]=sol_backward.u[end][2] #P
                uf[3]=sol_backward.u[end][3] #Q
                uf[4]=sol_backward.u[end][4] #Y
                dH=abs(H_test(uf)-H)
                # uf[1]=0
                # prob_test2= ODEProblem(Eq_of_M,u0,(0., -10.))
                # sol_test2=solve(prob_test2, RK4(),maxiters=1e10 )
                # S2=sign(sol_test.u[2][1])
                if  dH<1e-8 #&& S2<0
                    # prob_check= ODEProblem(Eq_of_M,uf,(0., t_end))
                    # sol_check=solve(prob_check, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_forward,save_start=true,save_end=true,save_everystep=false)
                    #return sol_backward.u[end][3],sol_backward.u[end][2],sol_backward.u[end][1],sol_backward.u[end][4],sol_backward.u[end][5],sol_check.u[end][3],sol_check.u[end][2],sol_check.u[end][1],sol_check.u[end][4],Y,sol_check.u[end][5]
                    # return Q,P,sol_backward.u[end][1],Y,sol_backward.u[end][5],sol_check.u[end][3],sol_check.u[end][2],sol_check.u[end][1],sol_check.u[end][4],Yi,sol_check.u[end][5]
                    return sol_backward.u[end][3],sol_backward.u[end][2],sol_backward.u[end][1],sol_backward.u[end][4],sol_backward.u[end][5],Y[1]
                end
            end
        end
    end
end
