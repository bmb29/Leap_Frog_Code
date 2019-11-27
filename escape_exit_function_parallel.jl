#
include("Yfind.jl")
include("is_it_exit2.jl")
include("is_it.jl")


@everywhere function escape_exit_function_parallel(mesh_list,t_end,Energy)
    Q=mesh_list[1]; P=mesh_list[2];
    H=(2*Energy)^2
    tol_dist=1e-5
    Y=Yfind(Q,P,H)
    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y[1] #Y
        u0[5]=0 #Y

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
        if dH<1e-5
            if is_it(uf,H,tol_dist)
                return sol.u[end][5]
            else
                return -1
            end
        end
    end
return -2
end
