
module escape
export escape_exit_function_parallel
using DifferentialEquations
using Roots
using ProgressMeter
using Printf
using MATLAB
include("is_it.jl")
include("leap_frog_definitions.jl")


condition(u,t,integrator)= u[5]>5
function condition2(u,t,integrator) # Event when event_f(u,t) == 0
    u[1]
end
affect!(integrator) = terminate!(integrator)

function affect2!(integrator)
    integrator.u[5]=integrator.u[5]+1
end

cb2=ContinuousCallback(condition2,affect2!,nothing)
cb1=DiscreteCallback(condition,affect!)
cb=CallbackSet(cb2,cb1)


function escape_exit_function_parallel(mesh_list,t_end,Energy_A)
    Q=mesh_list[1]; P=mesh_list[2];
    H=(2*Energy_A)^2
    tol_dist=1e-5
    Y=Yfind(Q,P,H)
    print(Y)
    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y #Y
        u0[5]=0 #Y

        prob = ODEProblem(Eq_of_M,u0,(0., t_end))
        # sol=solve(prob,RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb)
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb)
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


end
