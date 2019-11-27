
module escape_num_dimer
export escape_exit_num_dimer
using DifferentialEquations
using Roots
using Printf

# using MATLAB
# using ProgressMeter


include("leap_frog_definitions.jl")
max_hit_q1 = 10000
max_hit_p2 = 10000
barrier = 5


condition_max_hits(u, t, integrator) = u[3] > max_hit_q1 || u[6]>max_hit_p2 || maximum([abs(u[1]),abs(u[2]),abs(u[4]),abs(u[5])]) > barrier
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end

function condition_hits_PSS_p2(u, t, integrator) # Event when event_f(u,t) == 0
    u[5]
end

function affect_update_iterator_q1!(integrator)
    integrator.u[3] = integrator.u[3] + 1
end
function affect_update_iterator_p2!(integrator)
    integrator.u[6] = integrator.u[6] + 1
end


callback_max_hits = DiscreteCallback(condition_max_hits, affect_stop!)
callback_hits_PSS_q1 = ContinuousCallback(condition_hits_PSS_q1, affect_update_iterator_q1!, nothing)
callback_hits_PSS_p2 = ContinuousCallback(condition_hits_PSS_p2, affect_update_iterator_p2!, nothing)

cb = CallbackSet(callback_hits_PSS_q1, callback_hits_PSS_p2, callback_max_hits)

function escape_exit_num_dimer(mesh_list, t_end, H)
    Q2 = mesh_list[1]; P2 = mesh_list[2];
    P1 = P1_find_dimer(Q2, P2, H)

    if isempty(P1)
        P1 = P1_find_dimer_second(Q2, P2, H)
    end

    if ~isempty(P1)
        # println(Q1)
        q0, p0 = [zeros(3) for i in 1:2]
        q0[1] = 0
        q0[2] = Q2
        q0[3] = 0
        p0[1] = P1[1]
        p0[2] = P2
        p0[3] = 0
       #constructor for ODE
        prob = HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        sol = solve(prob, Vern9(),  reltol=1e-13,abstol=1e-15,maxiters = 1e10,  callback = cb, save_start = true, save_end = true, save_everystep = false)
        # uf=zeros(6)
        # uf[1]=sol[1,end]
        # uf[2]=sol[2,end]
        # uf[3]=sol[3,end]
        # uf[4]=sol[4,end]
        # uf[5]=sol[5,end]
        # uf[6]=sol[6,end]

        if sol.u[end][6] > max_hit_p2 - 1 ||  sol.u[end][3] > max_hit_q1 - 1
            return t_end
        else
            return sol.t[end]

        end
        # if dH<1e-5
        #     if is_it(uf,H,tol_dist)
        #         return sol.u[end][5]
        #     else
        #         return -1
        #     end
        # end
    else
        return NaN
    end
end


end
