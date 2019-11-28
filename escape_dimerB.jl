
module escape_dimerB
export escape_exit_num_dimer
using DifferentialEquations
using PolynomialRoots
using Printf
# using MATLAB
# using ProgressMeter

include("Dimer_Eq_of_M.jl")
include("leap_frog_definitions.jl")
H=.25
max_hit_q1 = 11
max_hit_p2 = 1000
barrier = 5
condition_max_hits_back_forward(u, t, integrator) = u[5] > max_hit_q1 || u[6] > max_hit_p2 || maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])]) > barrier
affect_stop!(integrator) = terminate!(integrator)

function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end
function condition_hits_PSS_p2(u, t, integrator) # Event when event_f(u,t) == 0
    u[4]
end


# function affect_update_iterator_q1!(integrator)
#     integrator.u[5] = integrator.u[5] + 1
# end

function affect_update_iterator_q1!(integrator)
    # integrator.u[5] = integrator.u[5] + 1
    q2 = integrator.u[2]
    p1 = integrator.u[3]
    p2 = integrator.u[4]
    P1 = P1_poly(q2, p2, H)
    if !isempty(P1) &&  abs(p1 - P1) < .5
        integrator.u[5] = integrator.u[5] + 1
    end
end

function affect_update_iterator_p2!(integrator)
    integrator.u[6] = integrator.u[6] + 1
end
callback_max_hits_bf = DiscreteCallback(condition_max_hits_back_forward, affect_stop!)
callback_hits_PSS_q1_forward = ContinuousCallback(condition_hits_PSS_q1, nothing,  affect_update_iterator_q1!, rootfind = true)
callback_hits_PSS_q1_backward = ContinuousCallback(condition_hits_PSS_q1,  affect_update_iterator_q1!, nothing, rootfind = true)


callback_hits_PSS_p2_forward = ContinuousCallback(condition_hits_PSS_p2, affect_update_iterator_p2!, rootfind = false)
callback_hits_PSS_p2_backward = ContinuousCallback(condition_hits_PSS_p2, affect_update_iterator_p2!, rootfind = false)


cb_forward = CallbackSet(callback_hits_PSS_q1_forward, callback_hits_PSS_p2_forward, callback_max_hits_bf)
cb_backward = CallbackSet(callback_hits_PSS_q1_backward, callback_hits_PSS_p2_backward, callback_max_hits_bf)

function escape_exit_num_dimer(mesh_list, t_end, H)
    Q2 = mesh_list[1]; P2 = mesh_list[2];
    P1 = P1_poly(Q2, P2, H)
    if ~isempty(P1)
        u0=zeros(6)
        u0[1]=0
        u0[2]=Q2
        u0[3]=P1
        u0[4]=P2
        u0[5]=0
        u0[6]=0

       #constructor for ODE
        # prob = HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        prob=ODEProblem(Dimer_Eq_of_M_hits,u0,(0., t_end) )
        sol_f = solve(prob, RK4(), maxiters = 1e20, reltol = 1e-8, abstol = 1e-10, callback = cb_forward, save_start = false, save_end = true, save_everystep = true)

        q1 = sol_f[1,end]
        q2 = sol_f[2,end]
        p1 = sol_f[3,end]
        p2 = sol_f[4,end]

        q = [q1, q2]
        p = [p1, p2]

        forward_hits = sol_f[5,end]

        # println(sol_f.t)

        diff = abs(Hamiltonian_Dimer(q, p, 0) - H)
        if length(sol_f[1,:]) > 1
            sign_sol = sign(sol_f[1,2])
        else
            sign_sol = 0
        end
        if maximum([abs(q1),abs(q2),abs(p1),abs(p2)]) > barrier
            return forward_hits
        else
            return -1
        end 
            # # return forward_hits
        #     q0, p0 = [zeros(3) for i in 1:2]
        #     q0[1] = 0
        #     q0[2] = Q2
        #     q0[3] = 0
        #     p0[1] = P1
        #     p0[2] = P2
        #     p0[3] = 0
        #     prob_b = HamiltonianProblem{true}(Hamiltonian_Dimer_Backwards, q0, p0, (0., t_end));
        #     #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        #     # sol_b=solve(prob_b,  RK4(), maxiters = 1e20,callback=cb_backward,save_start=false,save_end=true,save_everystep=false)
        #     # sol_b=solve(prob_b,  RK4(), maxiters = 1e20,reltol=1e-8,abstol=1e-10, callback=cb_backward,save_start=false,save_end=true,save_everystep=false)
        #     sol_b=solve(prob_b, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_backward,save_start=false,save_end=true,save_everystep=false)
    
        #     uf = zeros(6)
        #     uf[1] = sol_b[1,end]
        #     uf[2] = sol_b[2,end]
        #     uf[3] = sol_b[3,end]
        #     uf[4] = sol_b[4,end]
        #     uf[5] = sol_b[5,end]
        #     uf[6] = sol_b[6,end]
        #     q1 = uf[1]
        #     q2 = uf[2]
        #     p1 = uf[4]
        #     p2 = uf[5]
        #     q = [q1, q2]
        #     p = [p1, p2]
            
        #     backward_hits = uf[3]
        #     q = [q1, q2]
        #     p = [p1, p2]
        #     diff = abs(Hamiltonian_Dimer(q, p, 0) - H)
        #     if diff < 1e-10 && maximum([abs(q1),abs(q2),abs(p1),abs(p2)]) > barrier 
        #         return forward_hits+ backward_hits# sol_f.t, sol_b.t
        #     else
        #         return -1
        #     end  
        # else
        #     return -1  
        # end
        # if dH<1e-5
        #     if is_it(uf,H,tol_dist)
        #         return sol.u[end][5]
        #     else
        #         return -1
        #     end
        # end
    else
        return -1
    end

end
end
