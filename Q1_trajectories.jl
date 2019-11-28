using DifferentialEquations
using PolynomialRoots
using Printf
# using MATLAB
# using ProgressMeter


include("leap_frog_definitions.jl")
max_hit_q1 = 6
max_hit_p2 = 10000
barrier = 5
H=.25

condition_max_hits(u, t, integrator) = u[3] >max_hit_q1 || u[6] > max_hit_p2 || maximum([abs(u[1]),abs(u[2]),abs(u[4]),abs(u[5])]) > barrier
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end

function condition_hits_PSS_p2f(u, t, integrator) # Event when event_f(u,t) == 0
    u[5]+100
end

function affect_update_iterator_q1!(integrator)
    q2=integrator.u[2]
    p2=integrator.u[5]
    p1=integrator.u[4]
    P1=P1_poly(q2,p2, H)
    if !isempty(p1) && abs(p1-P1)<1e-5
        integrator.u[3] = integrator.u[3] + 1
    end
end
function affect_update_iterator_p2!(integrator)
    integrator.u[6] = integrator.u[6] + 1
end


callback_max_hits = DiscreteCallback(condition_max_hits, affect_stop!)
callback_hits_PSS_q1_forward = ContinuousCallback(condition_hits_PSS_q1, nothing,  affect_update_iterator_q1!,rootfind=true)
callback_hits_PSS_q1_backward = ContinuousCallback(condition_hits_PSS_q1,  affect_update_iterator_q1!,nothing,rootfind=true)

callback_hits_PSS_p2_forward = ContinuousCallback(condition_hits_PSS_p2f, affect_update_iterator_p2!)
callback_hits_PSS_p2_backward = ContinuousCallback(condition_hits_PSS_p2f, affect_update_iterator_p2!)

# callback_hits_PSS_q1_forward= ContinuousCallback(condition_hits_PSS_q1, affect_update_iterator_q1!)
# callback_hits_PSS_q1_backward= ContinuousCallback(condition_hits_PSS_q1,affect_update_iterator_q1!)

# callback_hits_PSS_p2_forward = ContinuousCallback(condition_hits_PSS_p2,  affect_update_iterator_p2!)
# callback_hits_PSS_p2_backward = ContinuousCallback(condition_hits_PSS_p2, affect_update_iterator_p2!)


cb_forward = CallbackSet(callback_hits_PSS_q1_forward, callback_hits_PSS_p2_forward, callback_max_hits)
cb_backward = CallbackSet(callback_hits_PSS_q1_backward, callback_hits_PSS_p2_backward, callback_max_hits)

function Q1_trajectories(Q2,P2, t_end, H,flag)
    P1 = P1_poly(Q2, P2, H)
    if ~isempty(P1)
        # println(Q1)
        q0, p0 = [zeros(3) for i in 1:2]
        q0[1] = 0
        q0[2] = Q2
        q0[3] = 0
        p0[1] = P1
        p0[2] = P2
        p0[3] = 0
       #constructor for ODE
        prob_f = HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        # sol_f = solve(prob_f,  RK4(),maxiters=1e20,callback = cb_forward, save_start = false, save_end = true, save_everystep = flag)
        sol_f = solve(prob_f, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback = cb_forward, save_start = false, save_end = true, save_everystep = flag)

        # sol_f=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_forward,save_start=false,save_end=true,save_everystep=false)


        prob_b = HamiltonianProblem{true}(Hamiltonian_Dimer_Backwards, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        # sol_b = solve(prob_b,  Vern9(), maxiters = 1e20, reltol = 1e-13, abstol = 1e-15,  callback = cb_backward, save_start = false, save_end = true, save_everystep = flag)
        # sol_b = solve(prob_b,  RK4(), maxiters = 1e20, callback = cb_backward, save_start = false, save_end = true, save_everystep = flag)
        sol_b=solve(prob_b, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_backward,save_start=false,save_end=true,save_everystep=flag)
        point=[sol_f[1,2],sol_f[2,2],sol_f[4,2],sol_f[5,2]]
        return sol_f.t, sol_f[1,:],sol_f[3,end],sol_b.t, sol_b[1,:],sol_b[3,end]

    else
        return -2
    end

end
