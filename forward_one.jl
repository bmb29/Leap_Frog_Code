
# module escape_dimer_one_hit
# export escape_exit_num_dimer
using DifferentialEquations
using Roots
using PolynomialRoots

using Printf
# using MATLAB
# using ProgressMeter


include("leap_frog_definitions.jl")
max_hit_q1_forward_one = 1
H=.25

condition_max_hits_forward_one(u, t, integrator) = u[3] == max_hit_q1_forward_one
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end
function affect_update_iterator_q1!(integrator)
    q2=integrator.u[2]
    p2=integrator.u[5]
    p1=integrator.u[4]
    P1=P1_poly(q2,p2, H)
    if !isempty(P1) && abs(p1-P1)<1e-5
        integrator.u[3] = integrator.u[3] + 1
    end
end


callback_max_hits_forward_one = DiscreteCallback(condition_max_hits_forward_one, affect_stop!)
callback_hits_PSS_q1_forward_one = ContinuousCallback(condition_hits_PSS_q1, nothing,  affect_update_iterator_q1!, rootfind = true)
cb_forward_one = CallbackSet(callback_hits_PSS_q1_forward_one, callback_max_hits_forward_one)
function forward_one(Q2, P2, t_end, H)


    P1 = P1_poly(Q2, P2, H)
    if ~isempty(P1) 
        q0, p0 = [zeros(3) for i in 1:2]
        q0[1] = 0
        q0[2] = Q2
        q0[3] = 0
        p0[1] = P1
        p0[2] = P2
        p0[3] = 0
       #constructor for ODE
        prob = HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        # sol = solve(prob,  Vern9(), maxiters = 1e20, reltol = 1e-13, abstol = 1e-15, callback = cb_forward_one, save_start = false, save_end = true, save_everystep = false)
        # sol=solve(prob,  RK4(), maxiters = 1e20, callback=cb_forward_one,save_start=false,save_end=true,save_everystep=false)
        # sol=solve(prob, RK4(),maxiters=1e20,reltol=1e-8,abstol=1e-10,callback=cb,save_start=true,save_end=true,save_everystep=false)
        sol = solve(prob, Vern9(), maxiters = 1e20, reltol = 1e-13, abstol = 1e-15, callback = cb_forward_one, save_start = true, save_end = false, save_everystep = false)
        # uf = zeros(6)
        # uf[1] = sol[1,3]
        # uf[2] = sol[2,3]
        # uf[3] = sol[3,3]
        # uf[4] = sol[4,3]
        # uf[5] = sol[5,3]
        # uf[6] = sol[6,3]
        # q1 = uf[1]
        # q2 = uf[2]
        # p1 = uf[4]
        # p2 = uf[5]
        if  sol[3,end] == 1 
            uf = zeros(6)
            uf[1] = sol[1,3]
            uf[2] = sol[2,3]
            uf[3] = sol[3,3]
            uf[4] = sol[4,3]
            uf[5] = sol[5,3]
            uf[6] = sol[6,3]
            q1 = uf[1]
            q2 = uf[2]
            p1 = uf[4]
            p2 = uf[5]
            q = [q1, q2]
            p = [p1, p2]
            diff = abs(Hamiltonian_Dimer(q, p, 0) - H)
            if diff < 1e-1  
                return q1, q2, p1, p2
            end
        end
        # if dH<1e-5
        #     if is_it(uf,H,tol_dist)
        #         return sol.u[end][5]
        #     else
        #         return -1
        #     end
        # end
    end
end
function unzipper2(array_of_points)
    N = length(array_of_points)
    Q = zeros(N)
    P = zeros(N)
    for (index, point) in enumerate(array_of_points)
        Q[index] = point[2]
        P[index] = point[4]
    end
    return Q, P
end
function unzipper3(array_of_points)
    N = length(array_of_points)
    Q0 = zeros(N)
    P0 = zeros(N)
    Q1 = zeros(N)
    P1 = zeros(N)
    for (index, point) in enumerate(array_of_points)
        Q0[index] = point[2]
        P0[index] = point[4]
        Q1[index] = point[6]
        P1[index] = point[8]
    end
    return Q0, P0, Q1, P1
end

function unzipper3(array_of_points)
    N = length(array_of_points)
    Q0 = zeros(N)
    P0 = zeros(N)
    Q1 = zeros(N)
    P1 = zeros(N)
    for (index, point) in enumerate(array_of_points)
        Q0[index] = point[2]
        P0[index] = point[4]
        Q1[index] = point[6]
        P1[index] = point[8]
    end
    return Q0, P0, Q1, P1
end
function unzipper_quad(array_of_points)
    N = length(array_of_points)
    Q = zeros(N)
    P = zeros(N)
    for (index, point) in enumerate(array_of_points)
        Q[index] = point[2]
        P[index] = point[4]
    end
    full_Q = vcat(Q, Q, -Q, -Q)
    full_P = vcat(P, -P, P, -P)
    return full_Q, full_P
end
function unzipper_half(array_of_points)
    N = length(array_of_points)
    Q = zeros(N)
    P = zeros(N)
    for (index, point) in enumerate(array_of_points)
        Q[index] = point[2]
        P[index] = point[4]
    end
    full_Q = vcat(Q, -Q)
    full_P = vcat(P, -P)
    return full_Q, full_P
end