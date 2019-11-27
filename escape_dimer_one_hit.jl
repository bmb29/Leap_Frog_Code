
# module escape_dimer_one_hit
# export escape_exit_num_dimer
using DifferentialEquations
using Roots
using Printf
# using MATLAB
# using ProgressMeter


include("leap_frog_definitions.jl")
max_hit_q1 = 1
barrier = 10


condition_max_hits(u, t, integrator) = u[3] > max_hit_q1
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end


function affect_update_iterator_q1!(integrator)
    integrator.u[3] = integrator.u[3] + 1
end


callback_max_hits = DiscreteCallback(condition_max_hits, affect_stop!)
callback_hits_PSS_q1 = ContinuousCallback(condition_hits_PSS_q1, nothing, affect_update_iterator_q1!, rootfind = true)

cb = CallbackSet(callback_hits_PSS_q1, callback_max_hits)

function escape_after_one_hit(point, t_end, H)
    Q2 = point[1]; P2 = point[2];
    P1 = P1_find_dimer(Q2, P2, H)
    Q1 = 0

    # if isempty(P1)
    #     P1 = P1_find_dimer_second(Q2, P2, H)
    #     t_end = 1
    # end

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
        sol=solve(prob, RK4(),maxiters=1e20,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol=solve(prob, RK4(),maxiters=1e20,reltol=1e-8,abstol=1e-10,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol = solve(prob, Vern9(), maxiters = 1e20, reltol = 1e-13, abstol = 1e-15, callback = cb, save_start = true, save_end = true, save_everystep = false)

        uf = zeros(6)
        uf[1] = sol[1,end]
        uf[2] = sol[2,end]
        uf[3] = sol[3,end]
        uf[4] = sol[4,end]
        uf[5] = sol[5,end]
        uf[6] = sol[6,end]
        if maximum([abs(uf[1]),abs(uf[2]),abs(uf[4]),abs(uf[5])]) > barrier && uf[3] == 1
            q1 = sol[1, 2]
            q2 = sol[2, 2]
            p1 = sol[4, 2]
            p2 = sol[5, 2]
            # q1 = sol[1, end-1]
            # q2 = sol[2, end-1]
            # p1 = sol[4, end-1]
            # p2 = sol[5, end-1]
            q = [q1, q2]
            p = [p1, p2]
            diff = abs(Hamiltonian_Dimer(q, p, 0) - H)

            # prob_sign = HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., .001));
            # sol_sign = solve(prob_sign, RK4(), maxiters = 1e20)
            # sign_orientation = sign(sol_sign[1,3])

            if diff < 1e-2
                # return Q1, Q2, P1[1], P2, q1, q2, p1, p2
                return Q1, Q2, P1[1], P2, q1, q2, p1, p2

                # return sol
            end
        end
        # return sol.t, sol[1,:], sol[5,:]
        # else

        # end
        # if uf[3]==1
        #     return uf[1], uf[2], uf[4], uf[5]
        # end


        # if dH<1e-5
        #     if is_it(uf,H,tol_dist)
        #         return sol.u[end][5]
        #     else
        #         return -1
        #     end
        # end1233
    # else
        # return NaN, NaN, NaN, NaN
    end
end

function add_four!(point, radius, one_hit_exit, t_end, energy)
    X = point[2]
    Y = point[4]
    X_plus = X + radius
    Y_plus = Y + radius
    X_minus = X - radius
    Y_minus = Y - radius
    
    
    point_2 = [X, Y_plus]
    
    point_4 = [X_minus, Y]
    
    point_5 = [X_plus, Y]
    
    point_7 = [X, Y_minus]
    
    
    new_points = [point_2, point_4, point_5, point_7]
    for point in new_points
        new_point = escape_after_one_hit(point, t_end, energy)
        if new_point != nothing
            push!(one_hit_exit, new_point)
        end
    end
end

function add_eight!(point, radius, one_hit_exit, t_end, energy)
    X = point[2]
    Y = point[4]
    
    X_plus = X + radius
    Y_plus = Y + radius
    X_minus = X - radius
    Y_minus = Y - radius
    
    
    point_1 = [X_minus, Y_plus]
    point_2 = [X, Y_plus]
    point_3 = [X_minus, Y_plus]
    
    point_4 = [X_minus, Y]
    point_5 = [X_plus, Y]
    
    point_6 = [X_minus, Y_minus]
    point_7 = [X, Y_minus]
    point_8 = [X_minus, Y_minus];
    
    
    new_points = [point_1, point_2, point_3, point_4, point_5, point_6, point_7, point_8]
    for point in new_points
        new_point = escape_after_one_hit(point, t_end, energy)
        if new_point != nothing
            push!(one_hit_exit, new_point)
        end
    end
end

function unzipper(array_of_points)
    N = length(array_of_points)
    Q = zeros(N)
    P = zeros(N)
    for (index, point) in enumerate(array_of_points)
        Q[index] = point[1]
        P[index] = point[2]
    end
    return Q, P
end

function unzipper3(array_of_points)
    N=length(array_of_points)
    Q0=zeros(N)
    P0=zeros(N)
    Q1=zeros(N)
    P1=zeros(N)
    for (index,point) in enumerate(array_of_points)
        Q0[index]=point[2]
        P0[index]=point[4]
        Q1[index]=point[6]
        P1[index]=point[8]
    end
    return Q0,P0,Q1,P1
end