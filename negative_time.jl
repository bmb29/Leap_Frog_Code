
# module escape_dimer_one_hit
# export escape_exit_num_dimer
using DifferentialEquations
using Roots
using Printf
# using MATLAB
# using ProgressMeter


include("leap_frog_definitions.jl")
max_hit_q1_back = 1
barrier = 5


condition_max_hits_back_one(u, t, integrator) = u[3] == max_hit_q1_back
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end

function affect_update_iterator_q1!(integrator)
    integrator.u[3] = integrator.u[3] + 1
end



callback_max_hits_back_one = DiscreteCallback(condition_max_hits_back_one, affect_stop!)
callback_hits_PSS_q1 = ContinuousCallback(condition_hits_PSS_q1,nothing,affect_update_iterator_q1!)
cb_neg = CallbackSet(callback_hits_PSS_q1, callback_max_hits_back_one)

function back_one(Q2,P2, t_end, H)
    P1 = P1_poly(Q2, P2, H)
    if ~isempty(P1) 
        q0, p0 = [zeros(3) for i in 1:2]
        q0[1] = 0
        q0[2] = Q2
        q0[3] = 0
        p0[1] = P1[1]
        p0[2] = P2
        p0[3] = 0
       #constructor for ODE
        prob = HamiltonianProblem{true}(Hamiltonian_Dimer_Backwards, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        # sol=solve(prob, RK4(),maxiters=1e20,callback=cb,save_start=false,save_end=true,save_everystep=false)

        # sol=solve(prob, RK4(),maxiters=1e20,reltol=1e-8,abstol=1e-10,callback=cb,save_start=true,save_end=true,save_everystep=false)
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb_neg,save_start=false,save_end=true,save_everystep=false)

        uf=zeros(6)
        uf[1]=sol[1,end]
        uf[2]=sol[2,end]
        uf[3]=sol[3,end]
        uf[4]=sol[4,end]
        uf[5]=sol[5,end]
        uf[6]=sol[6,end]
        q1=uf[1]
        q2=uf[2]
        p1=uf[4]
        p2=uf[5]
        q=[q1, q2]
        p=[p1, p2]
        diff=abs(Hamiltonian_Dimer(q,p,0)-H)


        if uf[3]==1 && diff<1e-10 
            return q1,q2,p1,p2 
       
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