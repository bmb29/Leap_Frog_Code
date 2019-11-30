include("leap_frog_definitions.jl")
using PolynomialRoots
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
callback_hits_PSS_q1 = ContinuousCallback(condition_hits_PSS_q1, nothing, affect_update_iterator_q1!, rootfind=true)
callback_hits_PSS_p2 = ContinuousCallback(condition_hits_PSS_p2,  affect_update_iterator_p2!, save_positions = (false,false))

cb = CallbackSet(callback_hits_PSS_q1, callback_hits_PSS_p2, callback_max_hits)
# cb = CallbackSet(callback_hits_PSS_q1, callback_max_hits)


function PSS_function(Q2,P2, H,  t_end)
    #for a given Q,P,H with X=0
    P1 = P1_poly(Q2, P2, H)
    if ~isempty(P1)
        # println(Q1)
        q0,p0=[zeros(3) for i in 1:2]
        q0[1]=0
        q0[2]=Q2
        q0[3]=0
        p0[1]=P1[1]
        p0[2]=P2
        p0[3]=0
       #constructor for ODE
        prob= HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        sol=solve(prob, Tsit5(), maxiters=1e20, reltol=1e-8, abstol=1e-11, callback=cb, save_start=true, save_end=true, save_everystep=false)
        # sol=solve(prob, RK4(),maxiters=1e20, reltol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)
        #output Q, P and dH
        Q1=sol[:,2:end-1][1,:]
        bool_filter=[abs(q1)<1e-6 for q1 in Q1]
        return sol[:,2:end-1][2,:][bool_filter],sol[:,2:end-1][5,:][bool_filter]
        # return sol[:,2:end-1][2,:],sol[:,2:end-1][5,:]

    else
    #need to return 3 values, dH=1 flags that there is no Y
        return 0, 0
    end
end
