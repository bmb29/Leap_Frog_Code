include("leap_frog_definitions.jl")
max_hit=10000
barrier=20


condition_max_hits(u,t,integrator)= u[6]>max_hit
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS(u,t,integrator) # Event when event_f(u,t) == 0
   u[4]
end
function affect_update_iterator!(integrator)
    integrator.u[6]=integrator.u[6]+1
end


callback_max_hits=DiscreteCallback(condition_max_hits,affect_stop!)
callback_hits_PSS=ContinuousCallback(condition_hits_PSS, affect_update_iterator!,nothing)
cb=CallbackSet(callback_hits_PSS, callback_max_hits)

function PSS_function(Q2,P2, H,  t_end)
    #for a given Q,P,H with X=0
    Q1=Q1_find_dimer(Q2,P2,H)
    if ~isempty(Q1)
        # println(Q1)
        q0,p0=[zeros(3) for i in 1:2]
        q0[1]=Q1[1]
        q0[2]=Q2
        q0[3]=0
        p0[1]=0
        p0[2]=P2
        p0[3]=0
       #constructor for ODE
        prob= HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol=solve(prob, RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-6,callback=cb,save_start=true,save_end=true,save_everystep=false)
        #output Q, P and dH
        return sol[:,2:end-1][2,:],sol[:,2:end-1][5,:]
    else
    #need to return 3 values, dH=1 flags that there is no Y
        return 0, 0
    end
end
