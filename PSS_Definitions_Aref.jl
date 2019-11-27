include("leap_frog_definitions.jl")
max_hit=500
barrier=10

condition_max_hits(u,t,integrator)= u[5]>max_hit || maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])])>barrier
affect_stop!(integrator) = terminate!(integrator)
function condition_hits_PSS(u,t,integrator) # Event when event_f(u,t) == 0
   u[1]
end
function affect_update_iterator!(integrator)
    integrator.u[5]=integrator.u[5]+1
end


callback_max_hits=DiscreteCallback(condition_max_hits,affect_stop!)
callback_hits_PSS=ContinuousCallback(condition_hits_PSS, affect_update_iterator!,nothing)
cb=CallbackSet(callback_hits_PSS, callback_max_hits)

function PSS_function(Q,P, H,  t_end)
    #setting tolerance for determing if vortices move off to infinity
    tol_dist=1e-3

    #for a given Q,P,H with X=0
    Y=Yfind_Aref(Q,P,H)
    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y #Y
        u0[5]=0 # counter variable

        #constructor for ODE
        prob = ODEProblem(Eq_of_M,u0,(0., t_end))

        #solve ode , save_everystep=false is important to prevent sol to include all points, not just event points
        # sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb,save_start=true,save_end=true,save_everystep=false)
        sol=solve(prob, RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)

        #find final value to check if energy was conserved
        uf=zeros(5)
        uf[1]=sol[1,end] #X
        uf[2]=sol[2,end] #P
        uf[3]=sol[3,end] #Q
        uf[4]=sol[4,end] #Y
        dH=abs(H_test(uf)-H)

        #output Q, P and dH
        return sol[:,2:end-1][3,:],sol[:,2:end-1][2,:], dH
    else
    #need to return 3 values, dH=1 flags that there is no Y
        return 0,0,1
    end
end
