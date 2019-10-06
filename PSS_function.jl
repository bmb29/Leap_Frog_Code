using Roots
using DifferentialEquations
using Printf
include("leap_frog_definitions.jl")
#equations used for describing motion and conserved quantities

#condition for when to stop integrating if we hit the PSS enough times
condition1(u,t,integrator)= u[5]>max_hit
#affect used when counter hits max
affect1!(integrator) = terminate!(integrator)

#condition for when we hit PSS, X=0
function condition2(u,t,integrator)
u[1]
end

#affect when we hit X=0, update counter
function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
end

#construct callback for when counter hits and terminates integration
cb1 = DiscreteCallback(condition1,affect1!)

#callback for when we hit PSS

# affect2, nothing  ===> upcrossing    neg->pos
# nothing, affect2  ===> downcrossing  pos ->neg
# affect2           ===> event happens at both
#save_postion in list of points only saves after the even rather than before and/or after
cb2 = ContinuousCallback(condition2,affect2!,nothing)

#creates callback set, put ContinuousCallback first since the update will affect DiscreteCallback
cb=CallbackSet(cb2,cb1)



function PSS_function(Q,P, Energy,  t_end,    max_hit)
    #setting tolerance for determing if vortices move off to infinity
    tol_dist=1e-3

    #converting energy into the hamiltonian used to find the Eq_of_M
    H=(2*Energy)^2

    #for a given Q,P,H with X=0
    Y=Yfind(Q,P,H)
    
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
        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb,save_start=true,save_end=true,save_everystep=false)
        # sol=solve(prob, RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)

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
