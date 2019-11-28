module Dimer_Lagrangian_Descriptor
export Dimer_Lagrangian_Descriptor_Function
using PolynomialRoots
using DifferentialEquations
using Printf
using LinearAlgebra
include("leap_frog_definitions.jl")

barrier=5

hit_max=100000
length_max=30000
condition_escape(u, t, integrator) =  u[5]>length_max || u[6]>hit_max || maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])]) > barrier 

function affect_stop!(integrator) 
    terminate!(integrator)
end

callback_escape = DiscreteCallback(condition_escape, affect_stop!)
function condition_hits_PSS_p2(u, t, integrator) # Event when event_f(u,t) == 0
    u[4] 
end
function affect_update_iterator_p2!(integrator)
    integrator.u[6] = integrator.u[6] + 1
end
callback_hits_PSS_p2 = ContinuousCallback(condition_hits_PSS_p2, affect_update_iterator_p2!, rootfind = false)

cb = CallbackSet(callback_hits_PSS_p2, callback_escape)
function Dimer_Lagrangian_Descriptor_Function(mesh, H,  t_end)
    Q2=mesh[1]
    P2=mesh[2]
     #for a given Q,P,H with X=0
    Point_1=[0,0]
    Point_2=[0,sqrt(2)]
    Point_3=[0,-sqrt(2)]
    Point=[.8*Q2,P2]
    if norm(Point_2-Point)>.8 && norm(Point_3-Point)>.8 && norm(Point_1-Point)>.235
        P1 = P1_poly(Q2, P2, H)
        if ~isempty(P1)
            u0=zeros(6)
            u0[1]=0 #q1
            u0[2]=Q2 #q2
            u0[3]=P1 #p1
            u0[4]=P2 #p2
            u0[5]=0 # lagrange descriptor variable
            u0[6]=0 # counter variable

        #constructor for ODE
  
            prob_b = ODEProblem(Dimer_Eq_of_M_Backwards,u0,(0., t_end))
            prob_f = ODEProblem(Dimer_Eq_of_M,u0,(0., t_end))
            # sol_f=solve(prob_f, Tsit5(),maxiters=1e20)

            # sol_f=solve(prob_f, Tsit5(),maxiters=1e20,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
            # sol_b=solve(prob_b, Tsit5(),maxiters=1e20,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
         
            # sol_f=solve(prob_f, Tsit5(),maxiters=1e20,reltol=1e-4,abstol=1e-7,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
            # sol_b=solve(prob_b, Tsit5(),maxiters=1e20,reltol=1e-4,abstol=1e-7,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
        
            # sol_f=solve(prob_f, Tsit5(),maxiters=1e20,reltol=1e-5,abstol=1e-8,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
            # sol_b=solve(prob_b, Tsit5(),maxiters=1e20,reltol=1e-5,abstol=1e-8,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)

            sol_f=solve(prob_f, Tsit5(),maxiters=1e20,reltol=1e-6,abstol=1e-9,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
            sol_b=solve(prob_b, Tsit5(),maxiters=1e20,reltol=1e-6,abstol=1e-9,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
        
            # sol_f=solve(prob_f, Tsit5(),maxiters=1e20,reltol=1e-7,abstol=1e-10,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
            # sol_b=solve(prob_b, Tsit5(),maxiters=1e20,reltol=1e-7,abstol=1e-10,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)


            # sol_f=solve(prob_f, Tsit5(),maxiters=1e20,reltol=1e-9,abstol=1e-12,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
            # sol_b=solve(prob_b, Tsit5(),maxiters=1e20,reltol=1e-9,abstol=1e-12,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)


            # sol_f=solve(prob_f, Tsit5(),maxiters=1e20,reltol=1e-6,abstol=1e-8,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
            # sol_b=solve(prob_b, Tsit5(),maxiters=1e20,reltol=1e-6,abstol=1e-8,callback=cb,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)

            # sol_f=solve(prob_f,  Vern9(),maxiters=1e20,reltol=1e-14,abstol=1e-14,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
        # sol_b=solve(prob_b,  Vern9(),maxiters=1e20,reltol=1e-14,abstol=1e-14,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol
        # sol_b=solve(prob_b, RK4(),maxiters=1e20,reltol=1e-6,abstol=1e-8)
        
        
        # sol_f=solve(prob_f, RK4(),maxiters=1e20,reltol=1e-8,abstol=1e-10)
        # sol_b=solve(prob_b, RK4(),maxiters=1e20,reltol=1e-8,abstol=1e-10)
        
        # sol_b=solve(prob_b, RK4(),maxiters=1e20)
        # sol_f=solve(prob_f, RK4(),maxiters=1e20)

        # sol_b=solve(prob_b, RK4(),maxiters=1e20,callback=callback_escape)
        # sol_f=solve(prob_f, RK4(),maxiters=1e20,callback=callback_escape)
      
        # sol_f=solve(prob_f, Vern9(),maxiters=1e20,reltol=1e-13,abstol=1e-15)
        # sol_b=solve(prob_b, Vern9(),maxiters=1e20,reltol=1e-13,abstol=1e-15)


        #find final value to check if energy was conserved
        # uf=zeros(5)
        # uf[1]=sol[1,end] #X
        # uf[2]=sol[2,end] #P
        # uf[3]=sol[3,end] #Q
        # uf[4]=sol[4,end] #Y
        # dH=abs(H_test(uf)-H)
        LD=sol_f[1,end]+sol_b[1,end]
           return LD
        else
            return NaN
        end
    else
        return NaN
    end
end

end