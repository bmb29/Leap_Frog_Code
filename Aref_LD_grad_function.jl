module Aref_LD_grad_function
export gradM
export final_T
export LG_Helper


using Roots
using DifferentialEquations
using Printf
include("leap_frog_definitions.jl")


barrier=5
condition_escape(u, t, integrator) =  maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])]) > barrier
affect_stop!(integrator) = terminate!(integrator)
callback_escape = DiscreteCallback(condition_escape, affect_stop!)
function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end
cb=ContinuousCallback(condition_hits_PSS_q1, affect_stop!, rootfind=true)


function gradM(mesh, H, t_end)
    delta=1e-7
    Q=mesh[2]
    P=mesh[1]

    mXb=[P-delta,Q]
    mXf=[P+delta,Q]
    mYb=[P,Q-delta]
    mYf=[P,Q+delta]

    ldXb=LD_Helper(mXb,H,t_end)
    ldXf=LD_Helper(mXf,H,t_end)
    # ldYb=LD_Helper(mYb,H,t_end)
    # ldYf=LD_Helper(mYf,H,t_end)
    dLDx=(ldXf-ldXb)/(2*delta)
    # dLDy=(ldYf-ldYb)/(2*delta)
    # M=(dLDx)^2+(dLDy)^2
    return maximum([log10(abs(dLDx)),4])

end

function LD_Helper(mesh, H,  t_end)
    Q=mesh[2]
    P=mesh[1]
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
  
        prob_b = ODEProblem(Eq_of_M_Lagrangian_Descriptors_Backwards,u0,(0., t_end))
        prob_f = ODEProblem(Eq_of_M_Lagrangian_Descriptors,u0,(0., t_end))
        

        sol_f=solve(prob_f, Tsit5(),maxiters=1e20,reltol=1e-10,abstol=1e-13,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
        # sol_b=solve(prob_b, Tsit5(),maxiters=1e20,reltol=1e-9,abstol=1e-12,save_idxs = [5],save_every_step=false, save_end=true, dense=false)#,abstol=1e-9)
        
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

        LD=sol_f[1,end]
        # LD=sol_b[1,end]
        # if maximum([abs(sol_f[1,end]),abs(sol_f[2,end]),abs(sol_f[3,)nd]),abs(sol_f[4,end]),abs(sol_b[1,end]),abs(sol_b[2,end]),abs(sol_b[3,end]),abs(sol_b[4,end])])<barrier
        return LD
    else
        return NaN
    end
end
function gradient_matrix(LD_Matrix,dx,dy)
    N,N=size(LD_Matrix)
    grad_X_Matrix=zeros(N-4,N-4)
    grad_Y_Matrix=zeros(N-4,N-4)
    for i=1:N-4
        for j=1:N-4
            grad_X_Matrix[i,j]=-LD_Matrix[i+4,j+1]+8*LD_Matrix[i+3,j+1]-8*LD_Matrix[i+1,j+1]+LD_Matrix[i,j+1]
            grad_Y_Matrix[i,j]=-LD_Matrix[i+1,j+4]+8*LD_Matrix[i+1,j+3]-8*LD_Matrix[i+1,j+1]+LD_Matrix[i+1,j]
        end
    end
    X= (grad_X_Matrix/(12*dx) ).^2
    Y= (grad_Y_Matrix/(12*dy) ).^2
    return X,Y
end

function final_T(H)
    Q=1e-8
    P=1e-8
    Y=Yfind_Aref(Q,P,H)
    if ~isempty(Y)
        u0=zeros(5)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y #Y
        u0[5]=0 # counter variable

        #constructor for ODE
  
        prob = ODEProblem(Eq_of_M_Lagrangian_Descriptors_Backwards,u0,(0., 20))
        sol=solve(prob, Vern9(),maxiters=1e20,reltol=1e-13,abstol=1e-15,callback=cb)
     
        return sol.t[end]
    end
end

function gradient_matrix_2(LD_Matrix,dx,dy)
    N,N=size(LD_Matrix)
    grad_X_Matrix=zeros(N-2,N-2)
    grad_Y_Matrix=zeros(N-2,N-2)
    for i=1:N-2
        for j=1:N-2
            grad_X_Matrix[i,j]=LD_Matrix[i+2,j+1]-LD_Matrix[i,j+1]
            grad_Y_Matrix[i,j]=LD_Matrix[i+1,j+2]-LD_Matrix[i+1,j]
        end
    end
    grad_M=( grad_X_Matrix/(2*dx) ).^2+( grad_Y_Matrix/(2*dy) ).^2
    return sqrt.(grad_M)
end

function gradient_matrix_4(LD_Matrix,dx,dy)
    N,N=size(LD_Matrix)
    grad_X_Matrix=zeros(N-4,N-4)
    grad_Y_Matrix=zeros(N-4,N-4)
    for i=1:N-4
        for j=1:N-4
            grad_X_Matrix[i,j]=-LD_Matrix[i+4,j+1]+8*LD_Matrix[i+3,j+1]-8*LD_Matrix[i+1,j+1]+LD_Matrix[i,j+1]
            grad_Y_Matrix[i,j]=-LD_Matrix[i+1,j+4]+8*LD_Matrix[i+1,j+3]-8*LD_Matrix[i+1,j+1]+LD_Matrix[i+1,j]
        end
    end
    grad_M=( grad_X_Matrix/(12*dx) ).^2+( grad_Y_Matrix/(12*dy) ).^2
    return sqrt.(grad_M)
end



function D_X(LD_Matrix,dx)
    N=size(LD_Matrix)
    N=N-6
    grad_X_Matrix=zeros(N,N)
    for i=1:N
        for j=1:N
            grad_X_Matrix[i,j]=-LD_Matrix[i+4,j+1]+8*LD_Matrix[i+3,j+1]-8*LD_Matrix[i+1,j+1]+LD_Matrix[i,j+1]
        end
    end
    grad_M=abs(grad_X_Matrix/(12*dx) )
    return grad_M
end


end