using DifferentialEquations
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

include("Dimer_Eq_of_M.jl")
Ham2ndSHO(p, q, t) = 8 * (p[1]^2 + p[2]^2) + (q[1]^2 + q[2]^2);
ω_0 = [1.,1e4]
tspan = (0., 100);
p0 = [1.0,0]
q0 = [0,1.0]
params_1 = [1,2]
params_2 = [1,sqrt(2)]
param = sqrt(2)
# prob = HamiltonianProblem{true}(Ham2ndSHO, p0, q0, tspan);
# sol = solve(prob, RK4(), maxiters = 1e20, reltol = 1e-8, abstol = 1e-10)

# prob_1 = HamiltonianProblem{true}(Ham2ndSHO, p0, q0, params_1, tspan);
# sol_1 = solve(prob_1, RK4(), maxiters = 1e20, reltol = 1e-8, abstol = 1e-10)

# prob_2 = HamiltonianProblem{true}(Ham2ndSHO, p0, q0, tspan);
# sol_2 = solve(prob_2, RK4(), maxiters = 1e20, reltol = 1e-8, abstol = 1e-10)


function condition_hits_PSS_q1(u, t, integrator) # Event when event_f(u,t) == 0
    u[1]
end


function affect_update_iterator_q1!(integrator)
    integrator.u[5] = integrator.u[5] + 1
end


u0=[0.,0,1,0]
callback_hits_PSS_q1_forward = ContinuousCallback(condition_hits_PSS_q1, nothing,  affect_update_iterator_q1!, rootfind = true)
cb_forward = CallbackSet(callback_hits_PSS_q1_forward)
prob=ODEProblem(Dimer_Eq_of_M,u0,tspan )
sol = solve(prob, RK4(), maxiters = 1e20, reltol = 1e-8, abstol = 1e-10, callback = cb_forward, save_start = false, save_end = true, save_everystep = false)
plot(sol[2,:],sol[4,:],".")