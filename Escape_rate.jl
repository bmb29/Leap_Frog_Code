using DifferentialEquations
using MATLAB
using ProgressMeter
using BSON: @save, @load


n_iter = 100
barrier = 2
H_start=.128
H_end=.13
HH = range(H_start, stop = H_end, length = n_iter)
exit_time = zeros(n_iter)
max_time = 1e6
tspan = (0.0, max_time)
length_max=100000

ODE1(z, w) = conj(im * w .* ( 1 ./ (z.^2 - w.^2) + 1 ./ (1 + w.^2) ))
ODE2(z, w) = conj(im * z .* ( 1 ./ (w.^2 - z.^2) + 1 ./ (1 + z.^2) ))
condition_escape(u, t, integrator) = u[5]>length_max || maximum([abs(u[1]),abs(u[2]),abs(u[3]),abs(u[4])]) > barrier 
function condition_hits_X(u, t, integrator) # Event when event_f(u,t) == 0
    u[4] 
end

function affect_update_iterator_X!(integrator)
    integrator.u[5] = integrator.u[5] + 1
end 
affect_stop!(integrator) = terminate!(integrator)


callback_escape = DiscreteCallback(condition_escape, affect_stop!)
callback_hits_PSS_X = ContinuousCallback(condition_hits_X, affect_update_iterator_X!, rootfind = false)
cb = CallbackSet(callback_hits_PSS_X, callback_escape)
function Eq_of_M(du, u, p, t)
    du[1] = real(ODE1(u[1] + im * u[2], u[3] + im * u[4]))
    du[2] = imag(ODE1(u[1] + im * u[2], u[3] + im * u[4]))
    du[3] = real(ODE2(u[1] + im * u[2], u[3] + im * u[4]))
    du[4] = imag(ODE2(u[1] + im * u[2], u[3] + im * u[4]))
    du[5] = 0
    return
end
# refernce
# z=X+iP W=Q+iY . Plane Q=0, P=0
# z=u[1]+i u[2] ; w=u[3]+i u[4]
function Yfind_Aref(Q, P, H_Aref)
    H = (2 * H_Aref)^2
    try
        Y = sqrt((1 + (-1) * H * ((-1) + P^2)^2)^(-1) * (P^2 + (-1) * Q^2 + H * ((-1) + P^2)^2 * ((-1) + Q^2) + sqrt(H * ((-1) + P^2)^4 + (-4) * ((-1) + H * ((-1) + P^2)^2)
        * ((-1) * P^2 + H * ((-1) + P^2)^2) * Q^2)))
    catch
        Y = zeros(0)
    end
end
n_iter_Q = 4;#50
n_iter_P = 4
N = n_iter_P * n_iter_Q
Q_end = 1e-4
P_end = 1e-5

ArrP = range(0, stop = P_end, length = n_iter_P)
ArrQ = range(0, stop = Q_end, length = n_iter_Q)
mesh = [(Q, P) for P in ArrP, Q in ArrQ]
mesh_list = reshape(mesh, 1, :)
hits=zeros(n_iter)


@showprogress for j = 1:n_iter
    H = HH[j]
    for k = 1:N
        Q,P=mesh_list[k]
        # println(mesh_list[k])
        Y_0 = Yfind_Aref(Q, P, H)
        u0 = [0;P;Q;Y_0;0]
        prob = ODEProblem(Eq_of_M, u0, tspan)
        t,A=solve(prob, Tsit5(),reltol=1e-8,abstol=1e-11,maxiters=1e15,callback=cb)
        # t, A = solve(prob, Tsit5(), maxiters = 1e15, callback = cb)
        if t.t[end]<max_time && maximum([abs(A[1,end]),abs(A[2,end]),abs(A[3,end]),abs(A[4,end])]) > barrier 
            hits[j]=hits[j]+1
        end
    end
end
mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
# logT = log.(exit_time)
normalized_hits=hits/N
# mat"plot($HH,$logT,'b.','MarkerSize',10)"
mat"plot($HH,$normalized_hits,'b.','MarkerSize',10)"

mat"axis([ $H_start, $H_end, 0,1])"