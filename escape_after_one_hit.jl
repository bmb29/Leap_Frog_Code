using ProgressMeter
using MATLAB
using Printf
using BSON: @save, @load
using PyCall
pygui(:qt)
using PyPlot
pygui(true)
include("escape_dimer_one_hit.jl")

H=.25

t_end = 1e3
width = 2.25; height = 4.5
n_iter_P = 1001;n_iter_Q =1001;
N = n_iter_P * n_iter_Q;
ArrQ = range(-width, stop = width, length = n_iter_P)
ArrP = range(0, stop = height, length = n_iter_Q)
mesh = [(Q, P) for Q in ArrQ, P in ArrP]
mesh_list = reshape(mesh, 1, :)
t_end_list = t_end * ones(N);
energy_list= H* ones(N);

h = replace(@sprintf("%.13f",H), "." => "_")
nQ = @sprintf("_%d",n_iter_Q)
nP = @sprintf("_%d",n_iter_P)
h_BSON="num_hhits_escape_data"*h*nQ*nP*".bson"
location = "/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
println(h_BSON)

h_title = @sprintf("h= %.13f",H)

num_until_exit = @showprogress map(escape_exit_num_dimer, mesh_list, t_end_list, energy_list)

one_hit_exit=[]
for i = 1:N
    Q, P = mesh_list[i]
    if num_until_exit[i] == 1
        push!(one_hit_exit,[Q,P])
    end
end
t_end=1e3
for k=1:3
old_list=copy(one_hit_exit)
energy=.25
radius=.01

@showprogress for (index,point) in enumerate(old_list)
    add_four!(point,radius, one_hit_exit, t_end, H)
end


Q,P=unzipper(one_hit_exit)
Q_one=vcat(Q,-Q); P_one=vcat(P,-P);

plot(Q_one,P_one, ".",markersize=1,c=:blue);
println(length(one_hit_exit))
SAVE_DATA=Dict("one_hit_exit"=>one_hit_exit)
@save h_BSON SAVE_DATA

end

