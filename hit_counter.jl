

# @everywhere module_dir = "/Users/brandonbehring/Desktop/Leap_Frog_2019"
# @everywhere push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
# # @everywhere module_dir ="/Users/brandonbehring/Desktop/Leap_Frog_2019"
# # @everywhere push!(LOAD_PATH, $module_dir)
# # @everywhere thisDir = dirname(@__FILE__())
# # @everywhere any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)
# include("escape.jl")
# @everywhere using .escape

using ProgressMeter
using Printf
using MATLAB
using BSON: @save, @load

include("escape_dimer.jl")


@everywhere begin
    include("escape_dimer.jl")
    # push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_2019")
    using ProgressMeter
    using Printf
    using MATLAB
    using .escape_dimer
    using BSON: @save, @load

    t_end = 1e5
    width = 2.25; height = 4.5
    n_iter_P = 4000;n_iter_Q =4000;
     N = n_iter_P * n_iter_Q;
    ArrP = range(-width, stop = width, length = n_iter_P)
    ArrQ = range(0, stop = height, length = n_iter_Q)
    mesh_P=[P for Q in ArrQ, P in ArrP]
    mesh_Q=[Q for Q in ArrQ, P in ArrP]
    mesh = [(P, Q) for Q in ArrQ, P in ArrP]
    mesh_list = reshape(mesh, 1, :)
    t_end = t_end * ones(n_iter_Q,n_iter_P);
 
    # H = range(.2, stop = 0.21, length = 11)
    H=.3
    location = "/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"

    count = 1
end

while count <= length(H)
    println(count)
    println(H[count])
    @everywhere Energy = H[count] * ones(N)
    @everywhere h = replace(@sprintf("%.13f",H[count]), "." => "_")
    @everywhere nQ = @sprintf("_%d",n_iter_Q)
    @everywhere nP = @sprintf("_%d",n_iter_P)
    @everywhere h_title = @sprintf("h= %.13f",H[count])

    @everywhere h_BSON="num_hhits_escape_data"*h*nQ*nP*".bson"
 

    @everywhere file_name = location*"HIT_Escape_" * h * ".fig"
    println(file_name)
    num_until_exit = @showprogress pmap(escape_dimer.escape_exit_num_dimer, mesh_list, t_end, Energy)
    SAVE_DATA=Dict("escape_times"=>num_until_exit, "mesh"=>mesh, "ArrP"=>ArrP, "ArrQ"=>ArrQ, "Energy"=>H)
    @save h_BSON SAVE_DATA

    @everywhere begin
        Q_0 = zeros(0);P_0 = zeros(0)
        Q_1 = zeros(0);P_1 = zeros(0)
        Q_2 = zeros(0);P_2 = zeros(0)
        Q_3 = zeros(0);P_3 = zeros(0)
        Q_4 = zeros(0);P_4 = zeros(0)
        Q_5 = zeros(0);P_5 = zeros(0)
        Q_b = zeros(0);P_b = zeros(0)

    end
    for i = 1:N
        Q, P = mesh_list[i]
        if num_until_exit[i] == -1
            push!(Q_b, Q)
            push!(P_b, P)
        elseif num_until_exit[i] == 0
            push!(Q_0, Q)
            push!(P_0, P)
        elseif num_until_exit[i] == 1
            push!(Q_1, Q)
            push!(P_1, P)
        elseif num_until_exit[i] == 2
            push!(Q_2, Q)
            push!(P_2, P)
        elseif num_until_exit[i] == 3
            push!(Q_3, Q)
            push!(P_3, P)
        elseif num_until_exit[i] == 4
            push!(Q_4, Q)
            push!(P_4, P)
        elseif num_until_exit[i] == 5
            push!(Q_5, Q)
            push!(P_5, P)
        end
    end
  


    mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
    mat"title($h_title)"
    mat"axis([ -$width,$width,-$height,$height ])"
    mat"plot($Q_b,$P_b ,'k.','MarkerSize',3)"
    mat"plot(-$Q_b,-$P_b ,'k.','MarkerSize',3)"

    mat"plot($Q_0,$P_0 ,'b.','MarkerSize',3)"
    mat"plot(-$Q_0,-$P_0 ,'b.','MarkerSize',3)"

    mat"plot($Q_1,$P_1 ,'r.','MarkerSize',3)"
    mat"plot(-$Q_1,-$P_1 ,'r.','MarkerSize',3)"

    mat"plot($Q_2,$P_2 ,'g.','MarkerSize',3)"
    mat"plot(-$Q_2,-$P_2 ,'g.','MarkerSize',3)"

    mat"plot($Q_3,$P_3 ,'m.','MarkerSize',3)"
    mat"plot(-$Q_3,-$P_3 ,'m.','MarkerSize',3)"

    mat"plot($Q_4,$P_4 ,'c.','MarkerSize',3)"
    mat"plot(-$Q_4,-$P_4 ,'c.','MarkerSize',3)"

    mat"plot($Q_5,$P_5 ,'y.','MarkerSize',3)"
    mat"plot(-$Q_5,-$P_5 ,'y.','MarkerSize',3)"

    mat"savefig($file_name)"
    @everywhere global count += 1
end
