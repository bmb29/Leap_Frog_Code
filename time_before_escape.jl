

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
include("escape_num.jl")


@everywhere begin
    include("escape_num.jl")
    # push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_2019")
    using ProgressMeter
    using Printf
    using MATLAB
    using .escape_num
    t_end = exp(10)
    # width = sqrt(3); height = 1
    width = 1; height =1
    n_iter_Q = 1000;n_iter_P = n_iter_Q; N = n_iter_P * n_iter_Q;
    ArrP = range(-width, stop = width, length = n_iter_P)
    ArrQ = range(-height, stop = height, length = n_iter_Q)
    mesh_P=[P for Q in ArrQ, P in ArrP]
    mesh_Q=[Q for Q in ArrQ, P in ArrP]
    mesh = [(P, Q) for Q in ArrQ, P in ArrP]
    mesh_list = reshape(mesh, 1, :)
    t_end = t_end * ones(n_iter_Q,n_iter_P);
    H=.22
    location = "/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
    # location = "/Users/brandonbehring/Desktop/"
end


    @everywhere Energy = H* ones(n_iter_Q,n_iter_P);
    # @everywhere Energy = H* ones(N);

    @everywhere h = replace(@sprintf("%.13f",H), "." => "_")
    @everywhere file_name = location * "time_before_Escape_" * h * ".fig"
    println(file_name)
    num_until_exit = @showprogress pmap(escape_num.escape_exit_num, mesh, t_end, Energy)

    logz=log1p.(num_until_exit)
    mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
    # # mat"plot($P_0,$Q_0 ,'b.','MarkerSize',3)"
    # # mat"plot($P_1,$Q_1 ,'r.','MarkerSize',3)"
    # # mat"plot($P_2,$Q_2 ,'g.','MarkerSize',3)"
    # # mat"plot($P_4,$Q_4 ,'c.','MarkerSize',3)"
    # # mat"plot($P_5,$Q_5 ,'y.','MarkerSize',3)"
    mat"imagesc([-$width,$width],[-$height,$height ],$logz)"
    mat"colorbar"

    mat"axis([ -$width,$width,-$height,$height ])"
    # mat"axis([ -1,1,-1,1])"
    mat"savefig($file_name)"

