

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
using LinearAlgebra
using BSON: @save, @load
include("escape_num_dimer.jl")
    # push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_2019")
using ProgressMeter
using Printf
using MATLAB
using .escape_num_dimer

Energy=range(.3,stop=.49,length=20);
t_end = exp(5)
    # width = sqrt(3); height = 1
width = 2.25; height = 3
n_iter_P = 375;n_iter_Q =500; N = n_iter_P * n_iter_Q;
ArrP = range(-width, stop = width, length = n_iter_P)
ArrQ = range(0, stop = height, length = n_iter_Q)
mesh_P = [P for Q in ArrQ, P in ArrP]
mesh_Q = [Q for Q in ArrQ, P in ArrP]
mesh = [(P,Q ) for Q in ArrQ, P in ArrP]



@showprogress for H in Energy
println(H)
location = "/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
h = replace(@sprintf("%.8f_",H), "." => "_")
nQ = @sprintf("%d",n_iter_Q)
file_name = location * "Dimer_time_before_Escape_" * h * nQ * ".fig"
h_BSON = location*"escape_data" * h * nQ * ".bson"
# num_until_exit = @showprogress pmap(escape_num_dimer.escape_exit_num_dimer, mesh, t_end, Energy)

# logz = log1p.(num_until_exit)

# SAVE_DATA = Dict("escape_times" => num_until_exit, "mesh" => mesh, "ArrP" => ArrP, "ArrQ" => ArrQ, "Energy" => H)
# @save h_BSON SAVE_DATA


SAVE_DATA=Dict()

# escape_time=zeros(n_iter_Q,n_iter_P)
escape_time=zeros(n_iter_Q,n_iter_P)

# figure()
# @showprogress for k=1:N
#     escape_time[k]=escape_num_dimer.escape_exit_num_dimer( mesh_list[k], t_end, H)
#     SAVE_DATA[mesh_list[k]]= escape_time[k]

# end
@showprogress for k=1:n_iter_Q
    for j=1:n_iter_P
        Q2=mesh[k,j][1]; P2=mesh[k,j][2];
        p1=[Q2,P2]
        p2=[0,1]
        p3=[0,-1]
        tol=.1
        if norm(p1)<tol|| norm(p1-p2)<tol || norm(p1-p3)<tol
            escape_time[k,j]=t_end
        end
        escape_time[k,j]=escape_num_dimer.escape_exit_num_dimer( mesh[k,j], t_end, H)
        SAVE_DATA[(ArrQ[k],ArrP[j])]=escape_time[k,j]
        # current_color=Brandons_Colors[mod(k,length(Brandons_Colors))+1]
    end
end
@save h_BSON SAVE_DATA

logz =log1p.(escape_time)
h_title=@sprintf("%.8f",H)

mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"
    # # # # mat"plot($P_0,$Q_0 ,'b.','MarkerSize',3)"
    # # # # mat"plot($P_1,$Q_1 ,'r.','MarkerSize',3)"
    # # # # mat"plot($P_2,$Q_2 ,'g.','MarkerSize',3)"
    # # # # mat"plot($P_4,$Q_4 ,'c.','MarkerSize',3)"
    # # # # mat"plot($P_5,$Q_5 ,'y.','MarkerSize',3)"
    # mat"imagesc([-$width,$width],[-$height,$height ],$logz)"
mat"imagesc([-$width,$width],[0,$height ],$logz)"
mat"imagesc([$width,-$width],[0,-$height ],$logz)"

mat"colorbar"

mat"axis([ -$width,$width,-$height,$height ])"
    # # # mat"axis([ -1,1,-1,1])"
mat"savefig($file_name)"


mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"
    # # # # mat"plot($P_0,$Q_0 ,'b.','MarkerSize',3)"
    # # # # mat"plot($P_1,$Q_1 ,'r.','MarkerSize',3)"
    # # # # mat"plot($P_2,$Q_2 ,'g.','MarkerSize',3)"
    # # # # mat"plot($P_4,$Q_4 ,'c.','MarkerSize',3)"
    # # # # mat"plot($P_5,$Q_5 ,'y.','MarkerSize',3)"
    # mat"imagesc([-$width,$width],[-$height,$height ],$logz)"
mat"contourf([-$width,$width],[0,$height ],$logz)"
mat"contourf([$width,-$width],[0,-$height ],$logz)"

mat"colorbar"

mat"axis([ -$width,$width,-$height,$height ])"
    # # # mat"axis([ -1,1,-1,1])"
end
