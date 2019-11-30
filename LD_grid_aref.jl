using ProgressMeter
using Printf
using MATLAB
using BSON: @save, @load

include("Aref_LD_grad_function.jl")


@everywhere begin
    include("Aref_LD_grad_function.jl")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_Code")
    using ProgressMeter
    using Printf
    using MATLAB
    using .Aref_LD_grad_function
    using BSON: @save, @load

    Energy=.129
    t_typical=Aref_LD_grad_function.final_T(Energy)
    n_iter_Q=3000;#50
    Q_end=.05
    Q_start=0
    n_iter_P=3000

    P_end=3.0*1e-3
    P_start=-P_end

    t_end=30*t_typical
    N=n_iter_Q*n_iter_P
    t_end_mesh = t_end * ones(n_iter_Q,n_iter_P)
    Energy_mesh = Energy* ones(n_iter_Q,n_iter_P)
    ArrP=range(P_start,stop=P_end,length=n_iter_P)
    ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)
    mesh = [(P, Q) for Q in ArrQ, P in ArrP]
    mesh_list = reshape(mesh, 1, :)
end
@everywhere location="/media/brandon_behring/Extra_Space/MATLAB_FIGURES/"
@everywhere h=replace(@sprintf("%.15f",Energy),"."=>"_")
@everywhere h=replace(@sprintf("%.15f",Energy),"."=>"_")
@everywhere nQ = @sprintf("_%d",n_iter_Q)
@everywhere tend= @sprintf("_%d",t_end)
@everywhere tend_T= @sprintf("%4f",t_end)
@everywhere h_BSON="lagrangian_descriptors_data"*h*nQ*tend*".bson"
@everywhere file_name=location*"MATLAB_LD_"*h*nQ*tend*".fig"
@everywhere h_title = @sprintf("h= %.6f",Energy)*"with tend="*tend_T


LD= @showprogress pmap(Aref_LD_grad_function.LD_Helper,mesh, Energy_mesh, t_end_mesh)
SAVE_DATA=Dict("LD"=>LD, "mesh"=>mesh, "ArrP"=>ArrP, "ArrQ"=>ArrQ, "Energy"=>Energy)
@save h_BSON SAVE_DATA

dx=(P_end-P_start)/N
dy=(Q_end-Q_start)/N

gradM=Aref_LD_grad_function.gradient_matrix_4(LD,dx,dy)
biggest=maximum(filter(!isnan,gradM))
smallest=minimum(filter(!isnan,gradM))


mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"

mat"imagesc([$P_start,$P_end ],[0,$Q_end],$LD)"
mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$LD)"

mat"colorbar"
mat"axis([-$P_end,$P_end,-$Q_end,$Q_end ])"


mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"

clims=[smallest biggest*.1]
mat"imagesc([$P_start,$P_end ],[0,$Q_end],$gradM, $clims)"
mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$gradM, $clims)"

mat"colorbar"
mat"axis([-$P_end,$P_end,-$Q_end,$Q_end ])"

mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"

clims=[smallest biggest*.01]
mat"imagesc([$P_start,$P_end ],[0,$Q_end],$gradM, $clims)"
mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$gradM, $clims)"

mat"colorbar"
mat"axis([-$P_end,$P_end,-$Q_end,$Q_end ])"


# mat"axis([ $P_start,$P_end,$Q_start,$Q_end ])"
mat"savefig($file_name)"


