using ProgressMeter
using Printf
using MATLAB
using BSON: @save, @load

include("Aref_Lagrangian_Descriptor.jl")


@everywhere begin
    include("Aref_Lagrangian_Descriptor.jl")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_Code")
    using ProgressMeter
    using Printf
    using MATLAB
    using .Aref_Lagrangian_Descriptor
    using BSON: @save, @load
    # H=range(.11,stop=0.145,length=107)
    # Energy=0.112012
    # Energy=.125789
    Energy=.129
    # Energy=.12501
    # Energy=.12501
    # Energy=.1249
    # Energy=.125
    t_typical=Aref_Lagrangian_Descriptor.final_T(Energy)
    # file_name=location*h*".fig"
    #defining colors for PSS
    n_iter_Q=1000;#50
    # Q_start=.20
    P_end=3.0*1e-3
    P_start=-P_end
    # Q_end=.7
    n_iter_P=1000
    # P_start=
    # P_end=.2
    Q_end=.035
    Q_start=-Q_end
    t_end=t_typical
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
@everywhere nQ = @sprintf("_%d",n_iter_Q)
@everywhere tend= @sprintf("_%d",t_end)
@everywhere tend_T= @sprintf("%4f",t_end)
@everywhere h_BSON="lagrangian_descriptors_data"*h*nQ*tend*".bson"
@everywhere file_name=location*"MATLAB_LD_"*h*nQ*tend*".fig"
@everywhere h_title = @sprintf("h= %.6f",Energy)*"with tend="*tend_T


LD= @showprogress pmap(Aref_Lagrangian_Descriptor.Aref_Lagrangian_Descriptor_Function,mesh, Energy_mesh, t_end_mesh)
SAVE_DATA=Dict("LD"=>LD, "mesh"=>mesh, "ArrP"=>ArrP, "ArrQ"=>ArrQ, "Energy"=>Energy)
@save h_BSON SAVE_DATA
mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"
biggest=maximum(filter(!isnan,LD))
smallest=minimum(filter(!isnan,LD))

# mat"plot($P_0,$Q_0 ,'b.','MarkerSize',3)"
# mat"plot($P_1,$Q_1 ,'r.','MarkerSize',3)"
# mat"plot($P_2,$Q_2 ,'g.','MarkerSize',3)"
# mat"plot($P_4,$Q_4 ,'c.','MarkerSize',3)"
# mat"plot($P_5,$Q_5 ,'y.','MarkerSize',3)"

# clims=[smallest 630]
# mat"imagesc([0,$P_end],[0,$Q_end ],$LD, $clims)"
# mat"imagesc([0,-$P_end],[0,$Q_end ],$LD, $clims)"
# mat"imagesc([0,$P_end],[0,-$Q_end ],$LD, $clims)"
# mat"imagesc([0,-$P_end],[0,-$Q_end ],$LD, $clims)"


mat"imagesc([$P_start,$P_end],[$Q_start,$Q_end ],$LD)"
# mat"imagesc([0,-$P_end],[0,$Q_end ],$LD)"
# mat"imagesc([0,$P_end],[0,-$Q_end ],$LD)"
# mat"imagesc([0,-$P_end],[0,-$Q_end ],$LD)"

# mat"imagesc([0,$P_end],[0,$Q_end ],$LD)"
# mat"imagesc([0,-$P_end],[0,$Q_end ],$LD)"
# mat"imagesc([0,$P_end],[0,-$Q_end ],$LD)"
# mat"imagesc([0,-$P_end],[0,-$Q_end ],$LD)"

mat"colorbar"
mat"axis([ -$P_end,$P_end,-$Q_end,$Q_end ])"
# mat"axis([ $P_start,$P_end,$Q_start,$Q_end ])"
mat"savefig($file_name)"
