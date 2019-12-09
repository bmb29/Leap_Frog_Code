using ProgressMeter
using Printf
using MATLAB
using BSON: @save, @load
using Dates
include("Dimer_Lagrangian_Descriptor.jl")


@everywhere begin
    include("Dimer_Lagrangian_Descriptor.jl")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_Code")
    using ProgressMeter
    using Printf
    using MATLAB
    using .Dimer_Lagrangian_Descriptor
    using BSON: @save, @load
    using Dates
    # H=range(.11,stop=0.145,length=107)
    # Energy=0.112012
    # Energy=.125789
    # Energy=.12501
    # Energy=.12501
    # Energy=.1249
    # Energy=.125
    Yfind(h)=sqrt(h/(2h+1));

    
    H = range(.13, stop = .1317, length = 7)
    # H = range(.11, stop = .21, length = 15)

    count = 1
    right_now = replace(replace(replace(string(Dates.now()),"."=>"_"),":"=>"_"),"-"=>"_")

end
while count <= length(H)
    # @everywhere t_end=15
      @everywhere t_end=120

    @everywhere P=Yfind(H[count])

    @everywhere     N=1000;#50
    # Q_start=.20
    # Q_end=2.5
    # @everywhere  Q_end=5e-2
    # @everywhere  Q_end=2.25
    @everywhere  Q_end=5e-2
    @everywhere  Q_start=-Q_end
    # @everywhere  Q_end=.00002
    @everywhere  Q_start=0.0
 
    # @everywhere  n_iter_P=2001
    # @everywhere  P_start=-.5
    @everywhere  n_iter_P=N
    @everywhere  n_iter_Q=N
    @everywhere  P_start=P-5e-3
    @everywhere  P_end=P+5e-3

    # @everywhere  P_start=0.0
    # @everywhere  P_end=.45
    # # P_end=1
    @everywhere  t_end_mesh = t_end * ones(n_iter_Q,n_iter_P)
    @everywhere  ArrP=range(P_start,stop=P_end,length=n_iter_P)
    @everywhere  ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)
    @everywhere location_fig="/media/brandon_behring/Extra_Space/MATLAB_FIGURES/"
    @everywhere location_bson="/media/brandon_behring/Extra_Space/BSON_FILES/"
    @everywhere  mesh = [(Q, P) for P in ArrP, Q in ArrQ]
    @everywhere  mesh_list = reshape(mesh, 1, :)
@everywhere Energy_mesh = H[count] * ones(n_iter_Q,n_iter_P)
@everywhere h=replace(@sprintf("%.15f",H[count]),"."=>"_")
@everywhere nQ = @sprintf("_%5d",n_iter_Q)
@everywhere Q_win = @sprintf("_%.3f",Q_end)
@everywhere P_win = @sprintf("_%.3f",P_end)
@everywhere tend= @sprintf("_%d",t_end)
@everywhere tend_T= @sprintf("%4f",t_end)
@everywhere h_BSON= location_bson * "BandF_Dimer_lagrangian_descriptors_data"*h*nQ*tend*Q_win*P_win*right_now*".bson"
@everywhere file_name1=location_fig* "Lagrangian_Descriptor_Dimer_1_"*h*nQ*tend*Q_win*P_win*right_now*".fig"
@everywhere file_name2=location_fig* "Lagrangian_Descriptor_Dimer_2_"*h*nQ*tend*Q_win*P_win*right_now*".fig"
@everywhere file_name3=location_fig* "Lagrangian_Descriptor_Dimer_3_"*h*nQ*tend*Q_win*P_win*right_now*".fig"


@everywhere h_title = @sprintf("h= %.6f",H[count])*" with tend="*tend_T*" and "*nQ
println(h_title)


LD= @showprogress pmap(Dimer_Lagrangian_Descriptor.Dimer_Lagrangian_Descriptor_Function, mesh, Energy_mesh, t_end_mesh)


dx=(Q_end-Q_start)/N
dy=(P_end-P_start)/N

gM=Dimer_Lagrangian_Descriptor.gradient_matrix_4(LD,dx,dy)
gradM_X,gradM_Y=Dimer_Lagrangian_Descriptor.gradient_matrix(LD,dx,dy)


Log_gradM_X=log10.(gradM_X)
Log_gradM_Y=log10.(gradM_Y)

biggest_X=maximum(filter(!isnan,Log_gradM_X))
biggest_Y=maximum(filter(!isnan,Log_gradM_Y))
smallest_X=minimum(filter(!isnan,Log_gradM_X))
smallest_Y=minimum(filter(!isnan,Log_gradM_Y))


gradM=log10.(gM)
# gradM=10 .^gradM
SAVE_DATA=Dict("LD"=>LD,"gradM"=>gradM,"Q_end"=>Q_end,"Q_start"=>Q_start, "P_start"=>P_start, "P_end"=>P_end,"N"=>N,"Energy"=>H[count],"ArrP"=>ArrP, "ArrQ"=>ArrQ,"t_end"=>t_end)
@save h_BSON SAVE_DATA

mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1000]); hold on;"
mat"title($h_title)"
mat"axis([ -$Q_end,$Q_end,$P_start,$P_end ])"

# mat"imagesc([$Q_start , $Q_end ],[$P_start, $P_end],$gradM, $clims)"
# mat"imagesc([-$Q_start, -$Q_end ],[-$P_start, -$P_end], $gradM, $clims)"
mat"imagesc([$Q_start , $Q_end ],[$P_start, $P_end],$LD)"
mat"imagesc([-$Q_start, -$Q_end ],[$P_start, $P_end], $LD)"
mat"imagesc([$Q_start, $Q_end ],[-$P_start, -$P_end], $LD)"
mat"imagesc([-$Q_start, -$Q_end ],[-$P_start, -$P_end], $LD)"
# mat"axis([ -$Q_end,$Q_end,-$P_end,$P_end ])"

# mat"imagesc([-$Q_start, -$Q_end ],[-$P_start, -$P_end], $gradM)"
# mat"imagesc([0,-$Q_end ],[$P_start, $P_end],$gradM, $clims)"
# mat"colormap winter"
# mat" cmocean('balance')"
mat"colormap(twilight)"

mat"colorbar"
Q_fix=sqrt(6)/3
mat"plot($Q_fix,0,'b.','MarkerSize',30)"
mat"plot(-$Q_fix,0,'b.','MarkerSize',30)"

mat"plot(0,$P,'r.','MarkerSize',10)"
mat"plot(0,-$P,'r.','MarkerSize',30)"
# mat"axis([ -$Q_end,$Q_end,$P_start,$P_end ])"
# mat"axis([ $Q_start,$Q_end,-$P_end,$P_end ])"
# mat"axis([ $Q_start, $Q_end,$P_start,$P_end ])"

mat"savefig($file_name1)"





mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1000]); hold on;"
mat"title($h_title)"
# mat"imagesc([$Q_start , $Q_end ],[$P_start, $P_end],$gradM, $clims)"
# mat"imagesc([-$Q_start, -$Q_end ],[-$P_start, -$P_end], $gradM, $clims)"
mat"imagesc([$Q_start , $Q_end ],[$P_start, $P_end],$LD)"
mat"imagesc([-$Q_start, -$Q_end ],[$P_start, $P_end], $LD)"
mat"imagesc([$Q_start, $Q_end ],[-$P_start, -$P_end], $LD)"
mat"imagesc([-$Q_start, -$Q_end ],[-$P_start, -$P_end], $LD)"
mat"axis([ -$Q_end,$Q_end,$P_start,$P_end ])"


# mat"imagesc([-$Q_start, -$Q_end ],[-$P_start, -$P_end], $gradM)"
# mat"imagesc([0,-$Q_end ],[$P_start, $P_end],$gradM, $clims)"
# mat"colormap winter"
# mat"colormap(twilight)"
mat"colormap(brewermap([],'PRGn'))"

mat"colorbar"
Q_fix=sqrt(6)/3
mat"plot($Q_fix,0,'b.','MarkerSize',30)"
mat"plot(-$Q_fix,0,'b.','MarkerSize',30)"

mat"plot(0,$P,'r.','MarkerSize',10)"
mat"plot(0,-$P,'r.','MarkerSize',30)"
# mat"axis([ $Q_start, $Q_end,$P_start,$P_end ])"


mat"savefig($file_name2)"


@everywhere global count += 1

end