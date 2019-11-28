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

    # t_end=100
    # Yfind(h)=sqrt(h/(2h+1));
    # P=Yfind(H)

    # n_iter_Q=500;#50
    # Q_start=0.0
    # Q_end=2.5
    # Q_end=2.5
    # n_iter_P=500
    # P_start=P-.00002
    # P_end=P+.00002
    # # P_end=1
    # N=n_iter_Q*n_iter_P
    # t_end_mesh = t_end * ones(n_iter_Q,n_iter_P)
    # ArrP=range(P_start,stop=P_end,length=n_iter_P)
    # ArrQ=range(-Q_end,stop=Q_end,length=n_iter_Q)
    # mesh = [(Q, P) for P in ArrP, Q in ArrQ]
    # mesh_list = reshape(mesh, 1, :)


    # H = range(.185, stop = 0.235, length = 51)
    # H=[.12, .13, .15,.17,.18, .19, .195,.20,.205,.21,.23,.25]
    # t_grid=[30,  40, 50 ,  60,  70,  80,  90, 100]
    # H=[.13,.13,.13, .13, .13, .13, .13, .13]

    # H=[.13,.135,.14,.145]
    # H=[.17, .18, .19, .195,.2, .205,.21]
    # H = range(.185, stop = 0.205, length = 21)
    # H = range(.185, stop = 0.205, length = 21)
    # H = range(.19, stop = 0.20, length = 11)
    # H=[.18, .185, .19, .195, .2, .205,.21]
    # H = range(.185, stop = 0.205, length = 5)
    # H=[.19,.2]
    H = range(.13, stop = .155, length = 6)
    H = range(.1285, stop = .129, length = 6)
    H=.12925

  
    count = 1
    location="/home/brandon_behring/Desktop/MATLAB_FIGURES/"
    location_bson="~/Desktop/BSON_FILES/"
end
while count <= length(H)
    @everywhere t_end=100
    @everywhere P=Yfind(H[count])

    @everywhere     n_iter_Q=1000;#50
    # Q_start=.20
    # Q_end=2.5
    @everywhere  Q_end=5e-2
    # @everywhere  Q_start=-Q_end
    # @everywhere  Q_end=.00002
    @everywhere  Q_start=0.0
    # @everywhere  Q_end=1.1
 
    # @everywhere  n_iter_P=2001
    # @everywhere  P_start=-.5
    @everywhere  n_iter_P=1000
    @everywhere  P_start=P-6e-4
    @everywhere  P_end=P+5e-4
    # P_end=1
    @everywhere  N=n_iter_Q*n_iter_P
    @everywhere  t_end_mesh = t_end * ones(n_iter_Q,n_iter_P)
    @everywhere  ArrP=range(P_start,stop=P_end,length=n_iter_P)
    @everywhere  ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)

    @everywhere  mesh = [(Q, P) for P in ArrP, Q in ArrQ]
    @everywhere  mesh_list = reshape(mesh, 1, :)
@everywhere Energy_mesh = H[count] * ones(n_iter_Q,n_iter_P)
@everywhere h=replace(@sprintf("%.15f",H[count]),"."=>"_")
@everywhere nQ = @sprintf("_%5d",n_iter_Q)
@everywhere Q_win = @sprintf("_%.3f",Q_end)
@everywhere P_win = @sprintf("_%.3f",P_end)
@everywhere tend= @sprintf("_%5d",t_end)
println(tend)
@everywhere tend_T= @sprintf("%4f",t_end)
@everywhere right_now = string(Dates.Time(Dates.now()))
@everywhere h_BSON=   "BandF_Dimer_lagrangian_descriptors_data"*h*nQ*tend*Q_win*P_win*right_now*".bson"
@everywhere file_name=location* "Lagrangian_Descriptor_Dimer_"*h*nQ*tend*Q_win*P_win*right_now*".fig"
@everywhere h_title = @sprintf("h= %.6f",H[count])*" with tend="*tend_T*" and "*nQ


LD= @showprogress pmap(Dimer_Lagrangian_Descriptor.Dimer_Lagrangian_Descriptor_Function, mesh, Energy_mesh, t_end_mesh)

SAVE_DATA=Dict("LD"=>LD, "mesh"=>mesh, "ArrP"=>ArrP, "ArrQ"=>ArrQ, "Energy"=>H[count])
@save h_BSON SAVE_DATA

mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"


# mat"plot($P_0,$Q_0 ,'b.','MarkerSize',3)"
# mat"plot($P_1,$Q_1 ,'r.','MarkerSize',3)"
# mat"plot($P_2,$Q_2 ,'g.','MarkerSize',3)"
# mat"plot($P_4,$Q_4 ,'c.','MarkerSize',3)"
# mat"plot($P_5,$Q_5 ,'y.','MarkerSize',3)"

clims=[900 1300]
# mat"imagesc([0,$Q_end],[0,$P_end ],$LD, $clims)"
# mat"imagesc([0,-$Q_end],[0,$P_end ],$LD, $clims)"
# mat"imagesc([0,$Q_end],[0,-$P_end ],$LD, $clims)"
# mat"imagesc([0,-$Q_end],[0,-$P_end ],$LD, $clims)"

# mat"imagesc([0, $Q_end ],[0, $P_end],$LD)"
# mat"imagesc([0,-$Q_end ],[0, $P_end],$LD)"
# mat"imagesc([0,-$Q_end ],[0,-$P_end],$LD)"
# mat"imagesc([0, $Q_end ],[0,-$P_end],$LD)"

mat"imagesc([0, $Q_end ],[$P_start, $P_end],$LD)"
mat"imagesc([0,-$Q_end ],[$P_start, $P_end],$LD)"

# mat"imagesc([0, $Q_end ],[$P_start, $P_end],$LD,$clims)"
# mat"imagesc([0,-$Q_end ],[$P_start, $P_end],$LD,$clims)"

# mat"imagesc([$Q_start, $Q_end ],[$P_start,$P_end],$LD)"
# mat"imagesc([$Q_start, $Q_end ],[$P_start,$P_end],$LD)"  

# mat"imagesc([$Q_start,$Q_end ],[-$P_end,0],$LD)"
# mat"imagesc([$Q_start,$Q_end ],[0,$P_end],$LD)"
# mat"imagesc([$Q_start,$Q_end ],[0,-$P_end],$LD)"


# mat"imagesc([$Q_start,$Q_end ],[$P_start,$P_end],$LD, $clims)"
# mat"imagesc([$Q_start,$Q_end ],[$P_start,$P_end],$LD, $clims)"

# mat"imagesc([$Q_start,$Q_end ],[0,$P_end],$LD, $clims)"
# mat"imagesc([$Q_end,$Q_start ],[0,-$P_end],$LD, $clims)"

# mat"imagesc([$Q_start,$Q_end ],[0,$P_end],$LD)"
# mat"imagesc([$Q_end,$Q_start ],[0,-$P_end],$LD)"


# mat"imagesc([$Q_start,$Q_end ],[0,$P_end],$LD)"
# mat"imagesc([$Q_start,$Q_end ],[0,-$P_end],$LD)"


mat"colorbar"
Q_fix=sqrt(6)/3
mat"plot($Q_fix,0,'b.','MarkerSize',30)"
mat"plot(-$Q_fix,0,'b.','MarkerSize',30)"

mat"plot(0,$P,'r.','MarkerSize',30)"
mat"plot(0,-$P,'r.','MarkerSize',30)"
# mat"axis([ -$Q_end,$Q_end,-$P_end,$P_end ])"

# mat"axis([ $Q_start, $Q_end,$P_start,$P_end ])"
# mat"axis([ -$Q_end,$Q_end,-$P_end,$P_end ])"
mat"axis([ -$Q_end,$Q_end,$P_start,$P_end ])"
# mat"axis([ -$Q_end,$Q_end,-$P_end,$P_end ])"
# mat"axis([ $Q_start,$Q_end,-$P_end,$P_end ])"
mat"savefig($file_name)"
# mat"axis([ $Q_start,$Q_end,-$P_end,$P_end ])"
# mat"axis([ $Q_start,$Q_end,$P_start,$P_end ])"
# mat"close"
@everywhere global count += 1


end