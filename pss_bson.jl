using DifferentialEquations
using ProgressMeter
using PyCall
using MATLAB
using Printf
using BSON: @save, @load
# pygui(:qt)
using PyPlot
include("PSS_Definitions_Dimer_X.jl")
# pygui(true);

H=range(.08,stop=0.15,length=25)
# H=.075
@showprogress for Energy in H
location_fig="/media/brandon_behring/Extra_Space/MATLAB_FIGURES/"
location_BSON="/media/brandon_behring/Extra_Space/BSON_FILES/"
h=replace(@sprintf("%.15f",Energy),"."=>"_")

t_end=1e5;
N_iter_Q=25 
Q_start=0.0
Q_end=.9
N_iter_P=25
P_start=0.0
P_end=.45

# t_end=1e5;
# N_iter_Q=50 
# Q_start=0.1
# Q_end=.3
# N_iter_P=50
# P_start=0.0
# P_end=.3


nQ=@sprintf("_%d",N_iter_Q)
nP=@sprintf("_%d",N_iter_P)
h_title=@sprintf("%.8f",Energy)
println(h_title)
file_name=location_fig*"DIMER_PSS_X"*h*nQ*nP*".fig"
h_BSON=location_BSON*"DIMER_PSS_X"*h*nQ*nP*".bson"

ArrP=range(P_start,stop=P_end,length=N_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=N_iter_Q)
Brandons_Colors=["#393b79" ,"#5254a3","#6b6ecf","#9c9ede" ,"#637939","#8ca252" ,"#b5cf6b" ,"#cedb9c" ,"#8c6d31","#bd9e39" ,"#e7ba52","#e7cb94","#843c39","#ad494a" ,"#d6616b","#e7969c" ,"#7b4173" ,"#a55194","#ce6dbd" ,"#de9ed6"];
SAVE_DATA=Dict()
# figure()
@showprogress for k=1:N_iter_Q
    for j=1:N_iter_P
        Q_PSS,P_PSS=PSS_function(ArrQ[k], ArrP[j], Energy, t_end);
        SAVE_DATA[(ArrQ[k],ArrP[j])]=Q_PSS, P_PSS
        # current_color=Brandons_Colors[mod(k,length(Brandons_Colors))+1]
    end
end
# @showprogress for k=1:N_iter_P
#     Q_PSS,P_PSS=PSS_function(0,ArrP[k], Energy, t_end);
#     current_color=Brandons_Colors[mod(k,length(Brandons_Colors))+1]
#     SAVE_DATA[k+N_iter_Q]=Q_PSS, P_PSS
#     if Q_PSS!=0
#         # mat"plot($Q_PSS,$P_PSS,'.','MarkerSize',3,'color',$current_color); hold on;"
#         # mat"plot(-$Q_PSS,-$P_PSS,'.','MarkerSize',3,'color',$current_color); hold on;"
#     end
# end

@save h_BSON SAVE_DATA
Yfind(h)=sqrt(h/(2h+1))
P=Yfind(Energy)
mat"figure();set(gcf, 'Position',  [0, 0, 2500, 2000]); hold on;"
mat"xlabel('Q')"
mat"ylabel('P')"
mat"title($h_title)"
mat"axis([-.3, .3, -.3, .3])"
# mat"axis([-1, 1, -.5, .5])"
@showprogress for k=1:N_iter_Q
    for j=1:N_iter_P
        Q_PSS,P_PSS= SAVE_DATA[(ArrQ[k],ArrP[j])]
        current_color=Brandons_Colors[mod(k*j,length(Brandons_Colors))+1]
        mat"plot($Q_PSS,$P_PSS,'.','MarkerSize',2,'color',$current_color); hold on;"
        # mat"plot(-$Q_PSS,$P_PSS,'.','MarkerSize',2,'color',$current_color); hold on;"
        # mat"plot($Q_PSS,-$P_PSS,'.','MarkerSize',2,'color',$current_color); hold on;"
        mat"plot(-$Q_PSS,-$P_PSS,'.','MarkerSize',2,'color',$current_color); hold on;"
    end
end
Q_fix=sqrt(6)/3
mat"plot($Q_fix,0,'b.','MarkerSize',30)"
mat"plot(-$Q_fix,0,'b.','MarkerSize',30)"
mat"plot(0,$P,'r.','MarkerSize',30)"
mat"plot(0,-$P,'r.','MarkerSize',30)"
mat"savefig($file_name)"
end
