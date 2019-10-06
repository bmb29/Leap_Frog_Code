using ProgressMeter
using MATLAB
using Printf
include("PSS_function.jl")



# H=range(.11,stop=0.145,length=107)
Energy=0.112012
Energy=.13
# @showprogress 1 "Computing..." for Energy in H

mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
# mat"""
# set(gca,'XTick',[])
# set(gca,'YTick',[])
# set(gca,'xticklabel',[])
# set(gca,'yticklabel',[])
# set(gca,'XColor','none')
# set(gca,'YColor','none')
# set(gca,'XColor','none')
# set(gca,'YColor','none')
# """
location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/"
h=replace(@sprintf("%.15f",Energy),"."=>"_")
file_name=location*h*".fig"
file_name=h*".fig"
COLOR=["#393b79" ,"#5254a3","#6b6ecf","#9c9ede" ,"#637939","#8ca252" ,"#b5cf6b" ,"#cedb9c" ,"#8c6d31","#bd9e39" ,"#e7ba52","#e7cb94","#843c39","#ad494a" ,"#d6616b","#e7969c" ,"#7b4173" ,"#a55194","#ce6dbd" ,"#de9ed6"]
println(file_name)
mat"axis([ -2.5, 2.5,-1.5, 1.5])"
#defining colors for PSS
N_iter_Q=10;#50
Q_start=-.3
Q_end=.3
N_iter_P=10
P_start=-.01
P_end=.01
t_end=1e3
max_hit=1000
mat"axis([ $P_start, $P_end, $Q_start, $Q_end])"
#create an empty list to store dH
H_differences=zeros(0)
#julia's version of linspace
ArrP=range(P_start,stop=P_end,length=N_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=N_iter_Q)

@showprogress 1 "Computing..." for k=1:N_iter_P
    for j=1:N_iter_Q
        Q1,P1,dH=PSS_function(ArrQ[j], ArrP[k], Energy, t_end, max_hit)
            # Q2,P2,dH=PSS_function(-ArrQ[k], ArrP[j], Energy, t_end, max_hit)
            # Q3,P3,dH=PSS_function(ArrQ[k], -ArrP[j], Energy, t_end, max_hit)
            # Q4,P4,dH=PSS_function(-ArrQ[k], -ArrP[j], Energy, t_end, max_hit)
            # Q=vcat(Q1,Q2)
            # P=vcat(P1,P2)
            # Q_2=vcat(Q3,Q4)
            # P_2=vcat(P3,P4)
            if dH<1e-10
                current_color=COLOR[mod(k,length(COLOR))+1]
                mat"plot($P1,$Q1,'.','MarkerSize',.1,'color',$current_color); hold on;"
            end
       
    end
end

#
mat"savefig($file_name)"
# mat"close"
# end
