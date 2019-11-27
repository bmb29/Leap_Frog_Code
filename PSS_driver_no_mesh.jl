using ProgressMeter
using MATLAB
using Printf
include("PSS_function.jl")



H=range(.11,stop=0.145,length=1001)
@showprogress 1 "Computing..." for Energy in H

mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;pbaspect([1.3 1 1]);"
location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/PSS_color_movie/"
h=replace(@sprintf("%.15f",Energy),"."=>"_")
file_name=location*h*".fig"
println(file_name)
mat"axis([ -.025,.025,-.3,.3 ])"
location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/"


n_iter_Q=51#50
Q_end=.3
n_iter_P=n_iter_Q;#51
P_end=.025
t_end=1e3
global max_hit=1000
#create an empty list to store dH
H_differences=zeros(0)
#julia's version of linspace
ArrP=range(0,stop=P_end,length=n_iter_P)
ArrQ=range(0,stop=Q_end,length=n_iter_Q)

#defining colors for PSS
COLOR=["#393b79" ,"#5254a3","#6b6ecf","#9c9ede" ,"#637939","#8ca252" ,"#b5cf6b" ,"#cedb9c" ,"#8c6d31","#bd9e39" ,"#e7ba52","#e7cb94","#843c39","#ad494a" ,"#d6616b","#e7969c" ,"#7b4173" ,"#a55194","#ce6dbd" ,"#de9ed6"]
#
for j=1:n_iter_P
    Q1,P1,dH=PSS_function(ArrQ[j],0, Energy, t_end, max_hit)
    Q2,P2,dH=PSS_function(-ArrQ[j],0, Energy, t_end, max_hit)
    Q3,P3,dH=PSS_function(0,ArrP[j], Energy, t_end, max_hit)
    Q4,P4,dH=PSS_function(0,-ArrP[j], Energy, t_end, max_hit)
    Q=vcat(Q1,Q2)
    P=vcat(P1,P2)
    current_color=COLOR[mod(2*j,length(COLOR))+1]
    Q_2=vcat(Q3,Q4)
    P_2=vcat(P3,P4)
    current_color_2=COLOR[mod(j+5,length(COLOR))+1]
    mat"plot($P_2,$Q_2,'.','MarkerSize',.1,'color',$current_color_2); hold on;"
    mat"plot($P,$Q,'.','MarkerSize',.1,'color',$current_color); hold on;"
end


n_iter_Q=13;#50
Q_start=-.3
Q_end=.3
n_iter_P=n_iter_Q;#51
P_start=-.025
P_end=.025
#create an empty list to store dH
H_differences=zeros(0)
#julia's version of linspace
ArrP=range(0,stop=P_end,length=n_iter_P)
ArrQ=range(0,stop=Q_end,length=n_iter_Q)

for j=1:n_iter_P
    for k=1:n_iter_Q
        Q1,P1,dH=PSS_function(ArrQ[k], ArrP[j], Energy, t_end, max_hit)
        Q2,P2,dH=PSS_function(-ArrQ[k], ArrP[j], Energy, t_end, max_hit)
        Q3,P3,dH=PSS_function(ArrQ[k], -ArrP[j], Energy, t_end, max_hit)
        Q4,P4,dH=PSS_function(-ArrQ[k], -ArrP[j], Energy, t_end, max_hit)
        Q=vcat(Q1,Q2)
        P=vcat(P1,P2)
        Q_2=vcat(Q3,Q4)
        P_2=vcat(P3,P4)
        current_color=COLOR[mod(j+k,length(COLOR))+1]
        current_color_2=COLOR[mod(j+5,length(COLOR))+1]
        #
        mat"plot($P,$Q,'.','MarkerSize',.1,'color',$current_color); hold on;"
        mat"plot($P_2,$Q_2,'.','MarkerSize',.1,'color',$current_color_2); hold on;"

        # mat"plot($P,$Q,'k.','MarkerSize',.1)"

    end
end
#
mat"savefig($file_name)"
mat"close"
end
