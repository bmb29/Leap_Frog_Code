using DifferentialEquations
using ProgressMeter
using PyCall
using MATLAB
using Printf

# pygui(:qt)
using PyPlot
include("PSS_Definitions_Dimer.jl")
# pygui(true);
using Roots
function Q1_find_dimer(Q2,P2,H_dimer)
    Q1_to_find(Q1)=Hamiltonian_Dimer([Q1,Q2],[0,P2],1)-H_dimer
    try
        Q1=find_zero(Q1_to_find,.01,maxeval=100,maxfnevals=300,tol=1e-15)
     catch
        Q1=zeros(0)
    end
end

H=range(.1,stop=0.25,length=14)

@showprogress 1 "Computing..." for Energy in H
location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/"
h=replace(@sprintf("%.15f",Energy),"."=>"_")
h_title=@sprintf("%.8f",Energy)
file_name=location*"DIMER_PSS"*h*".fig"
file_name=h*".fig"
mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"xlabel('Q')"
mat"ylabel('P')"
mat"title($h_title)"
mat"axis([-3, 3, -4, 4])"
# Q1=Q1_find_dimer(0,.4,H)
# t_end=1000.0;
# q0,p0=[zeros(2) for i in 1:2]
# q0[1]=Q1; 
# q0[2]=.2; 
# p0[1]=0; 
# p0[2]=.4;
# Hamiltonian_Dimer(q0,p0,1)
# prob= HamiltonianProblem{true}(Hamiltonian_Dimer, q0, p0, (0., t_end));
# @time t,A=solve(prob, Vern9(),reltol=1e-13,abstol=1e-15,maxiters=1e15);

# q1=A[1,:]
# q2=A[2,:]
# p1=A[3,:]
# p2=A[4,:];

# Q1=(q1+q2)/sqrt(2)
# Q2=(q1-q2)/sqrt(2)
# P1=(p1+p2)/sqrt(2)
# P2=(p1-p2)/sqrt(2);

# figure()
# plot(q1,p1,"b",linewidth=.1)
# plot(q2,p2,"r",linewidth=.1);

# figure()
# plot(Q1,P1,"b",linewidth=.1)
# plot(Q2,P2,"r",linewidth=.1);



t_end=1e5;
# figure()
# @time Q_PSS,P_PSS=PSS_function( q0[2], p0[2], H, t_end,count);
# plot(Q_PSS,P_PSS,".", markersize=4, c
#     ="k");

N_iter_Q=200;#50
Q_start=-2
Q_end=2
N_iter_P=200
P_start=-1.4
P_end=1.4
ArrP=range(P_start,stop=P_end,length=N_iter_P)
ArrQ=range(Q_start,stop=Q_end,length=N_iter_Q)
Brandons_Colors=["#393b79" ,"#5254a3","#6b6ecf","#9c9ede" ,"#637939","#8ca252" ,"#b5cf6b" ,"#cedb9c" ,"#8c6d31","#bd9e39" ,"#e7ba52","#e7cb94","#843c39","#ad494a" ,"#d6616b","#e7969c" ,"#7b4173" ,"#a55194","#ce6dbd" ,"#de9ed6"];

# figure()
@showprogress for k=1:N_iter_Q
    Q_PSS,P_PSS=PSS_function(ArrQ[k], 0, Energy, t_end);       
    current_color=Brandons_Colors[mod(k,length(Brandons_Colors))+1]
    if Q_PSS!=0
        mat"plot($Q_PSS,$P_PSS,'.','MarkerSize',1,'color',$current_color); hold on;"
        # plot(Q_PSS,P_PSS,".", markersize=1,c=current_color);
    end
end
@showprogress for k=1:N_iter_P
    Q_PSS,P_PSS=PSS_function(0,ArrP[k], Energy, t_end);       
    current_color=Brandons_Colors[mod(k,length(Brandons_Colors))+1]
    if Q_PSS!=0
        mat"plot($Q_PSS,$P_PSS,'.','MarkerSize',1,'color',$current_color); hold on;"
        #  plot(Q_PSS,P_PSS,".", markersize=1,c=current_color);
    end
end

mat"savefig($file_name)"

end
