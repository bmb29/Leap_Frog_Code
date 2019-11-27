using DifferentialEquations
using ProgressMeter
using Printf
using MATLAB
include("escape_exit_function.jl")

max_hit=5
t_end=1e3
n_iter_P=1000
n_iter_Q=2001
width=.7
height=1
ArrP=range(-width,stop=width,length=n_iter_P)
ArrQ=range(-height,stop=height,length=n_iter_Q)
location="/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
H=range(.21,stop=0.24,length=21)

@showprogress 1 "Computing..." for Energy in H

println(Energy)

Q_0=zeros(0)
P_0=zeros(0)

Q_1=zeros(0)
P_1=zeros(0)

Q_2=zeros(0)
P_2=zeros(0)

Q_3=zeros(0)
P_3=zeros(0)

Q_4=zeros(0)
P_4=zeros(0)

Q_5=zeros(0)
P_5=zeros(0)
#
# @showprogress 1 "Computing..." for Energy in H

h=replace(@sprintf("%.13f",Energy),"."=>"_")
file_name=location*"Escape_"*h*".fig"
for j=1:n_iter_P
    for k=1:n_iter_Q
        K=escape_exit_function(ArrQ[k],ArrP[j], Energy, t_end, max_hit)
        if K==0
            push!(Q_0,ArrQ[k])
            push!(P_0,ArrP[j])
        elseif K==1
            push!(Q_1,ArrQ[k])
            push!(P_1,ArrP[j])
        elseif K==2
            push!(Q_2,ArrQ[k])
            push!(P_2,ArrP[j])
        elseif K==3
            push!(Q_3,ArrQ[k])
            push!(P_3,ArrP[j])
        elseif K==4
            push!(Q_4,ArrQ[k])
            push!(P_4,ArrP[j])
        elseif K==5
            push!(Q_5,ArrQ[k])
            push!(P_5,ArrP[j])
        # elseif K==-1
        #     push!(Q_bound,ArrQ[k])
        #     push!(P_bound,ArrP[j])
        end
    end
end


size=3
#
# mat"axis([ -$width,$width,-$height,$height ])"
mat"figure(); hold on;"
mat"axis([ -$width,$width,-$height,$height ])"
# mat"plot($P_bound,$Q_bound ,'k.','MarkerSize',10)"
mat"plot($P_0,$Q_0 ,'b.','MarkerSize',$size)"
mat"plot($P_1,$Q_1 ,'r.','MarkerSize',$size)"
mat"plot($P_2,$Q_2 ,'g.','MarkerSize',$size)"
mat"plot($P_4,$Q_4 ,'c.','MarkerSize',$size)"
mat"plot($P_5,$Q_5 ,'y.','MarkerSize',$size)"
mat"savefig($file_name)"
mat"close"
 end
