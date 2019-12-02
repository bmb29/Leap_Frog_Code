using ProgressMeter
using Printf
using BSON: @save, @load
# using PyCall
# pygui(:qt)
# using PyPlot
# pygui(true)
using MATLAB
include("Dimer_Lagrangian_Descriptor.jl")


location="/media/brandon_behring/Extra_Space/BSON_FILES/"
file_name="BandF_Dimer_lagrangian_descriptors_data0_129250000000000_ 1000_2_2.250_4.5002019_12_01T12_53_22_675"



N=Int(1e5)
Energy=.129
t_end=5
Q_start=.2
Q_end=1.2

Q=range(Q_start, stop=Q_end, length=N)
mesh=[(q, 0.0) for q in Q]

LD=@showprogress map(Dimer_Lagrangian_Descriptor.Dimer_Lagrangian_Descriptor_Function, mesh, Energy*ones(N), t_end*ones(N))
dx=(Q_end-Q_start)/N
gradM=Dimer_Lagrangian_Descriptor.D_X(LD,dx)
gradMl=log10.(gradM)
# d2M=log10.(Dimer_Lagrangian_Descriptor.D_X(gradM,dx))
Q1=Q[4:N-3]
Q2=Q[7:N-6]

mat"figure();plot($Q, $LD,'k');hold on;"
mat"plot($Q1, $gradMl,'r')"
# mat"plot($Q2, $d2M,'b')"

# @load location*file_name*".bson" SAVE_DATA

# gradM=SAVE_DATA["gradM"]
# Energy=SAVE_DATA["Energy"]
# LD=SAVE_DATA["LD"]
# ArrP=SAVE_DATA["ArrP"]
# ArrQ=SAVE_DATA["ArrQ"]
# N=SAVE_DATA["N"]
# biggest=maximum(filter(!isnan,gradM))
# smallest=minimum(filter(!isnan,gradM))


# ArrP=ArrP[4:N-3]
# ArrQ=ArrQ[4:N-3]
# N=N-4
# # figure()
# # h_title = @sprintf("h= %.6f",Energy)
# # plot(ArrQ, LD[3,3:N-2], c=:blue)
# plot(ArrQ,gradM[1,:],c=:red)

# mat" imagesc($gradM) "
# grad=zeros(N,N)
# for i=1:N
#     for j=1:N
#        M=gradM[i,j]
#        if M<15.0
#             grad[i,j]=0.0
#        else
#             grad[i,j]=1.0
#        end
#     end
# end
# mat"figure"
# mat" imagesc($grad) "
        