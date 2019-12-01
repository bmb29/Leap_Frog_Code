using ProgressMeter
using Printf
using BSON: @save, @load
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

location="/media/brandon_behring/Extra_Space/BSON_FILES/"
file_name="BandF_Dimer_lagrangian_descriptors_data0_129250000000000_ 1000_2_2.250_4.5002019_12_01T12_53_22_675"

@load location*file_name*".bson" SAVE_DATA

gradM=SAVE_DATA["gradM"]
Energy=SAVE_DATA["Energy"]
LD=SAVE_DATA["LD"]
ArrP=SAVE_DATA["ArrP"]
ArrQ=SAVE_DATA["ArrQ"]
N=SAVE_DATA["N"]

ArrP=ArrP[3:N-2]
ArrQ=ArrQ[3:N-2]

figure()
h_title = @sprintf("h= %.6f",Energy)
plot(ArrQ, LD[3,3:N-2], c=:blue)
plot(ArrQ,gradM[1,:],c=:red)
# figure()
# pcolormesh(ArrP,ArrQ,log1p.(gradM))
# colorbar()


# xlabel("q2")
# ylabel("q2")
# title(h_title)
# xlim(-2.25,2.25) 
# ylim(-4.5,4.5) 