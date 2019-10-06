include("Yfind.jl")

# setprecision(100)
# Q=.2
# P=.3
#
#
# Energy=.25
# H=(2*Energy)^2
#
# Y=Yfind(Q,P,H)
#
# u0=zeros(4)
#
# u0[1]=0 #X
# u0[2]=P #P
# u0[3]=Q #Q
# u0[4]=Y #Y
#
# TEST1_0=H_test(u0)-H

Q=BigFloat(".2")
P=BigFloat(".3")


Energy=BigFloat(".25")
H=(2*Energy)^2

Y=Yfind(Q,P,H)

u0=zeros(BigFloat ,4)

u0[1]=0 #X
u0[2]=P #P
u0[3]=Q #Q
u0[4]=Y #Y

TEST1_0=H_test(u0)-H
