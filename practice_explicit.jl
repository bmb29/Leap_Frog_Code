using DifferentialEquations
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

# u[1]=q1, u[2]=q2, u[3]=p1, u[4]=p2, u[5]=q1 hit, u[6]=p2 hit, u[7]=H