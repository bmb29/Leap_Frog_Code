using DifferentialEquations
using MATLAB
using PyCall
pygui(:qt)
using PyPlot
pygui(true)
include("leap_frog_definitions.jl")
include("PSS_Definitions_Dimer_X.jl")
include("orbit_plots_leapfrog_dimer.jl")


H=.075
h=1/(2*H)
max_time=20
a=(2+h-2*sqrt(h+1))/h
Y=Yfind_Aref(1e-3,0,H);
Yfind(h)=sqrt(h/(2h+1));
P=Yfind(H)


q10=0
# q2=0
# p2=P

# q20=.3587
q20=.16
# q20=0.0
p20=0

p10=P1_poly(q20,p20,H)

Q10=(q10+q20)/sqrt(2)
Q20=(q10-q20)/sqrt(2)
P10=(p10[1]+p20)/sqrt(2)
P20=(p10[1]-p20)/sqrt(2);
println(p10)

u0=[Q10;P20;Q20;P10;0;0] # X P Q Y
P0=[Q10,Q20]; Q0=[P10,P20];
tspan = (0.0,max_time);

println(Aref_Hamiltonian(Q0,P0,0))
prob = ODEProblem(Eq_of_M_LAB_FRAME,u0,tspan)
t,A=solve(prob,Tsit5(),reltol=1e-14,abstol=1e-14,maxiters=1e15);


N=length(t.t)
Z=A[1,:]+im *A[2,:];
W=A[3,:]+im *A[4,:];

Q1=A[1,:]
Q2=A[3,:]
P1=A[4,:]
P2=A[2,:];

# figure()
# plot(Q1,P1,"b",linewidth=.1)
# plot(Q2,P2,"r",linewidth=.1)
# title("Aref coordinates");

q1=(Q1+Q2)/sqrt(2)
q2=(Q1-Q2)/sqrt(2)
p1=(P1+P2)/sqrt(2)
p2=(P1-P2)/sqrt(2)

figure()
plot(q1,p1,"b",linewidth=.1)
plot(q2,p2,"r",linewidth=.1)
title("Dimer coordinates");

Centeroid=A[5,:]+im *A[6,:];

lin_impulse=im*(1+a )*ones(N);
z_1pos=.5*(Centeroid+Z+lin_impulse+W);
z_1neg=.5*(Centeroid+Z-lin_impulse-W);
z_2pos=.5*(Centeroid-Z+lin_impulse-W);
z_2neg=.5*(Centeroid-Z-lin_impulse+W);

z_1posX=real(z_1pos)
z_1posY=imag(z_1pos)

z_2posX=real(z_2pos)
z_2posY=imag(z_2pos)

z_1negX=real(z_1neg)
z_1negY=imag(z_1neg)

z_2negX=real(z_2neg)
z_2negY=imag(z_2neg)

figure()
plot(z_1posX,z_1posY,c="b",linewidth=1)
plot(z_1negX,z_1negY,c="r",linewidth=1)
plot(z_2posX,z_2posY,c="b",linewidth=.5,linestyle="--")
plot(z_2negX,z_2negY,c="r",linewidth=.5,linestyle="--")
ylim(-3.2,3.2)
axis("equal")
title("Lab coordinates");

Q_PSS,P_PSS=PSS_function(q20,p20, H,  1e5)
mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"xlabel('Q')"
mat"ylabel('P')"
Q_fix=sqrt(6)/3

mat"axis([-3, 3, -3, 3])"
mat"plot($Q_fix,0,'b.','MarkerSize',30)"
mat"plot(-$Q_fix,0,'b.','MarkerSize',30)"
mat"plot(0,$P,'r.','MarkerSize',30)"
mat"plot(0,-$P,'r.','MarkerSize',30)"
mat"plot($Q_PSS,$P_PSS,'k.','MarkerSize',5); hold on;"
