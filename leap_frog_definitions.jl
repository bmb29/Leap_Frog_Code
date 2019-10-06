
Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+(1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
Centroid_ODE(z,w)=conj( 1 ./(1+z.^2)+1 ./(1+w.^2) )
Ham_numerator(Q,P)=((Q[1]+Q[2])^2+(P[1]+P[2])^2)*((Q[1]-Q[2])^2+(P[1]-P[2])^2)
Ham_denominator(Q,P)=(Q[1]^2+(P[2]-1)^2)*(Q[1]^2+(P[2]+1)^2)*(Q[2]^2+(P[1]-1)^2)*(Q[2]^2+(P[1]+1)^2)
Aref_Hamiltonian(Q,P,t)=1/2*sqrt(Ham_numerator(Q,P)/Ham_denominator(Q,P));
Hamiltonian_Dimer(q,p,t)=4sqrt(abs( (q[1]^2+p[1]^2)*(q[2]^2+p[2]^2) )/
(
        (2-2sqrt(2)*(p[1]+p[2])+(p[1]+p[2])^2+(q[1]-q[2])^2 )*
        (2+2sqrt(2)*(p[1]+p[2])+(p[1]+p[2])^2+(q[1]-q[2])^2 )*
        (2-2sqrt(2)*(p[1]-p[2])+(p[1]-p[2])^2+(q[1]+q[2])^2 )*
        (2+2sqrt(2)*(p[1]-p[2])+(p[1]-p[2])^2+(q[1]+q[2])^2 )

));
function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
    return
end
function Eq_of_M_LAB_FRAME(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=real(Centroid_ODE(u[1]+im*u[2],u[3]+im*u[4]))
    du[6]=imag(Centroid_ODE(u[1]+im*u[2],u[3]+im*u[4]))
    return
end;


function Yfind(Q,P,H)
    try
        Y=sqrt( (1+(-1)*H*((-1)+P^2)^2)^(-1)*(P^2+(-1)*Q^2+H*((-1)+P^2)
        ^2*((-1)+Q^2)+sqrt(H*((-1)+P^2)^4+(-4)*((-1)+H*((-1)+P^2)^2)
        *((-1)*P^2+H*((-1)+P^2)^2)*Q^2 ))   )
     catch
        Y=zeros(0)
    end
end
function Yfind_Aref(Q,P,H_Aref)
    H=(2*H_Aref)^2
    try
        Y=sqrt( (1+(-1)*H*((-1)+P^2)^2)^(-1)*(P^2+(-1)*Q^2+H*((-1)+P^2)
        ^2*((-1)+Q^2)+sqrt(H*((-1)+P^2)^4+(-4)*((-1)+H*((-1)+P^2)^2)
        *((-1)*P^2+H*((-1)+P^2)^2)*Q^2 ))   )
     catch
        Y=zeros(0)
    end
end
using Roots
function Q1_find_dimer(Q2,P2,H_dimer)
    Q1_to_find(Q1)=Hamiltonian_Dimer([Q1,Q2],[0,P2],1)-H_dimer
    try
        Q1=find_zeros(Q1_to_find,0,4,maxeval=100,maxfnevals=300,tol=1e-15)
     catch
        Q1=zeros(0)
    end
end
function P1_find_dimer(Q2,P2,H_dimer)
    P1_to_find(P1)=Hamiltonian_Dimer([0,Q2],[P1,P2],1)-H_dimer
    try
        P1=find_zeros(P1_to_find,0,4,maxeval=10000,maxfnevals=300,tol=1e-15)
     catch
        P1=zeros(0)
    end
end
function P1_find_dimer_second(Q2,P2,H_dimer)
    P1_to_find(P1)=Hamiltonian_Dimer([0,Q2],[P1,P2],1)-H_dimer
    try
        P1=find_zero(P1_to_find,1.4,maxeval=10000,maxfnevals=300,tol=1e-15)
     catch
        P1=zeros(0)
    end
end
