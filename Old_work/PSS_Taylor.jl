include("Yfind.jl")
include("is_it_taylor.jl")

H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))
ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
function Eq_of_M!(t,u,du)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    return nothing
end
function g(t,u,du)
    px_ = constant_term(u[4])
 # if px > 0...
 if px_ > zero(px_)
     # ...return x
     return u[1]
 else
     #otherwise, discard the crossing
     return nothing
 end
end

function PSS_function_taylor(Q,P, Energy,  t_end)
    Q_0=0
    P_0=0

    Q_1=0
    P_1=0


    Q_2=0
    P_2=0


    Q_3=0
    P_3=0


    tol_dist=1e-2
    H=(2*Energy)^2

    Tabs_Q=zeros(0)
    Tabs_P=zeros(0)
    Y=Yfind(Q,P,H)
    if ~isempty(Y) && Y>0
        u0=zeros(4)
        u0[1]=0 #X
        u0[2]=P #P
        u0[3]=Q #Q
        u0[4]=Y #Y
        # counter variable
        #are the initial conditions right?
        # prob = ODEProblem(Eq_of_M,u0,(0., t_end))
        tv_i, xv_i, tvS_i, xvS_i, gvS_i=taylorinteg(Eq_of_M!, g, u0, 0., t_end, 12, 1e-12,maxsteps=300000, newtoniter=30 )
        if is_it_taylor(xv_i,H,tol_dist)
            if isempty(xvS_i)
                Q_0=Q
                P_0=P
            else
                uf=zeros(5)
                uf[1]=xvS_i[end,1] #X
                uf[2]=xvS_i[end,2] #P
                uf[3]=xvS_i[end,3] #Q
                uf[4]=xvS_i[end,4] #Y
                dH=abs(H_test(uf)-H)
                if dH<1e-5
                    Tabs_Q=xvS_i[:,3]
                    Tabs_P=xvS_i[:,2]
                    if length(xvS_i[:1])<5
                        if length(Tabs_Q)==1
                            #never hits, exits
                            Q_0=Tabs_Q[end]
                            P_0=Tabs_P[end]
                        elseif length(Tabs_Q)==2
                            #hits once, exits
                            Q_0=Tabs_Q[end]
                            P_0=Tabs_P[end]
                            Q_1=Tabs_Q[end-1]
                            P_1=Tabs_P[end-1]
                        elseif length(Tabs_Q)==3
                            #hits twice, exits
                            Q_0=Tabs_Q[end]
                            P_0=Tabs_P[end]
                            Q_1=Tabs_Q[end-1]
                            P_1=Tabs_P[end-1]
                            Q_2=Tabs_Q[end-2]
                            P_2=Tabs_P[end-2]
                        elseif  length(Tabs_Q)==4
                            Q_0=Tabs_Q[end]
                            P_0=Tabs_P[end]
                            Q_1=Tabs_Q[end-1]
                            P_1=Tabs_P[end-1]
                            Q_2=Tabs_Q[end-2]
                            P_2=Tabs_P[end-2]
                            Q_3=Tabs_Q[end-3]
                            P_3=Tabs_P[end-3]
                        end
                    else
                    Q_0=Tabs_Q[end]
                    P_0=Tabs_P[end]
                    Q_1=Tabs_Q[end-1]
                    P_1=Tabs_P[end-1]
                    Q_2=Tabs_Q[end-2]
                    P_2=Tabs_P[end-2]
                    Q_3=Tabs_Q[end-3]
                    P_3=Tabs_P[end-3]
                    end
                end
            end
        end
    end
return     Q_0,P_0,Q_1,P_1,Q_2,P_2,Q_3,P_3
end
