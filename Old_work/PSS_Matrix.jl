
using Printf
using Formatting
using DifferentialEquations
using ProgressMeter
using JLD2, FileIO
using HDF5

include("Yfind.jl")
include("is_it.jl")
using PyCall
pygui(:qt)
using PyPlot
pygui(true)

tol_dist=1e-2
max_hit=100


condition(u,t,integrator)= u[5]>max_hit
function condition2(u,t,integrator) # Event when event_f(u,t) == 0
 u[1]
end
affect!(integrator) = terminate!(integrator)
function affect2!(integrator)
integrator.u[5]=integrator.u[5]+1
 end
cb2 = ContinuousCallback(condition2,affect2!,nothing,rootfind = true,   save_positions = (false,true))
cb1 = DiscreteCallback(condition,affect!)

cb=CallbackSet(cb2,cb1)


Energy=0.25
Energy=.1251
# Energy=0.1257895
H=(2*Energy)^2
n_iter=10
t_end=10000

# H1(p,q,t)= ( (q[2]-q[1])^2+(p[2]-p[1])^2 )*( (q[2]+q[1])^2+(p[2]+a[1])^2 )
# H2(p,q,t)= (p[2]^4+2*p[2]^2*(q[1]^2-1)+(1+q[1]^2)^2 )*(q[2]^4+2*q[2]^2*(p[1]^2+1)+(p[1]^2-1)^2 )
# HAM(p,q,t)=H1(p,q,t)/H2(p,q,t)
H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))


ODE1(z,w)=conj(  im * w.*( 1 ./(z.^2-w.^2)+1 ./(1+w.^2) ))
ODE2(z,w)=conj(  im * z.*( 1 ./(w.^2-z.^2)+1 ./(1+z.^2) ))
function Eq_of_M(du,u,p,t)
    du[1]=real(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[2]=imag(ODE1(u[1]+im*u[2],u[3]+im*u[4]))
    du[3]=real(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[4]=imag(ODE2(u[1]+im*u[2],u[3]+im*u[4]))
    du[5]=0
    return
end
# refernce
# z=X+iP W=Q+iY . Plane Q=0, P=0
# z=u[1]+i u[2] ; w=u[3]+i u[4]

#Q=0, have P find Y

ArrQ=range(1e-5,stop=.1,length=n_iter)
ArrP=range(1e-5,stop=.02,length=n_iter)
#
# ArrQ=range(1e-5,stop=2,length=n_iter)
# ArrP=range(1e-5,stop=.1,length=n_iter)
#initialize stuff
#initialize stuff

Last_Q=zeros(0)
Last_P=zeros(0)

Tabs_X=zeros((max_hit+2),n_iter^2)
Tabs_Y=zeros((max_hit+2),n_iter^2)
Tabs_Q=zeros((max_hit+2),n_iter^2)
Tabs_P=zeros((max_hit+2),n_iter^2)

#all colors but red, reserve for exit points
COLOR=["b", "g", "c", "m", "y", "k", "w"]
# figure()

# axis([ -.3,.3,-.03, .03]);# Efigure(dpi=3600)
axis([ -.025, .025,-.25, .25]) #around E=0.1257895
# axis([ -1.2, 1.2, -1.5, 1.5])
@showprogress 1 "Computing..."    for j=1:n_iter
        Q=ArrQ[j]
            for k=1:n_iter
                    P=ArrP[k]
                    Y=Yfind(Q,P,H)
                    if ~isempty(Y) && Y>0

                        u0=zeros(5)
                        u0[1]=0 #X
                        u0[2]=P #P
                        u0[3]=Q #Q
                        u0[4]=Y #Y
                        u0[5]=1 # counter variable (starts on PSS)
                        #are the initial conditions right?


                        prob = ODEProblem(Eq_of_M,u0,(0., t_end))
                        sol=solve(prob, Vern9(),maxiters=1e20, reltol=1e-13,abstol=1e-15,callback=cb,save_start=true,save_end=true,save_everystep=false)
                        # sol=solve(prob, RK4(),maxiters=1e20,reltol=1e-6,abstol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)
                        # plot(A[3,:],-A[2,:],color=current_color,",")
                        #testing if energy is conserved at end?
                       A =sol[:,:]
                        #values at end to check if energy is coservered

                        #save hit points
                        for i=1:length(A[1,:])-1
                            Tabs_P[i,(j-1)*k+k]=A[2,i]
                            Tabs_Q[i,(j-1)*k+k]=A[3,i]
                            Tabs_X[i,(j-1)*k+k]=A[1,i]
                            Tabs_Y[i,(j-1)*k+k]=A[4,i]
                        end
                        dH=zeros(length(A[1,:]))
                        umax=zeros(length(A[1,:]))
                        uf=zeros(4)
                        for i=1:length(A[1,:])
                        uf[1]=A[1,i] #X
                        uf[2]=A[2,i] #P
                        uf[3]=A[3,i]#Q
                        uf[4]=A[4,i] #Y
                        dH[i]=abs(H_test(uf)-H)
                        umax[i]=maximum(uf)
                        end
                        dH_max=maximum(dH)
                        #filter if energy is not conservered-toss
                        if dH_max<1e-5# t,A=solve(prob, RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)
                            current_color=COLOR[mod(k,length(COLOR))+1]
                            plot(A[2,2:end-1],A[3,2:end-1],color=current_color,".",markersize=.1, markeredgewidth=.1)
                            plot(-A[2,2:end-1],A[3,2:end-1],color=current_color,".",markersize=.1, markeredgewidth=.1)
                            plot(A[2,2:end-1],-A[3,2:end-1],color=current_color,".",markersize=.1, markeredgewidth=.1)
                            plot(-A[2,2:end-1],-A[3,2:end-1],color=current_color,".",markersize=.1, markeredgewidth=.1)
                            #testing if energy is conserved at end?

                            # add final values to list- I had to let the callback save the final value so
                            # I need to check if the final value is actually on the section X=0
                            # println(abs(A[1,end]))
                            # println(abs(A[1,end-1]))
                            uf=zeros(5)
                            uf[1]=A[1,end] #X
                            uf[2]=A[2,end] #P
                            uf[3]=A[3,end]#Q
                            uf[4]=A[4,end] #Y
                            uf[5]=0
                            if is_it_exit(uf,H,tol_dist)
                                push!(Last_Q,A[3,end-1])
                                push!(Last_P,A[2,end-1])

                            end
                        end
                    end
                end
            end


# titleS=
"P-Q Plane  PSS where X=0 at Energy "*string(Energy)*" with "*string(n_iter^2)*" iterations."
# filenameS="PSS_Red_ExitE_"*string(Energy)*"_n_iter_"*string(n_iter)*"_t_"*string(t_end)*"withprojection.png"
plot(Last_P,Last_Q,color="r",".",markersize=.2, markeredgewidth=.1)
plot(Last_P,-Last_Q,color="r",".",markersize=.2,  markeredgewidth=.1)
plot(-Last_P,Last_Q,color="r",".",markersize=.2,  markeredgewidth=.1)
plot(-Last_P,-Last_Q,color="r",".",markersize=.2, markeredgewidth=.1)
xlabel("P")
ylabel("Q")

# title(titleS)
# savefig(filenameS)
#
# H_vec=zeros(length(Last_P))
# u_max=zeros(length(Last_P))
# for j=1:length(Last_P)
#     uf[1]=Last_X[j]
#     uf[2]=Last_P[j]
#     uf[3]=Last_Q[j]
#     uf[4]=Last_Y[j]
#     u_max[j]=maximum(abs.([uf[1],uf[2],uf[3],uf[4]]))
#     H_vec[j]=H_test(uf)/H
# end
#
# figure()
# title("Is energy conserved?")
# plot(1:length(Last_P)-1,log.(H_vec[1:end-1]))
# figure()
# plot(1:length(Last_P),sort(Last_T) )
# figure()
# plot(sort(MAX))
