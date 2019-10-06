
using Printf
using Formatting
using DifferentialEquations
using ProgressMeter
include("Yfind.jl")
include("is_it.jl")
using PyCall
using PyPlot
PyPlot.svg(true)


function condition(u,t,integrator) # Event when event_f(u,t) == 0
 u[1]
end
function affect!(integrator)
end
# only keep plane X=0
cb2 = ContinuousCallback(condition,affect!,nothing)



cb=CallbackSet(cb2)
Energy=.45
Energy=0.125789
Energy=BigFloat("0.125789")


H=(2*Energy)^2
n_iter=10
t_end=300
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
    return
end
# refernce
# z=X+iP W=Q+iY . Plane Q=0, P=0
# z=u[1]+i u[2] ; w=u[3]+i u[4]

#Q=0, have P find Y

ArrQ=range(BigFloat("1e-5"),stop=BigFloat(".001"),length=n_iter)
ArrP=range(BigFloat("1e-5"),stop=BigFloat(".01"),length=n_iter)
#initialize stuff



Last_Q=zeros(BigFloat ,0)
Last_P=zeros(BigFloat ,0)
Last_X=zeros(BigFloat ,0)
Last_Y=zeros(BigFloat ,0)
Last_T=zeros(BigFloat ,0)
MAX=zeros(BigFloat ,0)

#all colors but red, reserve for exit points
COLOR=["b", "g", "c", "m", "y", "k", "w"]
# axis([ -.3,.3,-.03, .03]);# Efigure(dpi=3600)
axis([ -.025, .025,-.25, .25]) #around E=0.1257895
# axis([ -5, 5, -2, 2])
@showprogress 1 "Computing..." for j=1:n_iter
            Q=ArrQ[j]
            for k=1:n_iter

                    P=ArrP[k]
                    Y=Yfind(Q,P,H)
                    if ~isempty(Y) && Y>0
                        u0=zeros(4)
                        u0[1]=0 #X
                        u0[2]=P #P
                        u0[3]=Q #Q
                        u0[4]=Y #Y
                        #are the initial conditions right?
                        u0=convert(Array{Float64},u0)
                        H=convert(Float64,H)
                        prob = ODEProblem(Eq_of_M,u0,(0., t_end))
                        sol=solve(prob, Vern9(),maxiters=1e10, reltol=1e-13,abstol=1e-15,callback=cb,save_start=true,save_end=true,save_everystep=false)
                        # plot(A[3,:],-A[2,:],color=current_color,",")
                        #testing if energy is conserved at end?
                        A=sol[:,:]
                        #values at end to check if energy is coservered
                        dH=zeros(length(A[1,:]))
                        umax=zeros(length(A[1,:]))
                        for i=1:length(A[1,:])
                            uf=zeros(4)
                            uf[1]=A[1,i] #X
                            uf[2]=A[2,i] #P
                            uf[3]=A[3,i]#Q
                            uf[4]=A[4,i] #Y
                            dH[i]=abs(H_test(uf)-H)
                            umax[i]=maximum(uf)
                        end
                        dH_max=maximum(dH)
                        #filter if energy is not conservered-toss
                        if dH_max<1e-8# t,A=solve(prob, RK4(),maxiters=1e20, reltol=1e-6,abstol=1e-8,callback=cb,save_start=true,save_end=true,save_everystep=false)
                            current_color=COLOR[mod(k,length(COLOR))+1]

                            #convert to float64 before plotting
                            A=convert(Array{Float64},A)
                            plot(A[2,2:end-1],A[3,2:end-1],color=current_color,",")
                            plot(-A[2,2:end-1],A[3,2:end-1],color=current_color,",")
                            plot(A[2,2:end-1],-A[3,2:end-1],color=current_color,",")
                            plot(-A[2,2:end-1],-A[3,2:end-1],color=current_color,",")


                            how_far=maximum(umax)
                            push!(MAX,how_far)
                            # add final values to list- I had to let the callback save the final value so
                            # I need to check if the final value is actually on the section X=0
                            # println(abs(A[1,end]))
                            # println(abs(A[1,end-1]))
                            if abs(A[1,end])<1e-7
                                push!(Last_Q,A[3,end])
                                push!(Last_P,A[2,end])
                                push!(Last_X,A[1,end])
                                push!(Last_Y,A[4,end])
                                push!(Last_T,sol.t[end])
                            elseif abs(A[1,end-1])<1e-7
                                push!(Last_Q,A[3,end-1])
                                push!(Last_P,A[2,end-1])
                                push!(Last_X,A[1,end-1])
                                push!(Last_Y,A[4,end-1])
                                push!(Last_T,sol.t[end-1])
                            end
                        end
                end
            end
        end

H_vec=zeros(length(Last_P))
u_max=zeros(length(Last_P))
for j=1:length(Last_P)
    uf=zeros(BigFloat ,4)
    uf[1]=Last_X[j]
    uf[2]=Last_P[j]
    uf[3]=Last_Q[j]
    uf[4]=Last_Y[j]
    u_max[j]=maximum(abs.([uf[1],uf[2],uf[3],uf[4]]))
    H_vec[j]=H_test(uf)/H
end

H_vec=convert(Array{Float64},H_vec)
Last_P=convert(Array{Float64},Last_P)
Last_Q=convert(Array{Float64},Last_Q)
plot(Last_P,Last_Q,color="r",".",markersize=.5, markeredgewidth=.1)
plot(Last_P,-Last_Q,color="r",".",markersize=.5,  markeredgewidth=.1)
plot(-Last_P,Last_Q,color="r",".",markersize=.5,  markeredgewidth=.1)
plot(-Last_P,-Last_Q,color="r",".",markersize=.5, markeredgewidth=.1)
xlabel("P")
ylabel("Q")
savefig("work.svg")
# figure()
# title("Is energy conserved?")
# plot(1:length(Last_P)-1,log.(H_vec[1:end-1]))

# figure()
# plot(1:length(Last_P),sort(Last_T) )
# figure()
# # plot(sort(MAX))
# titleS="P-Q Plane PSS where X=0 at Energy "*string(Energy)*" with "*string(n_iter^2)*" iterations."
# filenameS="PSS_Red_ExitE_"*string(Energy)*"_n_iter_"*string(n_iter)*"_t_"*string(t_end)*"withprojection.png"
