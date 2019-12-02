using ProgressMeter
using Printf
using MATLAB
using BSON: @save, @load
using Dates
include("Dimer_Lagrangian_Descriptor.jl")


@everywhere begin
    include("Dimer_Lagrangian_Descriptor.jl")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_Code")
    using ProgressMeter
    using Printf
    using MATLAB
    using .Dimer_Lagrangian_Descriptor
    using BSON: @save, @load
    using Dates
    # H=range(.11,stop=0.145,length=107)
    # Energy=0.112012
    # Energy=.125789
    # Energy=.12501
    # Energy=.12501
    # Energy=.1249
    # Energy=.125

    N=Int(1e4)
    Energy=.129

    Q_start=.1
    Q_end=.5
    count = 1
    dp=1.0e-7
    Q=range(Q_start, stop=Q_end, length=N)
    P=range(-dp,stop=dp, length=5)
    P=0.0
    Nx=length(P)
    mesh=[(0.0,q) for q in Q]
    mesh = [(p, q) for q in Q, p in P]
    location_fig="/media/brandon_behring/Extra_Space/MATLAB_FIGURES/"
    right_now = replace(replace(replace(string(Dates.now()),"."=>"_"),":"=>"_"),"-"=>"_")
    t=1:50
end
while count <= length(t)
    @everywhere t_end=t[count]
    @everywhere h=replace(@sprintf("%.15f",Energy),"."=>"_")
    @everywhere Ns = @sprintf("_%5d",N)
    @everywhere tend= @sprintf("_%d",t_end)
    println(tend)
    @everywhere file_name=location_fig* "ld_grad_test_"*h*Ns*tend*right_now*".fig"
    @everywhere tend_T= @sprintf("%d",t_end)
    @everywhere h_title =h*" with tend="*tend_T*" and "*Ns
LD=@showprogress pmap(Dimer_Lagrangian_Descriptor.Dimer_Lagrangian_Descriptor_Function, mesh, Energy*ones(N,Nx), t_end*ones(N,Nx) )
dx=(Q_end-Q_start)/N
dy=2*dp/5
# gradM=Dimer_Lagrangian_Descriptor.gradient_matrix_4(LD,dy,dx)
# gradM2=Dimer_Lagrangian_Descriptor.D_X( LD[:,3],dx)
gradM2=Dimer_Lagrangian_Descriptor.D_X( LD,dx)
# gradMl=zeros(length(gradM))
gradMx=zeros(length(gradM2))
for i=1:length(gradM2)
    # gradMl[i]=log10(gradM[i])
    gradMx[i]=log10(gradM2[i])
end

# d2M=log10.(Dimer_Lagrangian_Descriptor.D_X(gradM,dx))
Q1=Q[3:N-2]
LDn=5*( LD/maximum(filter(!isnan,LD))-.7*ones(length(LD)) )
mat"figure();set(gcf, 'Position',  [0, 0, 2500, 1500]); hold on;"
mat"title($h_title)"
mat"plot($Q, $LDn,'b')"
# mat"plot($Q1, $gradMl,'r')"
mat"plot($Q1, $gradMx,'r')"
mat"savefig($file_name)"


@everywhere global count += 1


end