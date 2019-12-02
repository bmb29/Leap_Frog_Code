using ProgressMeter
using Printf
using MATLAB
using BSON: @save, @load

include("Aref_LD_grad_function.jl")


@everywhere begin
    include("Aref_LD_grad_function.jl")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_Code")
    using ProgressMeter
    using Printf
    using MATLAB
    using .Aref_LD_grad_function
    using BSON: @save, @load

    Energy=.129
    t_typical=Aref_LD_grad_function.final_T(Energy)
    N=1000;#50

    n_iter_Q=N;#50
    # Q_end=.05
    Q_end=.2
    # Q_end=-Q_start
    

    Q_start=0.0
    n_iter_P=N

    P_end=6.5*1e-3
    P_start=-P_end
    t_end=35*t_typical
    N=n_iter_Q*n_iter_P
    t_end_mesh = t_end * ones(n_iter_Q,n_iter_P)
    Energy_mesh = Energy* ones(n_iter_Q,n_iter_P)
    ArrP=range(P_start,stop=P_end,length=n_iter_P)
    ArrQ=range(Q_start,stop=Q_end,length=n_iter_Q)
    mesh = [(P, Q) for Q in ArrQ, P in ArrP]
    mesh_list = reshape(mesh, 1, :)
end
@everywhere location="/media/brandon_behring/Extra_Space/MATLAB_FIGURES/"
@everywhere h=replace(@sprintf("%.15f",Energy),"."=>"_")
@everywhere h=replace(@sprintf("%.15f",Energy),"."=>"_")
@everywhere nQ = @sprintf("_%d",n_iter_Q)
@everywhere tend= @sprintf("_%d",t_end)
@everywhere tend_T= @sprintf("%4f",t_end)
@everywhere h_BSON="lagrangian_descriptors_data"*h*nQ*tend*".bson"
@everywhere file_name=location*"MATLAB_LD_"*h*nQ*tend*".fig"
@everywhere h_title = @sprintf("h= %.6f",Energy)*"with tend="*tend_T


LD= @showprogress pmap(Aref_LD_grad_function.LD_Helper,mesh, Energy_mesh, t_end_mesh)
# LD= @showprogress pmap(Aref_LD_grad_function.gradM,mesh, Energy_mesh, t_end_mesh)

SAVE_DATA=Dict("LD"=>LD, "mesh"=>mesh, "ArrP"=>ArrP, "ArrQ"=>ArrQ, "Energy"=>Energy)
@save h_BSON SAVE_DATA

dx=(P_end-P_start)/N
dy=(Q_end-Q_start)/N

gradM_X,gradM_Y=Aref_LD_grad_function.gradient_matrix(LD,dx,dy)

gradM_X,gradM_Y=gradient_matrix(LD,dx,dy)

Log_gradM_X=log10.(gradM_X)
Log_gradM_Y=log10.(gradM_Y)

biggest_X=maximum(filter(!isnan,Log_gradM_X))
biggest_Y=maximum(filter(!isnan,Log_gradM_Y))
smallest_X=minimum(filter(!isnan,Log_gradM_X))
smallest_Y=minimum(filter(!isnan,Log_gradM_Y))



# gradMl=Aref_LD_grad_function.gradient_matrix_4(LD,dx,dy)

# gradMl2=Aref_LD_grad_function.gradient_matrix_4(gradMl,dx,dy)
# gradMl=log10.(gradMl2)


Nl,Nl=size(gradM_X)
gradM=zeros(Nl,Nl)
for i=1:Nl
    for j=1:Nl
        cutoff_X=biggest_X-2
        cutoff_Y=biggest_Y-2
        X=maximum([Log_gradM_X[i,j], cutoff_X])
        Y=maximum([Log_gradM_Y[i,j], cutoff_Y])
        if Log_gradM_X[i,j]>cutoff_X || Log_gradM_Y[i,j]>cutoff_Y
            Z=1
        else
            Z=0
        end
        # X=log10( gradM_X[i,j])
        # Y=log10( gradM_Y[i,j])
        # Z=maximum([X,Y])
        # Z=log10(10^(2*X)+10^(2*Y))
        gradM[i,j]=Z
        # if Z>cutoff
        #     gradM[i,j]=0
        # elseif Z<2.0
        #     gradM[i,j]=1
        # else
        #     gradM[i,j]=2
        # end
    end
end




Nl,Nl=size(gradM_X)
gradM2=zeros(Nl,Nl)
for i=1:Nl
    for j=1:Nl
        cutoff_X=biggest_X-4
        cutoff_Y=biggest_Y-4
        X=maximum([Log_gradM_X[i,j], cutoff_X])
        Y=maximum([Log_gradM_Y[i,j], cutoff_Y])

        if Log_gradM_X[i,j]>cutoff_X || Log_gradM_Y[i,j]>cutoff_Y
            Z=1
        else
            Z=0
        end
        # X=log10( gradM_X[i,j])
        # Y=log10( gradM_Y[i,j])
        # Z=maximum([X,Y])
        Z=log10(10^(2*X)+10^(2*Y))
        gradM2[i,j]=Z
        # if Z>cutoff
        #     gradM[i,j]=0
        # elseif Z<2.0
        #     gradM[i,j]=1
        # else
        #     gradM[i,j]=2
        # end
    end
end




mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"
# mat"imagesc([$P_start,$P_end ],[$Q_start,$Q_end],$LD)"

mat"imagesc([$P_start,$P_end ],[0,$Q_end],$LD)"
mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$LD)"
# mat"imagesc([0,$P_end],[0,$Q_end ],$LD)"
# mat"imagesc([0,-$P_end],[0,$Q_end ],$LD)"
# mat"imagesc([0,$P_end],[0,-$Q_end ],$LD)"
# mat"imagesc([0,-$P_end],[0,-$Q_end ],$LD)"


mat"colorbar"
mat"axis([-$P_end,$P_end,-$Q_end,$Q_end ])"


mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"
mat"imagesc([$P_start,$P_end ],[$Q_start,$Q_end],$gradM)"
mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$gradM)"

# mat"imagesc([$P_start,$P_end ],[0,$Q_end],$gradM)"
# mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$gradM)"
# mat"imagesc([0,$P_end],[0,$Q_end ],$gradM)"
# mat"imagesc([0,-$P_end],[0,$Q_end ],$gradM)"
# mat"imagesc([0,$P_end],[0,-$Q_end ],$gradM)"
# mat"imagesc([0,-$P_end],[0,-$Q_end ],$gradM)"

mat"colorbar"
mat"axis([-$P_end,$P_end,-$Q_end,$Q_end ])"

mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
mat"title($h_title)"
mat"imagesc([$P_start,$P_end ],[0,$Q_end],$gradM2)"
mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$gradM2)"

# mat"imagesc([$P_start,$P_end ],[0,$Q_end],$gradMl)"
# mat"imagesc([$P_end,$P_start ],[0,-$Q_end],$gradMl)"
# mat"imagesc([0,$P_end],[0,$Q_end ],$gradMl)"
# mat"imagesc([0,-$P_end],[0,$Q_end ],$gradMl)"
# mat"imagesc([0,$P_end],[0,-$Q_end ],$gradMl)"
# mat"imagesc([0,-$P_end],[0,-$Q_end ],$gradMl)"

mat"colorbar"
mat"axis([-$P_end,$P_end,-$Q_end,$Q_end ])"

mat"savefig($file_name)"


