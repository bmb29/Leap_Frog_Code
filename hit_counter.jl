

# @everywhere module_dir = "/Users/brandonbehring/Desktop/Leap_Frog_2019"
# @everywhere push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
# # @everywhere module_dir ="/Users/brandonbehring/Desktop/Leap_Frog_2019"
# # @everywhere push!(LOAD_PATH, $module_dir)
# # @everywhere thisDir = dirname(@__FILE__())
# # @everywhere any(path -> path == thisDir, LOAD_PATH) || push!(LOAD_PATH, thisDir)
# include("escape.jl")
# @everywhere using .escape

using ProgressMeter
using Printf
using MATLAB
using BSON: @save, @load
using PolynomialRoots
include("Dimer_Eq_of_M.jl")

include("escape_dimer.jl")


@everywhere begin
    include("escape_dimer.jl")
    # push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
    push!(LOAD_PATH, "/home/brandon_behring/Desktop/Leap_Frog_Code")
    using ProgressMeter
    using Printf
    using MATLAB
    using .escape_dimer
    using BSON: @save, @load

    t_end = 1e5
    width = 2.5; height = 4.5
    n_iter_Q = 1001;n_iter_P =Int(round(4.5/5*n_iter_Q));
    N = n_iter_P * n_iter_Q;
    ArrQ = range(-width, stop = width, length = n_iter_Q)
    ArrP = range(-height, stop = height, length = n_iter_P)
    mesh = [(Q, P) for P in ArrP, Q in ArrQ]
    mesh_list = reshape(mesh, 1, :)
    t_end = t_end * ones(n_iter_Q,n_iter_P);
 
    H = range(.15, stop = 0.41, length = 31)
    H=.2
    # H=.25
    # location = "/mnt/bdd38f66-9ece-451a-b915-952523c139d2/Escape/"
    # location_fig="~/Desktop/PSS_MOVIE/"

    
    count = 1
end

while count <= length(H)
    println(count)
    println(H[count])
    @everywhere Energy = H[count] * ones(N)
    @everywhere h = replace(@sprintf("%.13f",H[count]), "." => "_")
    @everywhere nQ = @sprintf("_%d",n_iter_Q)
    @everywhere nP = @sprintf("_%d",n_iter_P)
    @everywhere h_title = @sprintf("h= %.3f",H[count])

    @everywhere h_BSON="BACK_FORWARD_Hit_for_escape_data"*h*nQ*nP*".bson"
    @everywhere h_BSON_small="BACK_FORWARD_small_Hit_for_escape_data"*h*nQ*nP*".bson"


    @everywhere file_name = "HIT_Escape_" * h * ".fig"
    println(file_name)
    num_until_exit = @showprogress pmap(escape_dimer.escape_exit_num_dimer, mesh_list, t_end, Energy)
    SAVE_DATA=Dict("escape_times"=>num_until_exit, "mesh"=>mesh, "ArrP"=>ArrP, "ArrQ"=>ArrQ, "Energy"=>H)
    @save h_BSON SAVE_DATA

    @everywhere begin
        Q_0 = zeros(0);P_0 = zeros(0)
        Q_1 = zeros(0);P_1 = zeros(0)
        Q_2 = zeros(0);P_2 = zeros(0)
        Q_3 = zeros(0);P_3 = zeros(0)
        Q_4 = zeros(0);P_4 = zeros(0)
        Q_5 = zeros(0);P_5 = zeros(0)
        Q_6 = zeros(0);P_6 = zeros(0)
        Q_7 = zeros(0);P_7 = zeros(0)
        Q_8 = zeros(0);P_8 = zeros(0)
        Q_9 = zeros(0);P_9 = zeros(0)

        Q_b = zeros(0);P_b = zeros(0)
        Q_n = zeros(0);P_n = zeros(0)

    end

   

    for i = 1:N
        Q, P = mesh_list[i]
        if num_until_exit[i] == -2
            push!(Q_n, Q)
            push!(P_n, P)
        elseif num_until_exit[i] == -1
            push!(Q_b, Q)
            push!(P_b, P)
        elseif num_until_exit[i] == 0
            push!(Q_0, Q)
            push!(P_0, P)
        elseif num_until_exit[i] == 1
            push!(Q_1, Q)
            push!(P_1, P)
        elseif num_until_exit[i] == 2
            push!(Q_2, Q)
            push!(P_2, P)
        elseif num_until_exit[i] == 3
            push!(Q_3, Q)
            push!(P_3, P)
        elseif num_until_exit[i] == 4
            push!(Q_4, Q)
            push!(P_4, P)
        elseif num_until_exit[i] == 5
            push!(Q_5, Q)
            push!(P_5, P)
        elseif num_until_exit[i] == 6
            push!(Q_6, Q)
            push!(P_6, P)
        elseif num_until_exit[i] == 7
            push!(Q_7, Q)
            push!(P_7, P)
        elseif num_until_exit[i] == 8
            push!(Q_8, Q)
            push!(P_8, P)
        elseif num_until_exit[i] == 9
            push!(Q_9, Q)
            push!(P_9, P)
        end
    end
    SAVE_DATA_small=Dict("never"=>(Q_n,P_n), "bound"=> (Q_b,P_b), "zero"=>(Q_0,P_0), "one"=>(Q_1,P_1), "two"=>(Q_2,P_2), "three"=>(Q_3,P_3), "four"=>(Q_4,P_4), "five"=>(Q_5,P_5),"six"=>(Q_6,P_6),"seven"=>(Q_7,P_7),"eight"=>(Q_8,P_8),"nine"=>(Q_9,P_9))
    @save h_BSON_small SAVE_DATA_small

    Qn=vcat(Q_n,-Q_n)
    Pn=vcat(P_n,-P_n)
    
    Qb=vcat(Q_b,-Q_b)
    Pb=vcat(P_b,-P_b)
    
    Q0=vcat(Q_0,-Q_0)
    P0=vcat(P_0,-P_0)
    
    Q1=vcat(Q_1,-Q_1)
    P1=vcat(P_1,-P_1)
    
    Q2=vcat(Q_2,-Q_2)
    P2=vcat(P_2,-P_2)
    
    Q3=vcat(Q_3,-Q_3)
    P3=vcat(P_3,-P_3)
    
    Q4=vcat(Q_4,-Q_4)
    P4=vcat(P_4,-P_4)
    
    Q5=vcat(Q_5,-Q_5)
    P5=vcat(P_5,-P_5)

    Q6=vcat(Q_6,-Q_6)
    P6=vcat(P_6,-P_6)

    Q7=vcat(Q_7,-Q_7)
    P7=vcat(P_7,-P_7)

    Q8=vcat(Q_8,-Q_8)
    P8=vcat(P_8,-P_8)

    Q9=vcat(Q_9,-Q_9)
    P9=vcat(P_9,-P_9)

    # Qn=vcat(Q_n,Q_n,-Q_n,-Q_n)
    # Pn=vcat(P_n,-P_n,P_n,-P_n)
    
    # Qb=vcat(Q_b,Q_b,-Q_b,-Q_b)
    # Pb=vcat(P_b,-P_b,P_b,-P_b)
    
    # Q0=vcat(Q_0,Q_0,-Q_0,-Q_0)
    # P0=vcat(P_0,-P_0,P_0,-P_0)
    
    # Q1=vcat(Q_1,Q_1,-Q_1,-Q_1)
    # P1=vcat(P_1,-P_1,P_1,-P_1)
    
    # Q2=vcat(Q_2,Q_2,-Q_2,-Q_2)
    # P2=vcat(P_2,-P_2,P_2,-P_2)
    
    # Q3=vcat(Q_3,Q_3,-Q_3,-Q_3)
    # P3=vcat(P_3,-P_3,P_3,-P_3)
    
    # Q4=vcat(Q_4,Q_4,-Q_4,-Q_4)
    # P4=vcat(P_4,-P_4,P_4,-P_4)
    
    # Q5=vcat(Q_5,Q_5,-Q_5,-Q_5)
    # P5=vcat(P_5,-P_5,P_5,-P_5)

    # Q6=vcat(Q_6,Q_6,-Q_6,-Q_6)
    # P6=vcat(P_6,-P_6,P_6,-P_6)

    # Q7=vcat(Q_7,Q_7,-Q_7,-Q_7)
    # P7=vcat(P_7,-P_7,P_7,-P_7)

    # Q8=vcat(Q_8,Q_8,-Q_8,-Q_8)
    # P8=vcat(P_8,-P_8,P_8,-P_8)

    # Q9=vcat(Q_9,Q_9,-Q_9,-Q_9)
    # P9=vcat(P_9,-P_9,P_9,-P_9)

    SAVE_DATA=Dict("never"=>(Qn,Pn), "bound"=> (Qb,Pb), "zero"=>(Q0,P0), "one"=>(Q1,P1), "two"=>(Q2,P2), "three"=>(Q3,P3), "four"=>(Q4,P4), "five"=>(Q5,P5),"six"=>(Q_6,P_6),"seven"=>(Q_7,P_7),"eight"=>(Q_8,P_8),"nine"=>(Q_9,P_9))
    @save h_BSON SAVE_DATA

    COLOR=["#393b79" ,"#5254a3","#6b6ecf","#9c9ede" ,"#637939","#8ca252" ,"#b5cf6b" ,"#cedb9c" ,"#8c6d31","#bd9e39" ,"#e7ba52","#e7cb94","#843c39","#ad494a" ,"#d6616b","#e7969c" ,"#7b4173" ,"#a55194","#ce6dbd" ,"#de9ed6"]

    mat"""
    figure('position',[0, 0, 1500, 1500]);
    hold on;
    pbaspect([2 2 1])
    """
    mat"axis([-$width, $width, -$height, $height])"
    # mat"figure();set(gcf, 'Position',  [0, 0, 1500, 1500]); hold on;"
    mat"title($h_title)"

    # mat"plot($Qb,$Pb ,'k.','MarkerSize',8)"

    # mat"plot($Q0,$P0 ,'r.','MarkerSize',6)"

    # mat"plot($Q1,$P1 ,'b.','MarkerSize',6)"

    # mat"plot($Q2,$P2 ,'g.','MarkerSize',6)"

    # mat"plot($Q3,$P3 ,'m.','MarkerSize',6)"

    # mat"plot($Q4,$P4 ,'c.','MarkerSize',6)"

    # mat"plot($Q5,$P5 ,'y.','MarkerSize',6)"
    
    # colour=COLOR[2]
    # mat"plot($Q6,$P6,'.','MarkerSize',6,'color',$colour)"

    # colour=COLOR[4]
    # mat"plot($Q7,$P7,'.','MarkerSize',6,'color',$colour)"

    # colour=COLOR[6]
    # mat"plot($Q8,$P8,'.','MarkerSize',6,'color',$colour)"

    # colour=COLOR[8]
    # mat"plot($Q9,$P9,'.','MarkerSize',6,'color',$colour)"

  
    mat"plot($Q_b,$P_b ,'k.','MarkerSize',8)"

    mat"plot($Q_0,$P_0 ,'r.','MarkerSize',6)"

    mat"plot($Q_1,$P_1 ,'b.','MarkerSize',6)"

    mat"plot($Q_2,$P_2 ,'g.','MarkerSize',6)"

    mat"plot($Q_3,$P_3 ,'m.','MarkerSize',6)"

    mat"plot($Q_4,$P_4 ,'c.','MarkerSize',6)"

    mat"plot($Q_5,$P_5 ,'y.','MarkerSize',6)"



    mat"savefig($file_name)"
    # mat"close"
    @everywhere global count += 1
end
