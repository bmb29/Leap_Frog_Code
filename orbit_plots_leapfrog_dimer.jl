using DifferentialEquations
using PyCall
# pygui(:qt)
using PyPlot
# pygui(true)

include("leap_frog_definitions.jl")
include("PSS_Definitions_Dimer_X.jl")

function plot_orbits_leapfrog_aref(H, max_time, q2, p2)
    h = 1 / (2 * H)
    a = (2 + h - 2 * sqrt(h + 1)) / h

    q1 = 0
    p1 = P1_find_dimer(q2, p2, H)
    if isempty(p1)
        p1 = P1_find_dimer_second(q2, p2, H)
    end

    Q1 = (q1 + q2) / sqrt(2)
    Q2 = (q1 - q2) / sqrt(2)
    P1 = (p1[1] + p2) / sqrt(2)
    P2 = (p1[1] - p2) / sqrt(2);

    u0 = [Q1;P2;Q2;P1;0;0] # X P Q Y
    Q0 = [Q1,Q2]; P0 = [P1,P2];
    tspan = (0.0, max_time);

    prob = ODEProblem(Eq_of_M_LAB_FRAME, u0, tspan)
    t, A = solve(prob, Vern9(), reltol = 1e-13, abstol = 1e-15, maxiters = 1e15);
    N = length(t.t)
    Q1 = A[1,:]
    Q2 = A[3,:]
    P1 = A[4,:]
    P2 = A[2,:]

    q1 = (Q1 + Q2) / sqrt(2)
    q2 = (Q1 - Q2) / sqrt(2)
    p1 = (P1 + P2) / sqrt(2)
    p2 = (P1 - P2) / sqrt(2);

    figure(dpi = 250)
    plot(Q1, P1, "b", linewidth = .1)
    plot(Q2, P2, "r", linewidth = .1)
    title("Aref coordinates");
end


function plot_orbits_leapfrog_dimer(H, max_time, q2, p2)
    h = 1 / (2 * H)
    a = (2 + h - 2 * sqrt(h + 1)) / h
    
    q1 = 0
    p1 = P1_find_dimer(q2, p2, H)
    if isempty(p1)
        p1 = P1_find_dimer_second(q2, p2, H)
    end
    
    Q1 = (q1 + q2) / sqrt(2)
    Q2 = (q1 - q2) / sqrt(2)
    P1 = (p1[1] + p2) / sqrt(2)
    P2 = (p1[1] - p2) / sqrt(2);
    
    u0 = [Q1;P2;Q2;P1;0;0] # X P Q Y
    Q0 = [Q1,Q2]; P0 = [P1,P2];
    tspan = (0.0, max_time);
    
    prob = ODEProblem(Eq_of_M_LAB_FRAME, u0, tspan)
    t, A = solve(prob, Vern9(), reltol = 1e-13, abstol = 1e-15, maxiters = 1e15);
    N = length(t.t)
    Q1 = A[1,:]
    Q2 = A[3,:]
    P1 = A[4,:]
    P2 = A[2,:]
    
    q1 = (Q1 + Q2) / sqrt(2)
    q2 = (Q1 - Q2) / sqrt(2)
    p1 = (P1 + P2) / sqrt(2)
    p2 = (P1 - P2) / sqrt(2);
    
    figure(dpi = 250)
    plot(q1, p1, "b", linewidth = .1)
    plot(q2, p2, "r", linewidth = .1)
    title("Dimer coordinates");
end

function plot_orbits_leapfrog_lab(H, max_time, q2, p2)
    h = 1 / (2 * H)
    a = (2 + h - 2 * sqrt(h + 1)) / h
    
    q1 = 0
    p1 = P1_find_dimer(q2, p2, H)
    if isempty(p1)
        p1 = P1_find_dimer_second(q2, p2, H)
    end
    
    Q1 = (q1 + q2) / sqrt(2)
    Q2 = (q1 - q2) / sqrt(2)
    P1 = (p1[1] + p2) / sqrt(2)
    P2 = (p1[1] - p2) / sqrt(2);
    
    u0 = [Q1;P2;Q2;P1;0;0] # X P Q Y
    Q0 = [Q1,Q2]; P0 = [P1,P2];
    tspan = (0.0, max_time);
    
    prob = ODEProblem(Eq_of_M_LAB_FRAME, u0, tspan)
    t, A = solve(prob, Vern9(), reltol = 1e-13, abstol = 1e-15, maxiters = 1e15);
    N = length(t.t)
    Centeroid = A[5,:] + im * A[6,:];
    Z = A[1,:] + im * A[2,:];
    W = A[3,:] + im * A[4,:];
    lin_impulse = im * (1 + a ) * ones(N);
    z_1pos = .5 * (Centeroid + Z + lin_impulse + W);
    z_1neg = .5 * (Centeroid + Z - lin_impulse - W);
    z_2pos = .5 * (Centeroid - Z + lin_impulse - W);
    z_2neg = .5 * (Centeroid - Z - lin_impulse + W);

    z_1posX = real(z_1pos)
    z_1posY = imag(z_1pos)

    z_2posX = real(z_2pos)
    z_2posY = imag(z_2pos)

    z_1negX = real(z_1neg)
    z_1negY = imag(z_1neg)

    z_2negX = real(z_2neg)
    z_2negY = imag(z_2neg)

    figure(dpi = 250)
    plot(z_1posX, z_1posY, c = "b", linewidth = 1)
    plot(z_1negX, z_1negY, c = "r", linewidth = 1)
    plot(z_2posX, z_2posY, c = "b", linewidth = .5, linestyle = "--")
    plot(z_2negX, z_2negY, c = "r", linewidth = .5, linestyle = "--")
    ylim(-3.2, 3.2)
    axis("equal")
    title("Lab coordinates");
end