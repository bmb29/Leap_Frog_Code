# using Distributed
# addprocs(1)
 include("par_module.jl")


 @everywhere begin
    include("par_module.jl")
    push!(LOAD_PATH, "/Users/brandonbehring/Desktop/Leap_Frog_2019")
    using .par_module
    N = 1000
    E = ones(N) * .5
    A = 1:N;
    B = A + ones(N)
end

pmap(par_module.f, A, B, E)
