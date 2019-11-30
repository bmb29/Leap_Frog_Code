include("LD.jl")
using MATLAB
using ForwardDiff
using ProgressMeter
x=range(.25,stop=1.2,length=Int(1e4))
y1=@showprogress map(Dimer_Lagrangian_Descriptor_Function,x);
logz1=log1p.(abs.(y1))

mat"figure()"
mat"plot($x,$y1)"
# mat"plot($x,$logz1);hold on;"
delta=1e-12
g(x)=( Dimer_Lagrangian_Descriptor_Function(x+delta)-Dimer_Lagrangian_Descriptor_Function(x-delta))/(2*delta) 
f(x)=abs.(( Dimer_Lagrangian_Descriptor_Function(x-delta)-Dimer_Lagrangian_Descriptor_Function(x)+Dimer_Lagrangian_Descriptor_Function(x+delta))/(delta^2))
# g(x)=( Dimer_Lagrangian_Descriptor_Function(x-.001)-Dimer_Lagrangian_Descriptor_Function(x+.001))/(.002)
y2=@showprogress map(g,x);
logz2=log10.(abs.(y2))
mat"figure()"
mat"plot($x,$y2,'k')"

mat"figure()"
mat"plot($x,$logz2,'k.')"




# y3=@showprogress map(f,x)
# # logz=log1p.(abs.(y2))
# # logz3=log1p.(abs.(y3))

# mat"figure()"
# mat"plot($x,$y3)"
# # mat"plot($x,$logz3)"
