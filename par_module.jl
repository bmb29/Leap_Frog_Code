
module par_module
export f
using Roots
function f(k,b,c)
    b*[k k^2; k^2 k]+c*[1 1; 1 1]
end
end
