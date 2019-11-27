using Roots

H_test(u)=( (u[3]-u[1])^2+(u[2]-u[4])^2 )*( (u[3]+u[1])^2+(u[2]+u[4])^2 )/ ((u[2]^4+2*u[2]^2*(u[1]^2-1)+(1+u[1]^2)^2 )*(u[3]^4+2*u[3]^2*(u[4]^2+1)+(u[4]^2-1)^2 ))

function YfindNP(Q,P,H)
Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+( 1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))

        Y_find1(y)=Hamil(0,y,Q,P)-H
    try
        Y=find_zeros(Y_find1,0,10, maxeval=100,maxfnevals=300,tol=1e-15)
     catch
        Y=zeros(0)
    end
end
# function YfindNP(Q,P,H)
# Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+( 1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
#
#         Y_find1(y)=Hamil(0,y,Q,P)-H
#     try
#         Y=find_zero(Y_find1,0.01, maxeval=100,maxfnevals=300,tol=1e-15)
#      catch
#         Y=zeros(0)
#     end
# end
