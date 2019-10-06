using Roots


# @everywhere function Yfind(Q,P,H)
# @everywhere Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+( 1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
#
#     @everywhere Y_find1(y)=Hamil(BigInt(0),y,Q,P)-H
#     try
#         Y=find_zero(Y_find1,BigFloat(".1"))
#     catch
#         Y=zeros(BigInt,0)
#     end
# end
function YfindBIG(Q,P,H)
    Y_find1(y)=Hamil(BigFloat("0"),y,Q,P)-H
    try
        Y=find_zero(Y_find1,BigFloat(".1"), maxeval=100,maxfnevals=300,tol=1e-15)
     catch
        Y=zeros(0)
    end
end
