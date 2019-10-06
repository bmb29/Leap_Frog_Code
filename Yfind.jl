using Roots
function Yfind(Q,P,H)

    Hamil(XX,YY,QQ,PP)=( (QQ-XX)^2+(PP-YY)^2 )*( (QQ+XX)^2+(PP+YY)^2 )/((PP^4+2*PP^2*(XX^2-1)+( 1+XX^2)^2 )*(QQ^4+2*QQ^2*(YY^2+1)+(YY^2-1)^2 ))
    Y_find(y)=Hamil(0,y,Q,P)-H
    try
        Y=find_zero(Y_find,.01,maxeval=100,maxfnevals=300,tol=1e-15)
    catch
        Y=zeros(0)
    end
end
# @everywhere function Yfind(Q,P,H)
#     Y_find1(y)=Hamil(0,y,Q,P)-H
#         try
#             Y=find_zeros(Y_find1,0,6, maxeval=100,maxfnevals=300,tol=1e-12)
#          catch
#             Y=zeros(0)
#         end
#     end
