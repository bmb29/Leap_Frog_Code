@everywhere function is_it_exit2(A,H,tol)
    Z=A[1]+im *A[2]
    W=A[3]+im *A[4]
    h=1/(2*H)
    a=(2+h-2*sqrt(h+1))/h;
    lin_impulse=im*(1+a )

    z_1pos=.5*(Z+lin_impulse+W)
    z_1neg=.5*(Z-lin_impulse-W)
    z_2pos=.5*(-Z+lin_impulse-W)
    z_2neg=.5*(-Z-lin_impulse+W)

    #1pos and 1neg bound
    d_11=abs.(z_1pos-z_1neg)
    #1pos and 2neg bound
    d_12=abs.(z_1pos-z_2neg)

    delta_11=abs(d_11[1]-d_11[2])
    delta_12=abs(d_12[1]-d_12[2])

    return delta_11<tol || delta_12<tol
    # return delta_11, delta_12

end
