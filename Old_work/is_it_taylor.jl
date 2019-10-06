function is_it_taylor(A,H,tol)
    Z=A[:,1]+im *A[:,2]
    W=A[:,3]+im *A[:,4]
    h=1/(2*H)
    a=(2+h-2*sqrt(h+1))/h;
    N=length(A[:,1])
    A=zeros(0)
    lin_impulse=im*(1+a )*ones(N)

    z_1pos=.5*(Z+lin_impulse+W)
    z_1neg=.5*(Z-lin_impulse-W)
    z_2pos=.5*(-Z+lin_impulse-W)
    z_2neg=.5*(-Z-lin_impulse+W)

    #1pos and 1neg bound
    d_11=abs.(z_1pos-z_1neg)
    #1pos and 2neg bound
    d_12=abs.(z_1pos-z_2neg)

    delta_11=abs(d_11[end-10]-d_11[end])
    delta_12=abs(d_12[end-10]-d_12[end])

    Z=zeros(0)
    W=zeros(0)
    lin_impulse=zeros(0)
    z_1pos=zeros(0)
    z_1neg=zeros(0)
    z_2pos=zeros(0)
    z_2neg=zeros(0)
    d_11=zeros(0)
    d_12=zeros(0)


    return delta_11<tol || delta_12<tol
    # return delta_11, delta_12

end
