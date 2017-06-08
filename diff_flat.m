function MPC(x_init, mission, T_vec)
path=cat(1, x_init, mission);
T=[0 T_vec]';
X=x_init;
U=[];

for i=1:len(mission)
    x_cur=path(i, :);
    x_g=path(i+1, :);
    dt=T(i+1)-T(i);
    
    x_df=diff_flat(x_cur, x_g, dt);
    x=X(end, :)
    while not(x==x_g)
    [u, xnext, optval] = mpc_step(x,x_df,dt);
    U=cat(1, U, u);
    X=cat(1, X, x_next);
   
    x=x_next;
    end
end
end


function [fx,fy]=diff_flat(xi, xf, T)

A=[1,0,0,0, 0,0,0,0;...
    0,0,0,0, 1,0,0,0;...
    1,T,T^2,T^3, 0,0,0,0;...
    0,0,0,0, 1,T,T^2,T^3;...
    0,1,0,0, 0,0,0,0;...
    0,0,0,0, 0,1,0,0;...
    0,1,2*T,3*T^2, 0,0,0,0;...
    0,0,0,0, 0,1,2*T,3*T^2];

b=[xi(1),xi(2),xf(1),xf(2), xi(3),xi(4),xf(3),xf(4)];
coeff=inv(A)*b';

x_df=@(t) coeff'*[ones('like',t); t; t^2; t^3; ones('like',t); t; t^2; t^3];
end

