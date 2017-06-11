function [U, X]=MPC(x_init, mission, T_vec)
path=cat(1, x_init, mission);
T=[0 T_vec]';
X=x_init;
U=[];
[m, n]=size(mission);

for i=1:m
    x_cur=path(i, :);
    x_g=path(i+1, :);
    dt=T(i+1)-T(i);
    
    x_df=diff_flat(x_cur, x_g, dt);
    x=X(end, :);
    t=0;
    
    t1=tic;
    
    while not(norm(x(1:2)-x_g(1:2))<1.)
    [u, x_next, optval] = mpc_step(x,x_df,t);
    U=cat(1, U, u');
    X=cat(1, X, x_next);
    t=t+.1;
    x=x_next;
    end
    t2=toc(t1);
    time(i)=t2-t1;
    
    
end
end

function data=xydxdy(coeff, t)
x=coeff(1:4)'*[ones(size(t)); t; t.^2; t.^3];
y=coeff(5:8)'*[ones(size(t)); t; t.^2; t.^3];

dx=coeff(1:4)'*[zeros(size(t)); ones(size(t)); 2*t; 3*t.^2];
dy=coeff(5:8)'*[zeros(size(t)); ones(size(t)); 2*t; 3*t.^2];

data=[x', y', dx', dy'];
end


function x_df=diff_flat(xi, xf, T)

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

x_df=@(t) xydxdy(coeff, t);
end
