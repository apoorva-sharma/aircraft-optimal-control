function [u, xnext, optval] = mpc_step(x,x_df, t)
    Q = diag([1,1,0.1,0.1]);
    Qhalf = sqrtm(Q);
    
    v_max=15.5;
    % Turn radius = 20m
    w_max = v_max/20; % rad/sec
    m = 1; % kg
    f_max = w_max*m*v_max;
    
    N = 4;
    
    dt = .1;
    T=5;
    
    A = [1 0 dt 0;
         0 1 0 dt;
         0 0 1 0;
         0 0 0 1];
    B = [0 0;
         0 0;
         dt 0;
         0 dt];
    
    A_force = zeros(N,2);
    for n=1:N
        A_force(n,:) = [sin(2*pi*n/N) cos(2*pi*n/N)];
    end
    force_lim = f_max*ones(N,1);
    
    A_vel = zeros(N,4);
    for n=1:N
        A_vel(n,:) = [0 0 sin(2*pi*n/N) cos(2*pi*n/N)];
    end
    vel_lim = v_max*ones(N,1);

    desired_traj = x_df(t:dt:t+dt*T)'
    
    cvx_precision(max(min(abs(x))/35,1e-6))
    cvx_begin 
        variables X(4,T+1) U(2, T)
        
        max((A_force*U)')' <= force_lim;
        max((A_vel*X)')' <= vel_lim;
        
        X(:,2:T+1) == A*X(:,1:T)+B*U;
        X(:,1) == x';
        minimize ( norm(Qhalf*(X(:,1:T) - desired_traj(:,1:T)),'fro') )
    cvx_end
    
    u = U(:, 1);
    xnext = (A*x'+B*u)'+[normrnd(0, .1), normrnd(0, .1), 0, 0];
    
    optval = cvx_optval;
    if strcmp(cvx_status, 'Infeasible' ) == 1
        optval = Inf;
    end
end