function [] = MILP()
global T s_0 targets
num_target_tot = 10;
computation_time = zeros(1,num_target_tot);
for num_target = 1:num_target_tot
    mean_time = zeros(1,10);
    for i = 1:10
        tic
        define_parameters(num_target);
        [A_eq, B_eq] = generate_equality_constraints(T,num_target,s_0);
        [A_ineq,B_ineq] = generate_inequality_constraints(T,num_target,targets);
        f = cost_function();
        c_type = [repmat('C',1,6*T),repmat('B',1,num_target*T),'C'];
        [x,fval,exitflag,output] = cplexmilp(f,A_ineq,B_ineq,A_eq,B_eq,[],[],[],[],[],c_type);
        mean_time(i) = toc;
    end
    computation_time(num_target) = mean(mean_time);
end

plot(1:num_target_tot,computation_time);
% if exitflag == 1
%     pos_x = x(1:4:4*ceil((fval+dt)/dt));
%     pos_y = x(2:4:4*ceil((fval+dt)/dt));
%     plot(pos_x,pos_y); hold on
% end
% scatter(targets(1,:),targets(2,:));
% title('Trajectory');
% axis([min(pos_x)-1 max(pos_x)+1 min(pos_y)-1 max(pos_y)+1]);
end

function [] = define_parameters(num_target)
global v_max w_max m f_max N A B dt M T K targets s_0 max_dist
% initial state (x_pos,y_pos,x_vel,y_vel)
s_0 = [0;0;10;0];
% set of targets
max_dist = 50;
targets = -max_dist + 2*max_dist*rand(2,num_target);
K = num_target;
T = 100; % number of time steps
dt = 0.1; % time step
v_max = 15; % m/s
% Minimum turn radius = v^2/(g*sqrt(n^2-1));
R = v_max^2/(9.8*sqrt(4-1));
w_max = 15/R; % rad/sec
m = 1; % kg
f_max = w_max*m*v_max;
N = 4;
M = 2*max_dist + 1;
A = [1 0 dt 0;
    0 1 0 dt;
    0 0 1 0;
    0 0 0 1];
B = [0 0;
    0 0;
    dt 0;
    0 dt];
end

function [A_eq, B_eq] = generate_equality_constraints(T,K,init)
global A B
    A_eq_dyn = zeros(4*(T-1),6*T);
    B_eq_dyn = zeros(4*(T-1),1);
    for i = 0:(T-2)
        A_eq_dyn(4*i+1:4*i+4,4*i+1:4*i+4) = A;
        A_eq_dyn(4*i+1:4*i+4,4*i+5:4*i+8) = -eye(4);
        A_eq_dyn(4*i+1:4*i+4,4*T+2*i+1:4*T+2*i+2) = B;
    end
    A_eq_b = zeros(K, K*T+1);
    B_eq_b = ones(K,1);
    for i = 0:K-1
        A_eq_b(i+1,T*i+1:T*(i+1)) = ones(1,T);
    end
    
    A_eq_init = [eye(4), zeros(4,6*T+K*T+1-4)];
    B_eq_init = init;
    
    A_eq = [blkdiag(A_eq_dyn,A_eq_b);A_eq_init];
    B_eq = [B_eq_dyn; B_eq_b;B_eq_init];
end

function [A_ineq, B_ineq] = generate_inequality_constraints(T,K,targets)
global M N f_max v_max dt
A_ineq_state = zeros(4*T*K,(6+K)*T+1);
B_ineq_state = zeros(4*T*K,1);
for i = 0:T-1
    for k = 0:K-1
        A_ineq_state(4*K*i+4*k+1,i*4+1) = 1;
        A_ineq_state(4*K*i+4*k+1,6*T+k*T+i+1) = M;
        B_ineq_state(4*K*i+4*k+1,1) = M + targets(1,k+1);
        A_ineq_state(4*K*i+4*k+2,i*4+1) = -1;
        A_ineq_state(4*K*i+4*k+2,6*T+k*T+i+1) = M;
        B_ineq_state(4*K*i+4*k+2,1) = M - targets(1,k+1);
        A_ineq_state(4*K*i+4*k+3,i*4+2) = 1;
        A_ineq_state(4*K*i+4*k+3,6*T+k*T+i+1) = M;
        B_ineq_state(4*K*i+4*k+3,1) = M + targets(2,k+1);
        A_ineq_state(4*K*i+4*k+4,i*4+2) = -1;
        A_ineq_state(4*K*i+4*k+4,6*T+k*T+i+1) = M;
        B_ineq_state(4*K*i+4*k+4,1) = M - targets(2,k+1);
    end
end

A_ineq_force = zeros(T*N,(6+K)*T+1);
A_ineq_velocity = zeros(T*N,(6+K)*T+1);
B_ineq_force = f_max*ones(T*N,1);
B_ineq_velocity = v_max*ones(T*N,1);
for i =0:T-1
    for n = 1:N
        A_ineq_force(N*i+n,4*T+2*i+1) = sin(2*pi*n/N);
        A_ineq_force(N*i+n,4*T+2*i+2) = cos(2*pi*n/N);
        A_ineq_velocity(N*i+n,4*i+3) = sin(2*pi*n/N);
        A_ineq_velocity(N*i+n,4*i+4) = cos(2*pi*n/N);
    end
end

A_ineq_T = zeros(K,(6+K)*T+1);
B_ineq_T = zeros(K,1);

for k = 1:K
    A_ineq_T(k,6*T+(k-1)*T+1:6*T+k*T) = (0:T-1)*dt;
    A_ineq_T(k,(6+K)*T+1) = -1; 
end
A_ineq = [A_ineq_state;A_ineq_force;A_ineq_velocity;A_ineq_T];
B_ineq = [B_ineq_state;B_ineq_force;B_ineq_velocity;B_ineq_T];

end

function f = cost_function()
global T K
f = [zeros(1,(6+K)*T),1];
end



