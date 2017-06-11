clear all
x_i=[0,0,10,0];
mission=[40.49, 47.97, -6.52, 10.61;...
    -18.89,42.34, -5.88, -10.61 ;...
    -6.98, -31.52, -4.63, -7.55];
tf=[6.5, 12.5, 19.5]; %s
[U, X]=MPC(x_i, mission, tf);

figure(1)
plot(X(:,1), X(:, 2), 'k-', mission(:, 1), mission(:, 2), 'bo', x_i(1), x_i(2), 'r*')