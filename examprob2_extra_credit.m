% exam problem 2 extra credit
% all in meters here
clear all
clc

% ICs (velo)
x_dot_0 = 0.05;
y_dot_0 = -0.01;
v0 = [x_dot_0; y_dot_0];

% ICs (position
% at station so
r0 = [0; 0];

% constants
mu = 3.986e+14;
Re = 6378e+3;

% radius of orbit (assume circular)
r = 1.4*Re;

% mean motion
n = sqrt(mu/(r^3));

% period of orbit
period = (2*pi)/n;

%% find relative position and velocity 50 mins after deployment
t = 50*60; % 50 minutes

% Euler-hill matrices
A = [4 - 3*cos(n*t), 0;
    6*(sin(n*t) - n*t), 1];

B = [(1/n)*sin(n*t), (2/n)*(1 - cos(n*t));
    (-2/n)*(1 - cos(n*t)), (4/n)*sin(n*t) - 3*t];

A_prime = [3*n*sin(n*t), 0;
        6*n*cos(n*t) - 6*n, 0];

B_prime = [cos(n*t), 2*sin(n*t);
        -2*sin(n*t), 4*cos(n*t) - 3];

% position
r_50 = A*r0 + B*v0;
% velocity
v_50 = A_prime*r0 + B_prime*v0;

fprintf('position at t = 50 min is x = %.4f m, y = %.4f m, distance %.4f m \n', r_50(1), r_50(2), norm(r_50))
fprintf('velo at t = 50 min is vx = %.4f m/s, vy = %.4f m/s, magnitude %.4f m/s\n', v_50(1), v_50(2), norm(v_50))


%% plot the path relative to Omicron

xplot = [];
yplot = [];

tf = 120000*60; %plot for 150 minutes
tplotlist = 1:1:tf;  % use to find 50 minute mark

idx_of_50min = find(tplotlist == 50*60);

% plot every second along period
for t = 1:120:tf
    % rebuild matrices
    A = [4 - 3*cos(n*t), 0;
        6*(sin(n*t) - n*t), 1];

    B = [(1/n)*sin(n*t), (2/n)*(1 - cos(n*t));
        (-2/n)*(1 - cos(n*t)), (4/n)*sin(n*t) - 3*t];
    
    r = A*r0 + B*v0;
    
    % append
    xplot = [xplot r(1)];
    yplot = [yplot r(2)];
end

endidx = length(xplot)-1;
% switch axes
yplot = -1*yplot;
figure1 = figure;
axes1 = axes('Parent', figure1);
hold(axes1, 'on');
plot(yplot, xplot, '-p', 'MarkerIndices', [idx_of_50min idx_of_50min], 'MarkerFaceColor', 'red', 'MarkerSize', 15);
title('Satellite Deployment Path, relative to Omicron')
xlabel('$\hat{y}$ (m)','Interpreter','latex');
ylabel('$\hat{x}$ (m)','Interpreter','latex');
box(axes1,'on');
startstr = 'Start';
endstr = 'End';
text(yplot(1), xplot(1), startstr)
text(yplot(endidx), xplot(endidx), endstr)
axis equal
grid on

%% use derivatives wrt t to find velo on arrival
%want t arrival so
t = tf;

A_prime = [3*n*sin(n*t), 0;
        6*n*cos(n*t) - 6*n, 0];

B_prime = [cos(n*t), 2*sin(n*t);
        -2*sin(n*t), 4*cos(n*t) - 3];


v_arrival = A_prime*r0 + B_prime*v0;

fprintf('arrival velo is vxf = %.4f m/s, vyf = %.4f m/s, magnitude %.4f m/s\n', v_arrival(1), v_arrival(2), norm(v_arrival))