% exam problem 2
% all in meters here
clear all
clc

% altitude of orbit
alt = 200e+3;

% ICs (position)
x0 = -100;
y0 = -100;

% constants
mu = 3.986e+14;
Re = 6378e+3;

% radius of orbit (assume circular)
r = Re + alt;

% mean motion
n = sqrt(mu/(r^3));

% period of orbit
period = (2*pi)/n;

% end conditions
tf = 0.2*period;
xf = 0;
yf = 0;

%% euler-hill equations at tf

t = tf;
A = [4 - 3*cos(n*t), 0;
    6*(sin(n*t) - n*t), 1];

B = [(1/n)*sin(n*t), (2/n)*(1 - cos(n*t));
    (-2/n)*(1 - cos(n*t)), (4/n)*sin(n*t) - 3*t];

r0 = [x0; y0];

v0 = linsolve(B, -(A*r0));

fprintf('initial velo needed is vx0 = %.4f m/s, vy0 = %.4f m/s, magnitude %.4f m/s\n', v0(1), v0(2), norm(v0))

%% plot the path relative to the shuttle

xplot = [];
yplot = [];

% plot every second along period
for t = 1:1:tf
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
plot(yplot, xplot)%, '-p', 'MarkerIndices', [1, endidx], 'MarkerFaceColor', 'red', 'MarkerSize', 15);
title('Approach Path, relative to Shuttle')
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

