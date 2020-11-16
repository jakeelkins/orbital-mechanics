% Elliptic orbital transfer solver
% Jake Elkins, AEM569
clear all
clc
format long

Re = 6378;
mu = 0.399e+6;

% input: known parameters of transfer
%---beginning conditions---
a1 = 1.5*Re;
e1 = 0.0;
i1 = 35*(pi/180);
RAAN1 = 45*(pi/180);
arg_peri1 = 0*(pi/180);
% assume on ECI x-axis
theta1 = 0;
TA1 = 0;  % would need edited if init not a circle
%---end conditions---
a2 = 5.3*Re;
e2 = 0.25;
i2 = 35*(pi/180);
RAAN2 = 45*(pi/180);
arg_peri2 = 60*(pi/180);

theta2 = 260*(pi/180);
TA2 = theta2 - arg_peri2;
% -----

% desired TOF of transfer, for Lamberts algorithm
TOF_des = 15*60*60;

%% transfer type determination

% find transfer angle in transfer plane by theta of end condition
transfer_angle = theta2 - theta1;

%if TA > 180, type II, else, type I
if transfer_angle > pi
    type = 2;
    fprintf('transfer angle is %.2f rad, > pi, so type II transfer\n', transfer_angle)
else
    type = 1;
    fprintf('transfer angle is %.2f rad, < pi, so type I transfer\n', transfer_angle)
end

% find r1, r2:
if e1 == 0.0
    r1 = a1;
else
    P1 = a1*(1 - e1^2);
    r1 = P1/(1 + e1*cos(TA1));
end

if e2 == 0.0
    r2 = a2;
else
    P2 = a2*(1 - e2^2);
    r2 = P2/(1 + e2*cos(TA2));
end

DCM1 = [cos(RAAN1)*cos(theta1) - sin(RAAN1)*cos(i1)*sin(theta1), -cos(RAAN1)*sin(theta1) - sin(RAAN1)*cos(i1)*cos(theta1), sin(RAAN1)*sin(i1);
    sin(RAAN1)*cos(theta1) + cos(RAAN1)*cos(i1)*sin(theta1), -sin(RAAN1)*sin(theta1) + cos(RAAN1)*cos(i1)*cos(theta1), -cos(RAAN1)*sin(i1);
    sin(i1)*sin(theta1), sin(i1)*cos(theta1), cos(i1)];

% form DCM for when we need it and for transfer angle calcs
DCM2 = [cos(RAAN2)*cos(theta2) - sin(RAAN2)*cos(i2)*sin(theta2), -cos(RAAN2)*sin(theta2) - sin(RAAN2)*cos(i2)*cos(theta2), sin(RAAN2)*sin(i2);
    sin(RAAN2)*cos(theta2) + cos(RAAN2)*cos(i2)*sin(theta2), -sin(RAAN2)*sin(theta2) + cos(RAAN2)*cos(i2)*cos(theta2), -cos(RAAN2)*sin(i2);
    sin(i2)*sin(theta2), sin(i2)*cos(theta2), cos(i2)];

pos_vec1 = DCM1*[r1;0;0];
pos_vec2 = DCM2*[r2;0;0];
transfer_angle = acos(dot(pos_vec1, pos_vec2)/(norm(pos_vec1)*norm(pos_vec2)));

if type==2
    transfer_angle = 2*pi - transfer_angle;
end

%% solve space triangle

% law of cosines for chord
if transfer_angle > pi
    loc_arg = 2*pi - transfer_angle;
else
    loc_arg = transfer_angle;
end

c = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(loc_arg));

% semiparameter of space triangle
s = 0.5*(r1 + r2 + c);

% find TOF of parabolic transfer to compare:
if type==1
    TOF_para = (1/3)*(sqrt(2/mu))*(s^(3/2) - (s - c)^(3/2));
else
    TOF_para = (1/3)*(sqrt(2/mu))*(s^(3/2) + (s - c)^(3/2));
end

% now compared desired TOF to TOF parabolic. if greater, need an ellipse
% (wanna go slower). Else, need hyperbola (go faster).
if TOF_des > TOF_para
    shape = 'elliptic';
    fprintf('transfer is elliptic based on TOF para\n')
else
    shape = 'hyperbolic';
    fprintf('transfer is hyperbolic based on TOF para\n')
end

% -now compute minimum energy transfer TOF to find secondary type-
a_min = s/2;

% prinicpal vals for TOF calc
alpha_0 = 2*asin(sqrt(s/(2*a_min)));
beta_0 = 2*asin(sqrt((s-c)/(2*a_min)));

if type==1
    alpha = alpha_0;
    beta = beta_0;
else
    alpha = alpha_0;
    beta = -beta_0;
end

TOF_min = sqrt((a_min^3)/mu)*(alpha - sin(alpha) - (beta - sin(beta)));

% compare. if TOF_des > TOF_min, type B. else, type A.
if TOF_des > TOF_min
    sec_type = 'B';
    fprintf('transfer is sec type B based on TOF min\n')
else
    sec_type = 'A';
    fprintf('transfer is sec type A based on TOF min\n')
end

fprintf('final transfer type calculated: %d, %s, %s \n', type, shape, sec_type)

%% finding orbital parameters of the transfer
% first: iterate on 'a' to find SMA of transfer orbit

% pick an initial guess 'a'. try a_min
a = a_min;

% prinicpal vals for TOF calc
alpha_0 = 2*asin(sqrt(s/(2*a)));
beta_0 = 2*asin(sqrt((s-c)/(2*a)));

% find alphas for the transfer types
if type==1 && strcmp(sec_type,'A')
    alpha = alpha_0;
    beta = beta_0;
end

if type==1 && strcmp(sec_type,'B')
    alpha = 2*pi - alpha_0;
    beta = beta_0;
end

if type==2 && strcmp(sec_type,'A')
    alpha = alpha_0;
    beta = -beta_0;
end

if type==2 && strcmp(sec_type,'B')
    alpha = 2*pi - alpha_0;
    beta = -beta_0;
end

% Lambert's solved for a. we'll iterate with this
a2 = ((sqrt(mu)*TOF_des)/(alpha - beta - (sin(alpha) - sin(beta))))^(2/3);

while abs(a-a2) > 1
    a = a + 0.1;  % adjust these for fidelity to converge
    
    % prinicpal vals for TOF calc
    alpha_0 = 2*asin(sqrt(s/(2*a)));
    beta_0 = 2*asin(sqrt((s-c)/(2*a)));
    
    if type==1 && strcmp(sec_type,'A')
        alpha = alpha_0;
        beta = beta_0;
    end

    if type==1 && strcmp(sec_type,'B')
        alpha = 2*pi - alpha_0;
        beta = beta_0;
    end

    if type==2 && strcmp(sec_type,'A')
        alpha = alpha_0;
        beta = -beta_0;
    end

    if type==2 && strcmp(sec_type,'B')
        alpha = 2*pi - alpha_0;
        beta = -beta_0;
    end

    % Lambert's solved for a. we'll iterate with this
    a2 = ((sqrt(mu)*TOF_des)/(alpha - beta - (sin(alpha) - sin(beta))))^(2/3);
    
end

fprintf('converged on SMA val of %.2f km \n', a)

%% now run thru the other orbital parameter calculations for departure

%need type of orbit to get which P (semilatus rectum) to grab
P_guess1 = ((4*a*(s - r1)*(s - r2))/(c^2))*(sin((alpha+beta)/2))^2;
P_guess2 = ((4*a*(s - r1)*(s - r2))/(c^2))*(sin((alpha-beta)/2))^2;

if (type==1 && strcmp(sec_type,'A')) || (type==2 && strcmp(sec_type,'B'))
    P = max(P_guess1, P_guess2);
end

if (type==1 && strcmp(sec_type,'B')) || (type==2 && strcmp(sec_type,'A'))
    P = min(P_guess1, P_guess2);
end

% NOTE: all this is for elliptic. would need to rewrite for hyperbolic

% eccentricity
e = sqrt(1 - (P/a));

% radius at departure
r_D = r1;

% velocity before first maneuver (v1) using orbital energy
v1 = sqrt(2*mu*((1/r1) - (1/(2*a1))));

% velocity at departure (after maneuver) on transfer ellipse now
v_D = sqrt(2*mu*((1/r_D) - (1/(2*a))));

% [!] How do we resolve quadrant ambiguities? [!]
% assume that the first burn is 'up' always if raising orbit
% flight path angle at departure using specific ang mom
FPA_D = acos((sqrt(mu*P))/(r_D*v_D));

% [!] How do we resolve quadrant ambiguities? [!]
% true anomaly at departure on transfer ellipse
rvarg = ((r_D*v_D^2)/mu);
TA_D = atan((rvarg*cos(FPA_D)*sin(FPA_D))/(rvarg*(cos(FPA_D))^2 - 1));

if TA_D < 0
    TA_D = TA_D + pi;
end

TA_D_guess1 = acos((1/e)*((P/r_D) - 1));
if TA_D_guess1 < 0
    TA_D_guess1 = 2*pi + TA_D_guess1;
end

if abs(TA_D - TA_D_guess1) < 0.1
    fprintf('true anomaly difference at dep. is good, delta: %.2f rad \n', abs(TA_D - TA_D_guess1))
else
    fprintf('[!] CHECK: true anomaly difference at dep. [!]\n')
end

%% arrival

% we know the transfer angle and true anom at departure, true anom at
% arrival should be:
TA_A_guess1 = TA_D + transfer_angle;

% radius at arrival, already did that way up top, should be same
r_A = r2;

% velocity at arrival (before maneuver) on transfer ellipse
v_A = sqrt(2*mu*((1/r_A) - (1/(2*a))));

% flight path angle at arrival using specific ang mom
FPA_A = acos((sqrt(mu*P))/(r_A*v_A));

% check if we passed apoapsis. if so, we need the negative FPA, bc
% descending
if TA_A_guess1 > pi && TA_D < pi
    % we passed it. take negative
    FPA_A = -FPA_A;
end

% now calc what true anom at arrival would be with equation and compare
rvarg = ((r_A*v_A^2)/mu);
TA_A = atan((rvarg*cos(FPA_A)*sin(FPA_A))/(rvarg*(cos(FPA_A))^2 - 1));

% resolve the ambiguity
if TA_A_guess1 > pi
    TA_A = TA_A + pi;
end

% correct TA_A_guess1
if TA_A_guess1 > 2*pi
    TA_A_guess1 = TA_A_guess1 - 2*pi;
end

TA_A_guess2 = acos((1/e)*((P/r_A) - 1));
TA_A_guess3 = 2*pi - TA_A_guess2;

%r_A - P/(1 + e*cos(TA_A_guess1))
% now we pick the true anomaly that is closest
% this part loops thru TA_A vals, finds abs of difference, saves the two
% that give that, and avergaes them for the TA_A val we use.
minTAidx1 = 0;
minTAidx2 = 0;
minabsval = 2*pi;
TA_A_array = [TA_A, TA_A_guess1, TA_A_guess2, TA_A_guess3];

for i=1:length(TA_A_array)
    for j=1:length(TA_A_array)
        if i ~= j
            currabsval = abs(TA_A_array(i) - TA_A_array(j));
            if currabsval < minabsval
                minTAidx1 = i;
                minTAidx2 = j;
                minabsval = currabsval;
            end
        end

    end
end

% calc mean of the two most-similar guesses to use of true anom arrival
TA_A = mean([TA_A_array(minTAidx1), TA_A_array(minTAidx2)]);
fprintf('true anom at arrival found, discrep. = %.2f \n', abs(TA_A_array(minTAidx1) - TA_A_array(minTAidx2)))


%% find inclination, RAAN, arg of periapsis of the transfer orbit

% get ang momentum unit vector for inclination
h_hat = cross(pos_vec2, pos_vec1)/norm(cross(pos_vec2, pos_vec1));

% z coord of h hat is cos of inclin: (take positive)
i_transfer = acos(h_hat(3));

% now we can get RAAN
RAAN_transfer1 = real(acos(-h_hat(2)/sin(i_transfer)));
RAAN_transfer2 = asin(h_hat(1)/sin(i_transfer));

% resolve ambiguities
minRAANidx1 = 0;
minRAANidx2 = 0;
minabsval = 2*pi;
RAAN_array = [RAAN_transfer1, pi-RAAN_transfer1, RAAN_transfer2, 2*pi - RAAN_transfer2];

for i=1:length(RAAN_array)
    for j=1:length(RAAN_array)
        if i ~= j
            currabsval = abs(RAAN_array(i) - RAAN_array(j));
            if currabsval < minabsval
                minRAANidx1 = i;
                minRAANidx2 = j;
                minabsval = currabsval;
            end
        end

    end
end

% calc mean of the two most-similar guesses to use of RAAN
RAAN_transfer = mean([RAAN_array(minRAANidx1), RAAN_array(minRAANidx2)]);
fprintf('RAAN found, discrep. = %.2f \n', abs(RAAN_array(minRAANidx1) - RAAN_array(minRAANidx2)))

% now find arg of peri from position at arrival
r_hat2 = pos_vec2/norm(pos_vec2);

% get arg of TA at arrival
theta2_transfer = asin(r_hat2(3)/sin(i_transfer));

if theta2_transfer < 0
    theta2_transfer = 2*pi + theta2_transfer;
end

% and subtract to get arg of peri using arrival (todo: departure too)
arg_peri_transfer = theta2_transfer - TA_A;

%% delta-v calcs at departure

% assume circ orbit at first, so no delta FPA (FPA1 = 0)
% get departure velo in VNC frame
v_D_vnc = v_D.*[cos(i_transfer-i1)*cos(FPA_D), sin(i_transfer-i1), cos(i_transfer-i1)*sin(FPA_D)];

% now just subtract the vectors
delta_v_D = v_D_vnc - ([v1, 0, 0]);

% magnitude
dv_D_magn = norm(delta_v_D);

fprintf('delta_V cost of departure: %.4f km/s \n', dv_D_magn)

%% delta-V calcs at arrival

% need the FPA2 (arrival, of s/c)
% flight path angle at arrival using specific ang mom
P2 = a2*(1 - e2^2);
r2 = P2/(1 + e2*cos(TA2));
% velocity before first maneuver (v1) using orbital energy
v2 = sqrt(2*mu*((1/r2) - (1/(2*a2))));

FPA_2 = acos((sqrt(mu*P2))/(r2*v2));

% resolve the correct FPA since we know TA of the orbit
if (TA2 > pi) && (FPA_2 > pi/2)
   FPA_2 = FPA_2 - pi;
end
if (TA2 > pi) && (FPA_2 < pi/2)
   FPA_2 = -FPA_2;
end
if (TA2 < pi) && (FPA_2 > pi/2)
   FPA_2 = pi - FPA_2 ;
end

% get arrival velo in VNC frame (use deltas of FPA, i)
v_A_vnc = v_A.*[cos(i2-i_transfer)*cos(FPA_2-FPA_A), sin(i2-i_transfer), cos(i2-i_transfer)*sin(FPA_2-FPA_A)];

% now just subtract the vectors
delta_v_A = v_A_vnc - ([v2, 0, 0]);

% magnitude
dv_A_magn = norm(delta_v_A);

fprintf('delta_V cost of arrival: %.4f km/s \n', dv_A_magn)

%% print out total cost, relevant quantities

fprintf('total delta_V cost of maneuver: %.4f km/s \n', dv_A_magn+dv_D_magn)



