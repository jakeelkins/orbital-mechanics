% Lambert's iterator for elliptic orbits

clc
clear all

mu = 0.399e+6;

c = 44101.032;
s = 47544.803;

TOF = 3.26*30.4167*24*3600;

r1 = 9567;
r2 = 41421.574;

% pick transfer type
type = '2'; % 1b or 2
shape = 'elliptic'; % elliptic or hyperbolic

% pick an initial guess 'a'
a = r1+50;

if strcmp(shape,'elliptic')
    alpha_0 = real(2*asin(sqrt(s/(2*a))));
    beta_0 = real(2*asin(sqrt((s-c)/(2*a))));
end
 
if strcmp(shape,'elliptic') && strcmp(type,'2')
    a2 = (((sqrt(mu)*TOF)/((2*pi) - (alpha_0 - sin(alpha_0)) + (beta_0 - sin(beta_0))))^(2))^(1/3);
end

if strcmp(shape,'elliptic') && strcmp(type,'1b')
    a2 = (((sqrt(mu)*TOF)/((2*pi) - (alpha_0 - sin(alpha_0)) - (beta_0 - sin(beta_0))))^(2))^(1/3);
end

while abs(a-a2) > 1
    a = a + 1;
    if strcmp(shape,'elliptic')
    alpha_0 = real(2*asin(sqrt(s/(2*a))));
    beta_0 = real(2*asin(sqrt((s-c)/(2*a))));
    end
 
    if strcmp(shape,'elliptic') && strcmp(type,'2')
    a2 = (((sqrt(mu)*TOF)/((2*pi) - (alpha_0 - sin(alpha_0)) + (beta_0 - sin(beta_0))))^(2))^(1/3);
    end
    
    if strcmp(shape,'elliptic') && strcmp(type,'1b')
    a2 = (((sqrt(mu)*TOF)/((2*pi) - (alpha_0 - sin(alpha_0)) - (beta_0 - sin(beta_0))))^(2))^(1/3);
    end
    
end

abs(a-a2)


(4*a*(s-r1)*(s-r2))/(c^2)
