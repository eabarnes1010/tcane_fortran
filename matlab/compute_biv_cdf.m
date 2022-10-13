u       = 4.5;
v       = 1.5;
mu_u    = 2.5;
mu_v    = 1.0;
sigma_u = 2.0;
sigma_v = 3.0;
rho     = 0.5;

X = [u; v];
Mu = [mu_u; mu_v];
Sigma = [sigma_u^2, rho*sigma_u*sigma_v; rho*sigma_u*sigma_v, sigma_v^2];

r_sqr = (X - Mu)' * inv(Sigma) * (X - Mu);
F = 1 - exp(-r_sqr/2);
disp(F)

u       = [ 8.0, 9.0, 2.0, 9.0 ];
v       = [ 4.0, 9.0, 8.0, 9.0 ];
mu_u    = [ 4.0, 4.0, 7.0, 8.0 ];
mu_v    = [ 2.0, 5.0, 5.0, 6.0 ];
sigma_u = [ 9.0, 4.0, 6.0, 3.0 ];
sigma_v = [ 6.0, 2.0, 2.0, 5.0 ];
rho     = [ 0.2, 0.3, 0.8, 0.3 ];

for i = 1:4
    X = [u(i); v(i)];
    Mu = [mu_u(i); mu_v(i)];
    Sigma = [sigma_u(i)^2, rho(i)*sigma_u(i)*sigma_v(i); rho(i)*sigma_u(i)*sigma_v(i), sigma_v(i)^2];
    
    r_sqr = (X - Mu)' * inv(Sigma) * (X - Mu);
    F(i) = 1 - exp(-r_sqr/2);
end
F'

u       = [ 8.0, 9.0, 2.0, 9.0 ];
v       = [ 4.0, 9.0, 8.0, 9.0 ];
mu_u    = 4.0;
mu_v    = 2.0;
sigma_u = 9.0;
sigma_v = 6.0;
rho     = 0.2;

for i = 1:4
    X = [u(i); v(i)];
    Mu = [mu_u; mu_v];
    Sigma = [sigma_u^2, rho*sigma_u*sigma_v; rho*sigma_u*sigma_v, sigma_v^2];
    
    r_sqr = (X - Mu)' * inv(Sigma) * (X - Mu);
    F(i) = 1 - exp(-r_sqr/2);
end
F'