function [E_t] = OE_scl_ptb(E_0, t)
    
    a0 = E_0(1); e0 = E_0(2); i0 = E_0(3); Omega0 = E_0(4); omega0 = E_0(5); M0 = E_0(6);
    J_2 = 1.082627e-3;
    R_E = 6378.137;
    mu = 398600.4418;
    C_J2 = 1.5*J_2*R_E^2*sqrt(mu)*a0^(-3.5);
    
    Omega_t = Omega0 - C_J2*cos(i0)/(1-e0^2)^2*t;
    omega_t = omega0 + C_J2*(2-2.5*sin(i0)^2)/(1-e0^2)^2*t;
    M_t = M0 + sqrt(mu/a0^3)*(1+0.5*C_J2*(1-e0^2)^(-1.5)*(3*cos(i0)^2-1))*t;
    E_t = [a0, e0, i0, Omega_t, omega_t, M_t];
end