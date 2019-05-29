function[DV, vf, psif] = get_DV(setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TARGET RADIUS AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ee1 = setup.ee1;
ee2 = setup.ee2;
theta0 = setup.theta0;
thetaA = setup.thetaA;
thetaB = setup.thetaB;
r0 = setup.r0;
v0 = setup.v0;
psi0 = setup.psi0;
thetaf = setup.thetaf;
t0 = 0;
%
%------------------------------------------------------------------
% FIRST SPIRAL ARC
%------------------------------------------------------------------
%
eei = ee1;
K1i = v0^2 - 2 / r0 * (1 - eei);
K2i = r0 * v0^2 * sin(psi0);
%
regime0 = sign(cos(psi0));
angle   = linspace(theta0,thetaA,100);
%
[r_all, v_all, ~ a, e, w, ~ psi_all, DV1, ~] = propagate_spirals_DV(K1i, K2i, eei, regime0, r0, theta0, t0, angle);
%
%------------------------------------------------------------------
% COAST ARC / KEPLERIAN ORBIT
%------------------------------------------------------------------
%
vk = thetaB - w(end);
%
r = a(end) * (1 - e(end)^2) ./ (1 + e(end) * cos(vk));
%
psi  = atan2 (1 + e(end) * cos( vk ), e(end) * sin( vk )) ;
%
%------------------------------------------------------------------
% SECOND SPIRAL ARC
%------------------------------------------------------------------
%
eei = ee2;
%
K1i = (2 * a(end) * eei - r) / (r * a(end));
K2i = sqrt(1+2*e(end)*cos(vk)+e(end)^2);
%
regime0 = sign(cos(psi));
angle   = linspace(thetaB,thetaf,100);
%
[~ v, ~ ~ ~ ~ ~ psi, DV2, ~] = propagate_spirals_DV(K1i, K2i, eei, regime0, r(end), thetaB, t0, angle);
%
%------------------------------------------------------------------
% COMPUTE TOTAL DV
%------------------------------------------------------------------
%
DV = abs(DV1) + abs(DV2);
%
psif = psi(end);
%
vf = v(end);
%
