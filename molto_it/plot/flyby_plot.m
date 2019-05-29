%%--------------------------------------------------------------------------
function [DV, aux] = flyby_plot(x, setup)
%%--------------------------------------------------------------------------
%
lc = 149597870.700e03;
mu = 132712440018e09;
tc = sqrt(lc^3/mu);
vc = lc / tc;
%
% Initial spiral position
%
theta0 = setup.theta0;
r0 = setup.r0;
v0 = setup.v0;
psi0 = setup.psi0;
%
% Initial planet state
%
planet1 = setup.planet1;
vfb = setup.vfb;
psifb = setup.psifb;
%
% Perform Flyby
%
[v0, psi0] = flyby(v0, psi0, vfb, psifb, mu, planet1, x(4));
%
% Arrival planet position
%
ToF = setup.ToF + x(2);
[rf, thetaf, vf, psif] = oe2polar(setup.oe, ToF);
%
% Plot Arrival Planet
%
theta_init = linspace(0, 2*pi, 1000);
%
r = setup.oe.SMA * (1 - setup.oe.ECC^2) ./ (1 + setup.oe.ECC * cos(-setup.oe.LNODE-setup.oe.ARGP+theta_init));
%
[xx, yy] = pol2cart(theta_init, r);
plot(xx, yy, 'black -.')
hold on
grid on
axis equal
%

%
% Ensure thetaf > theta0
%
while (thetaf < theta0)
    thetaf = thetaf + 2 * pi;
end
%
thetaf = thetaf + 2 * pi * setup.n;
%
ee1 = x(1);
ee2 = 0;
xA = x(3);
%
thetaA = xA * (thetaf - theta0) + theta0;
thetaB = thetaf;
%
if setup.type == 1 || setup.type == 2
    ee2 = x(5);
    xB = x(6);
    thetaB = xB * (thetaf - thetaA) + thetaA;
end
%
setup.ee1 = ee1;
setup.ee2 = ee2;
setup.theta0 = theta0;
setup.thetaA = thetaA;
setup.thetaB = thetaB;
setup.r0 = r0;
setup.v0 = v0;
setup.psi0 = psi0;
setup.rf = rf;
setup.vf = vf;
setup.thetaf = thetaf;
%
% Obtain DV
%
[DV, vsp, psisp, at, t, data_all] = get_DV_plot(setup);
%
% Relative Velocity
%
v1 = vsp(end) * cos(psisp(end)) - vf * cos(psif);
v2 = vsp(end) * sin(psisp(end)) - vf * sin(psif);
%
aux.thetasp = thetaf;
aux.rsp = rf;
aux.vsp = vsp;
aux.psisp = psisp;
aux.vp = vf;
aux.psip = psif;
aux.ToF = ToF;
aux.ToF_year = ToF * tc / (365 * 24 * 3600);
aux.at = at;
aux.t = t;
aux.vrel = sqrt(v1^2+v2^2) * vc / 1e3;
aux.data_all = data_all;
setup.n
%
