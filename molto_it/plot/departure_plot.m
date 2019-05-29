%%--------------------------------------------------------------------------
function [DV, aux] = departure_plot(x, setup)
%%--------------------------------------------------------------------------
%
lc = 149597870.700e03;
mu = 132712440018e09;
tc = sqrt(lc^3/mu);
vc = lc / tc;
%
%Compute Initial date and Flight Time
%
et0_min = setup.initial_date;
et0_max = setup.final_date;
et0_factor = setup.xind + x(4);
et02 = et0_min + (et0_max - et0_min) * et0_factor;
ToF1 = (et02 - setup.et0) / tc;
ToF = setup.ToF + ToF1 + x(2);
%
% Departure planet position
%
[r0, theta0, v0, psi0] = oe2polar(setup.oe0, ToF1);
%
% Plot Departure Planet
%
theta_init = linspace(0, 2*pi, 1000);
%
r = setup.oe0.SMA * (1 - setup.oe0.ECC^2) ./ (1 + setup.oe0.ECC * cos(-setup.oe0.LNODE-setup.oe0.ARGP+theta_init));
%
[xx, yy] = pol2cart(theta_init, r);
figure
plot(xx, yy, 'black -.')
hold on
grid on
axis equal
%
% Arrival planet position
%
[rf, thetaf, vf, psif] = oe2polar(setup.oe, ToF);

%
% Plot Arrival Planet
%
theta_end = linspace(0, 2*pi, 1000);
%
r = setup.oe.SMA * (1 - setup.oe.ECC^2) ./ (1 + setup.oe.ECC * cos(-setup.oe.LNODE-setup.oe.ARGP+theta_end));
%
[xx, yy] = pol2cart(theta_end, r);
plot(xx, yy, 'black -.')
%
% Add Launch velocity
%
%if strcmp(setup.planet1,setup.planet2)
%     %
v1 = v0 * cos(psi0) + setup.vinf0 * 1e3 / vc * cos(x(5));
v2 = v0 * sin(psi0) + setup.vinf0 * 1e3 / vc * sin(x(5));
v0 = norm([v1, v2]);
psi0 = atan2(v2, v1);
%
%elseif rf > r0
%
%    v0   = v0 + setup.vinf0*1e3/vc;
%
%else
%
%    v0   = v0 - setup.vinf0*1e3/vc;
%
%end
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
if setup.type == 1
    ee2 = x(6);
    xB = x(7);
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
setup.thetaf = thetaf;
%
% Obtain DV
%
[DV, vsp, psisp, at, t, data_all] = get_DV_plot(setup);
%
aux.thetasp = thetaf;
aux.rsp = rf;
aux.vsp = vsp;
aux.psisp = psisp;
aux.vp = vf;
aux.psip = psif;
aux.ToF = ToF;
aux.ToF_year = ToF * tc / (3600 * 24 * 365);
aux.et = et02;
aux.et_fact = et0_factor;
aux.at = at;
aux.t = t;
aux.data_all = data_all;
%%--------------------------------------------------------------------------
