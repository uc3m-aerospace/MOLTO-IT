%--------------------------------------------------------------------------
function [DV, aux] = departure_obj(x,setup)
%--------------------------------------------------------------------------
%	MOLTO-IT Software Computation Core
%
%	This program is developed at the Universidad Carlos III de Madrid,
%   as part of a PhD program.
%
%   The software and its components are developed by David Morante González
%
%   The program is released under the MIT License
%
%   Copyright (c) 2019 David Morante González
%
%--------------------------------------------------------------------------
%
%    Function that provides the total DV consumed during a departure leg type
%
%--------------------------------------------------------------------------
%
tc  = setup.tc;
%
%Compute Initial date and Flight Time
%
et0_min    = setup.initial_date;
et0_max    = setup.final_date;
et0_factor = setup.xind + x(4);
et02       = et0_min + (et0_max - et0_min)*et0_factor;
ToF1       = (et02 - setup.et0)/tc;
ToF        = setup.ToF + ToF1+ x(2);
%
% Departure planet position
%
[r0, theta0, v0, psi0] = oe2polar(setup.oe0, ToF1);
%
% Arrival planet position
%
[rf, thetaf, vf, psif] = oe2polar(setup.oe, ToF);
%
% Add Launch velocity
%
%if strcmp(setup.planet1,setup.planet2)
%
v1   = v0*cos(psi0) + setup.vinf0*cos(x(5));
v2   = v0*sin(psi0) + setup.vinf0*sin(x(5));
v0   = norm([v1,v2]);
psi0 = atan2(v2,v1);
%
% Ensure thetaf > theta0
%
while( thetaf < theta0 )
    thetaf = thetaf + 2*pi;
end
%
thetaf = thetaf + 2*pi*setup.n;
%
ee1 = x(1);
ee2 = 0;
xA  = x(3);
%
thetaA  = xA*(thetaf-theta0) + theta0;
thetaB  = thetaf;
%
if setup.type == 1
    ee2     = x(6);
    xB      = x(7);
    thetaB  = xB*(thetaf-thetaA) + thetaA;
end
%
setup.ee1    = ee1;
setup.ee2    = ee2;
setup.theta0 = theta0;
setup.thetaA = thetaA;
setup.thetaB = thetaB;
setup.r0     = r0;
setup.v0     = v0;
setup.psi0   = psi0;
setup.thetaf = thetaf;
%
% Obtain DV
%
[DV, vsp, psisp] = get_DV(setup);
%
aux.thetasp = thetaf;
aux.rsp     = rf;
aux.vsp     = vsp;
aux.psisp   = psisp;
aux.vp      = vf;
aux.psip    = psif;
aux.ToF     = ToF;
aux.et_fact = et0_factor;
aux.et      = et02;
%%--------------------------------------------------------------------------













