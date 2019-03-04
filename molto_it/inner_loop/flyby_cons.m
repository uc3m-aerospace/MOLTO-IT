%--------------------------------------------------------------------------
function [c,ceq] = flyby_cons(x,setup)
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
%    Function that provides the infeasability of the constraints (equality
%    and inequality) for  a flyby type leg
%
%--------------------------------------------------------------------------
mu = setup.mu;
%
% Initial spiral position
%
theta0 = setup.theta0;
r0     = setup.r0;
v0     = setup.v0;
psi0   = setup.psi0;
%
% Initial planet state
%
planet1 = setup.planet1;
vfb     = setup.vfb;
psifb   = setup.psifb;
%
% Perform Flyby
%
[v0,psi0] = flyby(v0,psi0,vfb,psifb,mu,planet1,x(4));
%
% Arrival planet position
%
ToF                   = setup.ToF + x(2);
[rf, thetaf,vf, psif] = oe2polar(setup.oe, ToF);
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
if numel(x)>5
    ee2     = x(5);
    xB      = x(6);
    thetaB  = xB*(thetaf-thetaA) + thetaA;
end
%
%------------------------------------------------------------------
% FIRST SPIRAL ARC
%------------------------------------------------------------------
%
[t1, v, r, theta, psi ] = propagate_spirals(v0,r0,theta0,psi0,thetaA,ee1) ;
%
coast_arc()
%
if numel(x)>5
    %
    [t2, v, r, ~, psi ] = propagate_spirals( v,r, thetaB, psi, thetaf, ee2);
    tk = tk + t2;
    %
end
%
%------------------------------------------------------------------
% COMPUTE CONSTRAINTS
%------------------------------------------------------------------
%
ToF_spiral = t1 + tk;
%
c(1) = -ToF;
%
if setup.type == 2
    v1   = v(end)*cos(psi(end)) - vf*cos(psif);
    v2   = v(end)*sin(psi(end)) - vf*sin(psif);
    c(2) =  -setup.vinff_max^2 +v1^2+v2^2;
end
%
ceq(1) = (r(end)     - rf)/1;
ceq(2) = (ToF_spiral - ToF)/1;
%
if setup.type == 1
    %
    ceq(3) = (v(end) - vf);
    ceq(4) = (psi(end) - psif)/1;
    %
end
%
if ~isreal(ToF_spiral) || isnan(ToF_spiral)
    %
    ceq = NaN *ones(size(ceq));
    %
end













