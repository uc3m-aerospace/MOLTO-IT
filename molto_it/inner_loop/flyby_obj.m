%--------------------------------------------------------------------------
function [DV,aux] = flyby_obj(x,setup)
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
%    Function that provides the total DV consumed during a flyby leg type
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
    ee2      = x(5);
    xB       = x(6);
    thetaB   = xB*(thetaf-thetaA) + thetaA;
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
setup.rf     = rf;
setup.vf     = vf;
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
%













