function [Obj,Cons] = fitness_nsga2(x,setup,~)
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
%   Function that computes the Flight time and opropellant mass fraction
%   given an input gen x
%
%%-------------------------------------------------------------------------
% DEFINE AUXILIAR PARAMETERS
%%-------------------------------------------------------------------------
%
lc  = 149597870.700e03;
mu  = 132712440018e09;
tc  = sqrt(lc^3/mu);
vc  = lc/tc;
ac  = vc/tc;
%
data.lc = lc;
data.mu = mu;
data.tc = tc;
data.vc = vc;
data.ac = ac;
%
%%-------------------------------------------------------------------------
% DECODE INPUT VECTOR
%%-------------------------------------------------------------------------
%
ind  = setup.ind;
%
%%-------------------------------------------------------------------------
% GET FLYBY SEQUENCE
%%-------------------------------------------------------------------------
%
fbb          = x(ind.fbb(1) : ind.fbb(2));
planets      = setup.planets;
n_fb         = length(fbb);
planet_dep   = setup.planet_dep ;
planet_arr   = setup.planet_arr ;
planet_flyby = cell(1,n_fb);
planet_avai  = setup.planet_avai;
nfb_r = n_fb;
%
for i = 1:n_fb
    %
    if fbb(i) > planet_avai
        planet_flyby(i) = {'Null'};
        nfb_r = nfb_r-1;
    else
        planet_flyby(i) = planets((fbb(i)));
    end
    %
end
%
%
planet_seq = [planet_dep,planet_flyby,planet_arr];
%
%%--------------------------------------------------------------------------
% GET TIME OF FLIGHTS and Revolutions Number
%%--------------------------------------------------------------------------
%
ToFg = x(ind.ToF(1):ind.ToF(2));
rev  = x(ind.rev(1):ind.rev(2));
%
%%--------------------------------------------------------------------------
% GET INITIAL TIME
%%--------------------------------------------------------------------------
%
et0_min           = cspice_str2et(setup.t0_min);
et0_max           = cspice_str2et(setup.t0_max);
et0_factor        = x(ind.t0);
et0               = et0_min + (et0_max - et0_min)*et0_factor;
data.initial_date = et0_min;
data.final_date   = et0_max;
data.xind         = x(ind.t0);
data.vinf0        = setup.vinf0*1e3/vc;
data.plot         = setup.plot;
transfer_type     = setup.type;
%
%%--------------------------------------------------------------------------
% GET IhYPERBOLIC AVAILABLE VELOCITY (IF DEFINED)
%%--------------------------------------------------------------------------
%
if isfield(setup,'vinff_max')
    data.vinff_max      = setup.vinff_max*1e3/vc;
end
%
%%--------------------------------------------------------------------------
% INITIALIZE VARIABLES
%%--------------------------------------------------------------------------
%
et        = zeros(1,n_fb+2);
r         = zeros(1,n_fb+2);
v         = zeros(1,n_fb+2);
theta     = zeros(1,n_fb+2);
psi       = zeros(1,n_fb+2);
type      = zeros(1,n_fb+1);
type(end) = transfer_type;
et(1)     = et0;
%
%%--------------------------------------------------------------------------
% PROPAGATE SPIRALS
%%--------------------------------------------------------------------------
%
DV    = zeros(n_fb+1,1);
ToF   = zeros(n_fb+2,1);
flag  = 0;
%
j     = 0;
t     = 0;
flagT = n_fb+1;
%
for i = 1 : n_fb + 1
    %
    data.n = rev(i);
    
    %%--------------------------------------------------------------------------
    % Null Flyby
    %%--------------------------------------------------------------------------
    %
    if strcmp(planet_seq(i+1),'Null')
        %
        planet_seq(i+1) = planet_seq(i);
        et(i+1)    = et(i);
        r(i+1)     = r(i);
        v(i+1)     = v(i);
        theta(i+1) = theta(i);
        psi(i+1)   = psi(i);
        ToF(i+1)     = 0;
        flagT      = flagT -1;
        %
    else
        %
        %%--------------------------------------------------------------------------
        % Compute Leg Spiral
        %%--------------------------------------------------------------------------
        %
        if  j ==0
            %
            % Departure
            %
            [DV(i),out,flag] = departure_spiral(planet_seq(1),planet_seq(i+1),et(1),ToFg(t+1),type(i),data);
            %
            % Update Initial Launch date
            %
            et(i) = out.et;
            %
        else
            %
            % From Flyby
            %
            [DV(i),out,flag] = flyby_spiral(planet_seq(i),planet_seq(i+1),et(i),ToFg(t+1),type(i),data);
            %
        end
        %
        j = j+1;
        %
        %%--------------------------------------------------------------------------
        % Break loop is an unfeasible trajectory is found
        %%--------------------------------------------------------------------------
        %
        if flag == -1 || DV(i) > 500
            break
        end
        %
        % Update initial values for the following leg
        %
        ToF(i+1)        = out.ToF;
        data.vfb        = out.vp;
        data.psifb      = out.psip;
        data.theta0     = out.thetasp;
        data.r0         = out.rsp;
        data.v0         = out.vsp;
        data.psi0       = out.psisp;
        et(i+1)         = ToF(i+1)*tc + et(i);
        %
        flagT = flagT - 1;
        %
    end
    t = t+1;
    
end
%
%--------------------------------------------------------------------------
% COMPUTE TOTAL DV FOR THE SPIRAL
%%--------------------------------------------------------------------------
%
if flag == 0
    %
    % Feasible Trajectory
    %
    DV_main  =  (sum(DV)) *  vc / 1000;
    DV_total =  1 - exp ( - ( DV_main /3 )  / (9.81) );
    Time = ( sum(ToF) )*tc/(3600*24*365) ;
    %
    % Penalty for high DV trajectories
    %
    if DV_total > 0.5
        flagT = 0.5;
    end
else
    %
    % Unfeasible Trajectory
    %
    DV_total = 1;
    Time = 10;
    %
end
%
%--------------------------------------------------------------------------
% COMPUTE OBJECTIVE FUNCTION
%--------------------------------------------------------------------------
%
Obj  = [Time, DV_total, -nfb_r]
%
Cons = flagT
%
%--------------------------------------------------------------------------









