function [Obj, Cons, outdata] = fitness_plot(x, setup, ~)
%
%%-------------------------------------------------------------------------
% AUXILIAR PARAMETERS
%%-------------------------------------------------------------------------
%
lc = 149597870.700e03;
mu = 132712440018e09;
tc = sqrt(lc^3/mu);
vc = lc / tc;
ac = lc / tc^2;
%
%%-------------------------------------------------------------------------
% DECODE INPUT VECTOR
%%-------------------------------------------------------------------------
%
ind = setup.ind;
%
%%-------------------------------------------------------------------------
% GET FLYBY SEQUENCE
%%-------------------------------------------------------------------------
%
fbb = x(ind.fbb(1):ind.fbb(2));
planets = setup.planets;
n_fb = length(fbb);
planet_dep = setup.planet_dep;
planet_arr = setup.planet_arr;
planet_flyby = cell(1, n_fb);
planet_available = setup.planet_available;
%
% Order Flyby Sequence
%
fbb_new = fbb;
j = 0;
k = 1;
for i = 1:n_fb
    %
    if fbb(i) > planet_available
        fbb_new(end-j) = fbb(i);
        j = j + 1;
    else
        fbb_new(k) = fbb(i);
        k = k + 1;
    end
    %
end
%
%
x(ind.fbb(1):ind.fbb(2)) = fbb_new;
%
%
fbb_new2 = [fbb_new(fbb_new > planet_available), fbb_new(fbb_new <= planet_available)];
%
fbb = fbb_new2;
%
%
for i = 1:n_fb
    %
    if fbb(i) > planet_available
        planet_flyby(i) = {'Null'};
    else
        planet_flyby(i) = planets(ind.planets(fbb(i)));
    end
    %
end
%
%
planet_seq = [planet_dep, planet_flyby, planet_arr];
%
%%--------------------------------------------------------------------------
% GET THE TIME OF FLIGHTS
%%--------------------------------------------------------------------------
%
ToFg = x(ind.ToF(1):ind.ToF(2));
%
%%--------------------------------------------------------------------------
% GET THE INITIAL TIME
%%--------------------------------------------------------------------------
%
et0_min = cspice_str2et(setup.t0_min);
et0_max = cspice_str2et(setup.t0_max);
et0_factor = x(ind.t0);
et0 = et0_min + (et0_max - et0_min) * et0_factor;
data.initial_date = et0_min;
data.final_date = et0_max;
data.xind = x(ind.t0);
data.vinf0 = setup.vinf0;
data.plot = setup.plot;
transfer_type = setup.type;

if isfield(setup, 'vinff_max')
    data.vinff_max = setup.vinff_max * 1e3 / vc;
end
%
%%--------------------------------------------------------------------------
% INITIALIZE VARIABLES
%%--------------------------------------------------------------------------
%
et = zeros(1, n_fb+2);
r = zeros(1, n_fb+2);
v = zeros(1, n_fb+2);
theta = zeros(1, n_fb+2);
psi       = zeros(1,n_fb+2);
type      = zeros(1,n_fb+1);
type(end) = transfer_type;
et(1) = et0;
%
%%--------------------------------------------------------------------------
% PROPAGATE SPIRALS
%%--------------------------------------------------------------------------
%
DV = zeros(n_fb+1, 1);
ToF = zeros(n_fb+2, 1);
flag  = 0;
%
j     = 0;
t = 0;
flagT = n_fb + 1;
outdata.time = [];
outdata.at = [];
%
for i = 1:n_fb + 1
    %
    if strcmp(planet_seq(i+1), 'Null')
        %
        et(i+1) = et(i);
        r(i+1) = r(i);
        v(i+1) = v(i);
        theta(i+1) = theta(i);
        psi(i+1) = psi(i);
        ToF(i+1) = ToF(i);
        ToF(i) = 0;
        flagT = flagT - 1;
        %
    else
        %
        %%--------------------------------------------------------------------------
        % Compute Leg Spiral
        %%--------------------------------------------------------------------------
        %
        if j == 0
            %
            % Departure
            %
            [DV(i), out, flag] = departure_spiral(planet_seq(1), planet_seq(i+1), et(1), ToFg(t+1), type(i), data);
            %
            et(i) = out.et;
            if flag == 0
                x(ind.t0) = out.et_fact;
            end
            %
            outdata.date1(i, :) = out.date1;
            outdata.date2(i, :) = out.date2;
            outdata.all(i).data = out.data_all;
            outdata.ToF(i) = out.ToF_year * 365;
            outdata.time = [outdata.time, out.t * tc / (3600 * 24 * 365)];
            outdata.at = [outdata.at, out.at * ac];
            %
        else
            %
            % From Flyby
            %
            [DV(i), out, flag] = flyby_spiral(planet_seq(i), planet_seq(i+1), et(i), ToFg(t+1), type(i), data);
            %
            outdata.date1(i, :) = out.date1;
            outdata.date2(i, :) = out.date2;
            outdata.ToF(i) = out.ToF_year * 365;
            outdata.all(i).data = out.data_all;
            outdata.time = [outdata.time, outdata.time(end) + out.t * tc / (3600 * 24 * 365)];
            outdata.at = [outdata.at, out.at * ac];
            %
            outdata.vrel = out.vrel;
            %
        end
        %
        j = j + 1;
        %
        %%--------------------------------------------------------------------------
        % Determine if a feasible trajectory was found
        %%--------------------------------------------------------------------------
        %
        if flag == -1 || DV(i) > 500
            break
        end
        %
        ToF(i+1) = out.ToF;
        x(ind.ToF(1)+t) = ToF(i+1) * tc / (3600 * 24);
        data.vfb = out.vp;
        data.psifb = out.psip;
        data.theta0 = out.thetasp;
        data.r0 = out.rsp;
        data.v0 = out.vsp;
        data.psi0 = out.psisp;
        et(i+1) = ToF(i+1) * tc + et(i);
        %
        flagT = flagT - 1;
        %
    end
    %
    t = t + 1;
    %
end
%
%--------------------------------------------------------------------------
% COMPUTE TOTAL DV FOR THE SPIRAL
%%--------------------------------------------------------------------------
%
if flag == 0
    %
    DV_main = (sum(DV)) * vc / 1000;
    DV_total = 1 - exp(-(DV_main / 4.19)/(9.81));
    Time = (sum(ToF)) * tc / (3600 * 24 * 365);
    %
else
    %
    DV_total = 1;
    Time = 10;
    %
end
%
%if Time < 2
%    flagT = 1;
%end
%
%--------------------------------------------------------------------------
% COMPUTE OBJECTIVE FUNCTION
%--------------------------------------------------------------------------
%
Obj = [Time, DV_total]
%
Cons = flagT
%--------------------------------------------------------------------------
