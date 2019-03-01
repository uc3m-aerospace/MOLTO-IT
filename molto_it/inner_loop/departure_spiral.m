%--------------------------------------------------------------------------
function[DV,out,flag] = departure_spiral(planet1,planet2,et0,ToF,type,aux)
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
%    Function that solve the low-thrust interplanetary transfer problem
%    between two given planets. The condition in planet1 is before launch
%    and is to flyby or rendezvous planet2
%
%--------------------------------------------------------------------------
tc  = aux.tc;
lc  = aux.lc;
mu  = aux.mu;
%%--------------------------------------------------------------------------
aux.planet1      = planet1;
aux.planet2      = planet2;
aux.et0          = et0;
aux.ToF          = ToF/(tc)*(3600*24);
aux.type         = type;
aux.n            = 0;
%--------------------------------------------------------------------------
% Orbital Elements for Departure planet
%--------------------------------------------------------------------------
S        = cspice_spkezr(planet1, et0, 'ECLIPJ2000', 'NONE', '0');
O        = cspice_oscelt(S, et0, mu/1e9 );
%
oe0.ECC   = O(2);
oe0.INC   = O(3);
oe0.LNODE = O(4);
oe0.ARGP  = O(5);
oe0.SMA   = O(1)/(1-oe0.ECC)*1e3/lc;
oe0.M0    = O(6);
aux.oe0   = oe0;
%--------------------------------------------------------------------------
% Orbital Elements for Arrival planet
%--------------------------------------------------------------------------
S        = cspice_spkezr(planet2, et0, 'ECLIPJ2000', 'NONE', '0');
O        = cspice_oscelt(S, et0, mu/1e9 );
%
oe.ECC   = O(2);
oe.INC   = O(3);
oe.LNODE = O(4);
oe.ARGP  = O(5);
oe.SMA   = O(1)/(1-oe.ECC)*1e3/lc;
oe.M0    = O(6);
aux.oe = oe;
%
% Guess for time
%
tlb =  - 200/(tc)*(3600*24);
tub =  + 200/(tc)*(3600*24);
%
% Guess for psi0
%
if oe.SMA > oe0.SMA
    psi0 = pi/2;
elseif oe.SMA < oe0.SMA
    psi0 = -pi/2;
else
    psi0 = 0;
end
%
%
%%--------------------------------------------------------------------------
% Solve Minimization Problem
%%--------------------------------------------------------------------------
%
%
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'Display','off','MaxFunctionEvaluations',9000,'FunctionTolerance',1e-5,'ConstraintTolerance',1e-3);
%
if aux.plot == 1
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'Display','iter','MaxFunctionEvaluations',2000,'FunctionTolerance',1e-5,'ConstraintTolerance',1e-10);
end
%
try
    %
    lb = [0, tlb, 0, -0.1, -pi];
    ub = [1, tub, 1,  0.1, pi];
    x0 = [0.3, 0, 0.1, 0, psi0];
    %
    if type == 2
        aux.vinff_max = aux.vinff_max;
    end
    %
    if type == 1 || type == 2
        
        lb = [lb 0 0];
        ub = [ub 1 1];
        x0 = [x0 0.4 0.99];
    end
    %
    NONLCON = @(x)departure_cons(x,aux);
    %
    [x,~,exitflag] = fmincon(@(x)departure_obj(x,aux),x0,[],[],[],[],lb,ub,NONLCON,options);
catch
    %
    x = x0;
    exitflag = -1;
    %
end
%
[c,ceq]  = departure_cons(x,aux);
%
if exitflag <= 0 ||(norm(ceq) > 1e-2)
    try
        %
        lb = [0, tlb, 0, -0.1, -pi];
        ub = [1, tub, 1,  0.1, pi];
        x0 = [0.5, 0, 0.2, 0, psi0 ];
        %
        if type == 2
            aux.vinff_max = aux.vinff_max;
        end
        %
        if type == 1 || type == 2
            lb = [lb 0 0];
            ub = [ub 1 1];
            x0 = [x0 0.4 0.99];
        end
        %
        NONLCON = @(x)departure_cons(x,aux);
        %
        [x,~,exitflag] = fmincon(@(x)departure_obj(x,aux),x0,[],[],[],[],lb,ub,NONLCON,options);
    catch
        %
        x = x0;
        exitflag = -1;
        %
    end
    
end
%
[c,ceq]  = departure_cons(x,aux);
%
if exitflag <= 0 ||(norm(ceq) > 1e-2)
    %
    flag = -1;
    out  = [];
    out.et = 0;
    DV   = NaN;
    %
else
    %
    flag     = 0;
    if aux.plot == 0
        [DV,out] = departure_obj(x,aux);
    else
        departure_cons(x,aux)
        [DV,out] = departure_plot(x,aux);
        [DV,out] = departure_obj(x,aux);
        out.et_fact
        cspice_et2utc(out.et,'C',0)
    end
    %
end
%
end




