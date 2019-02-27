%%--------------------------------------------------------------------------
function[DV,out,flag] = flyby_spiral(planet1,planet2,et0,ToF,type,aux)
%%--------------------------------------------------------------------------
lc  = 149597870.700e03;
mu  = 132712440018e09;
tc  = sqrt(lc^3/mu);
%%--------------------------------------------------------------------------
setup.planet1      = planet1;
setup.planet2      = planet2;
setup.et0          = et0;
setup.ToF          = ToF/(tc)*(3600*24);
setup.vfb          = aux.vfb;
setup.psifb        = aux.psifb;
setup.type         = type;
setup.theta0       = aux.theta0;
setup.r0           = aux.r0;
setup.v0           = aux.v0;
setup.psi0         = aux.psi0;
setup.n            = 0;
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
setup.oe0 = oe0;
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
setup.oe = oe;
%
% Determine pseudoperiod
%
ra1 = oe0.SMA;%*(1+oe0.ECC );
p1  = 2*pi * sqrt(ra1^3);
%
ra2 = oe.SMA;%*(1+oe.ECC );
p2  = 2*pi * sqrt(ra2^3);
%amax = max(oe0.SMA,oe.SMA);
%
% if strcmp(planet1,planet2)
%     tmin = 0.1*p1;
%     tmax = 4*p1;
% elseif amax <2
%     tmin = 0.1*min(p1,p2);
%     tmax = 2*max(p1,p2);
% else
%     tmin = 0.1*min(p1,p2);
%     tmax = max(p1,p2);
% end
%
% expeted_n_min = setup.ToF/(max(p1,p2));
% expeted_n_max = 0.7*setup.ToF/(min(p1,p2));
% n_min = floor(expeted_n_min);
% n_max = floor(expeted_n_max);
% if n_min == n_max
%     setup.n = n_min;
% elseif n_min > n_max
%     setup.n = floor(expeted_n_min);
% else
%     setup.n = n_max;
% end
% %
% if strcmp(planet1,planet2)
%     setup.n = floor(expeted_n_min) ;
% end
%
%setup.n = n_min;
%
% if strcmp(planet1,'Mercury')
%         setup.n = floor(expeted_n_min) -1;
% end
%
a_med = (ra1+ra2)/2;
pmed  = 2*pi * sqrt(a_med^3);

setup.n = floor(setup.ToF/pmed);
% setup.n = 0;
% if setup.ToF > (tmin+tmax)/2
% setup.n = 1;
% end

%
tlb =  - 150/(tc)*(3600*24);
tub =  + 150/(tc)*(3600*24);
%
% if tlb > 600/(tc)*(3600*24)
%     tlb = 600/(tc)*(3600*24);
% end
% %
% if tub < 1000/(tc)*(3600*24)
%     tub = 1000/(tc)*(3600*24);
% end
% 
% if setup.ToF < 1.1*tlb
%     tg = 1.1*tlb;
% elseif setup.ToF > 0.9*tub
%     tg = 0.9*tub;
% else 
%     tg = setup.ToF;
% end
%
tg = setup.ToF;
%%--------------------------------------------------------------------------
% Solve Minimization Problem
%%--------------------------------------------------------------------------
%
flag = 0;
options = [];
options.Algorithm         = 'sqp';
options.Display           = 'iter';
%  
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'Display','off','MaxFunctionEvaluations',9000,'FunctionTolerance',1e-5,'ConstraintTolerance',1e-3);

if aux.plot == 1
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'Display','iter','MaxFunctionEvaluations',9000,'FunctionTolerance',1e-5,'ConstraintTolerance',1e-3);
end
%
try
    %
    lb = [0, tlb, 0, -1];
    ub = [1, tub,  1,  1];
    x0 = [0.5, 0, 0.01, 0 ];
    %
    if type == 2
       setup.vinff_max = aux.vinff_max; 
    end
    if type>0

        lb = [lb 0 0];
        ub = [ub 1 1];
        x0 = [x0 0.5 0.99];
    end
    %
    NONLCON = @(x)flyby_cons(x,setup);
    %
    [x,~,exitflag] = fmincon(@(x)flyby_obj(x,setup),x0,[],[],[],[],lb,ub,NONLCON,options);
catch
    %
    x = x0;
    exitflag = -1;
    %
end
%
if exitflag <=0
try
    %
    lb = [0, tlb, 0, -1];
    ub = [1, tub,  1,  1];
    x0 = [0.5, 0, 0.01, 0 ];
    %
    if type == 2
       setup.vinff_max = aux.vinff_max; 
    end
    %
    if type >0
        
        lb = [lb 0 0];
        ub = [ub 1 1];
        x0 = [x0 0.2 0.99];
    end
    %
    NONLCON = @(x)flyby_cons(x,setup);
    %
    [x,~,exitflag] = fmincon(@(x)flyby_obj(x,setup),x0,[],[],[],[],lb,ub,NONLCON,options);
catch
    %
    x = x0;
    exitflag = -1;
    %
end
%
end
%
[c,ceq]  = flyby_cons(x,setup);

if exitflag >0 && (norm(ceq) > 1e-2)
    
    try
    NONLCON = @(x)flyby_cons(x,setup);
    %
    [x,~,exitflag] = fmincon(@(x)flyby_obj(x,setup),x,[],[],[],[],lb,ub,NONLCON,options);
    catch
      exitflag  = -1;  
    end
end
%
[c,ceq]  = flyby_cons(x,setup);
%
if exitflag <= 0||(norm(ceq) > 1e-2)
    %
    flag = -1;
    DV   = NaN;
    out  = [];
    out.date1 = '0000 000 00 00:00:00';
    out.date2 = '0000 000 00 00:00:00';
    out.ToF_year =0;
    out.t = [];
    out.at = [];
    out.vrel = [];
    out.data_all = [];
    %
else
    %
    if aux.plot == 0
        %
        [DV,out] = flyby_obj(x,setup);
        %
    else
        %
        [DV,out] = flyby_plot(x,setup);
        out.date1 = cspice_et2utc(et0,'C',0);
        out.date2 = cspice_et2utc(et0 + out.ToF_year*365*24*3600,'C',0);
        %
    end
    %
end
%
end




