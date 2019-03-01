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
tlb =  - 200/(tc)*(3600*24);
tub =  + 200/(tc)*(3600*24);
%
%%--------------------------------------------------------------------------
% Solve Minimization Problem
%%--------------------------------------------------------------------------
%
options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'Display','off','MaxFunctionEvaluations',9000,'FunctionTolerance',1e-5,'ConstraintTolerance',1e-3);
%
if aux.plot == 1
    options = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',1500,'Display','iter','MaxFunctionEvaluations',2000,'FunctionTolerance',1e-5,'ConstraintTolerance',1e-10);
end
%
try
    %
    lb = [0, tlb, 0, -1];
    ub = [1, tub,  1,  1];
    x0 = [0.4, 0, 0.01, 0 ];
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
        x0 = [0.6, 0, 0.01, 0 ];
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
if exitflag <= 0 ||(norm(ceq) > 1e-2)
    %
    flag = -1;
    DV   = NaN;
    out  = [];
    %
else
    %
    flag = 0;
    %
    if aux.plot == 0
        [DV,out] = flyby_obj(x,setup);
    else
        [c,ceq]  = flyby_cons(x,setup)
        [DV,out] = flyby_plot(x,setup);
        [DV,out] = flyby_obj(x,setup);
        %    expeted_n_min
        %    expeted_n_max
    end
    %
end
%
end




