function T_available = thrust_model(r)

P0  = 30;   %KW
m0  = 3000; %Kg
P   = (P0./ r.^2-1)*0.85/2;
P   = P.*(P>0);
%
Pa  = 7;
P   =  Pa*(P>Pa) + P.*(P<=Pa) ;
%
% Thrust
%
T_available =   2*(1.19388817E-02 + 1.60989424E-02*P + 1.14181412E-02*P.^2 -2.04053417E-03*P.^3 + 1.19388817E-04*P.^4)/m0.* (P>0);