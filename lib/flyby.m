function [vf,psif, hpm, v_inf] = flyby(v0,psi0,vfb,psifb,mu_cb,fb_b,d_factor)
%--------------------------------------------------------------------------
% PERFORMED FLYBY IN MARS
%--------------------------------------------------------------------------
%
lc  = 149597870.700e03;
mu  = 132712440018e09;
tc  = sqrt(lc^3/mu);
vc  = lc/tc;
%
%id      = cspice_bodn2c(fb_b);

id      = str2double(fb_b);
id      = id*100+99;
R_flyby = mean( cspice_bodvcd(id,'RADII',3) )*1e3/lc;
mum     = cspice_bodvcd(id,'GM',3)*1e9/mu_cb;
hpmin   = 200*1e3/lc; 
%--------------------------------------------------------------------------
% project velocity onto the orbitalframe
%--------------------------------------------------------------------------
%
vrr1         =  v0*cos(psi0);
voo1         =  v0*sin(psi0);
%
vrr_planet   =  vfb*cos(psifb);
voo_planet   =  vfb*sin(psifb);
%
v_inf_rr1    =  vrr1-vrr_planet;
v_inf_oo1    =  voo1-voo_planet;
v_inf        =  sqrt ( v_inf_rr1^2+v_inf_oo1^2 );
%v_inf*vc/1e3
%
angle_v_inf1 =  atan2(v_inf_oo1,v_inf_rr1); 
%
delta_max    =  2*asin( 1 ./(  1 + v_inf^2 *(R_flyby + hpmin)/mum  ) ) ;
%
delta        = d_factor*delta_max;
%
hpm          = (  ( 1/ sin(abs(delta/2))-1 )*mum/v_inf^2 - R_flyby )*lc/1e3;
%
angle_v_inf2 =  angle_v_inf1 + delta;
%
v_inf_rr2    =  v_inf*cos(angle_v_inf2);
v_inf_oo2    =  v_inf*sin(angle_v_inf2);
%
vrr2         =  vrr_planet + v_inf_rr2;
voo2         =  voo_planet + v_inf_oo2;
%
psif         =  atan2(voo2,vrr2);
vf           =  sqrt( vrr2^2 + voo2^2 );
%



