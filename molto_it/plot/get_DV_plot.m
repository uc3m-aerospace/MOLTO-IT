function[DV,vf,psif,data_all] = get_DV_plot(setup)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TARGET RADIUS AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ee1    = setup.ee1;
ee2    = setup.ee2;
theta0 = setup.theta0;
thetaA = setup.thetaA;
thetaB = setup.thetaB;
r0     = setup.r0;
v0     = setup.v0;
psi0   = setup.psi0;
thetaf = setup.thetaf;
t0     = 0;
%
%------------------------------------------------------------------
% FIRST SPIRAL ARC
%------------------------------------------------------------------
%
eei = ee1;
K1i = v0^2 - 2 / r0 * ( 1- eei );
K2i = r0 * v0^2 * sin( psi0 );
%
regime0 = sign(cos(psi0));
angle   = linspace(theta0,thetaA,60);
%
[ r , v , t1 , a , e , w , at1 , psi ,alpha,DV1, ~] = propagate_spirals_DV_plot( K1i, K2i, eei, regime0, r0, theta0, t0, angle);
%
data_all = [t1'*setup.tc,r',v'*setup.vc,psi',at1'*setup.ac,alpha',angle'];
%
%------------------------------------------------------------------
% COAST ARC / KEPLERIAN ORBIT
%------------------------------------------------------------------
%
angle = linspace(thetaA,thetaB,100);
%
vk    = angle - w(end);
vv0   = thetaA - w(end);
%
r   = a(end)*( 1- e(end)^2 ) ./ ( 1 + e(end)*cos(vk) );
%
psi    = atan2 (1 + e(end) * cos( vk ), e(end) * sin( vk )) ;
%
v   = sqrt( 2./ r - 1/ abs(a(end)) );
%
n   = floor((angle-thetaA)/(2*pi));
%
% TRANFER TIME IN THE COAST ARC
%
if e(end) < (1 - eps) % ELLIPTICAL KEPLERIAN ARC
    
    EE0   =  wrapTo2Pi(atan2( sqrt(1-e(end)^2)*sin(vv0),(e(end)+cos(vv0))));
    EEf   =  wrapTo2Pi(atan2( sqrt(1-e(end)^2)*sin(vk),(e(end)+cos(vk))));
    EEf   =  EEf + 2*pi*(EEf<EE0);
    EEf   =  EEf + 2*pi*n;
    tk    =  sqrt(a(end)^3)*(EEf - e(end)*sin(EEf) - EE0 + e(end)*sin(EE0));
    
elseif ((1 + eps) < e(end)) %  HIPERBOLIC KEPLERIAN ARC
    
    FF0   =  2 * atanh( sqrt( (e(end)-1)/(e(end)+1) ) * tan(vv0/2) );
    FFf   =  2 * atanh( sqrt( (e(end)-1)/(e(end)+1) ) * tan(vk/2) );
    FFf   =  FFf + 2*pi*(FFf < FF0);
    FFf   =  FFf + 2*pi*n;
    tk    =  sqrt(abs(a(end))^3)*(e(end)*sinh(FFf) -FFf - e(end)*sinh(FF0) + FF0);
    
else % PARABOLIC KEPLERIAN ARC
    
    DD0   =  tan(vv0/2);
    DDf   =  tan(vk/2);
    DDf   =  DDf + 2*pi*(DDf < DD0);
    DDf   =  DDf + 2*pi*n;
    tk    =  sqrt(a(end)^3)*(DDf + 1/3*DDf^3 - DD0 - 1/3*DD0^3);
    
end
%
t2      = t1(end)+tk;
at2     = zeros(size(t2));
alpha   = zeros(size(t2));
%
data_all = [data_all;t2(2:end)'*setup.tc,r(2:end)',v(2:end)'*setup.vc,psi(2:end)',at2(2:end)'*setup.ac,alpha(2:end)',angle(2:end)'];
%
%
%------------------------------------------------------------------
% SECOND SPIRAL ARC
%------------------------------------------------------------------
%
if setup.type >0
    %
    eei = ee2;
    %
    K1i = (2*a(end)*eei-r(end))/(r(end)*a(end));
    K2i = sqrt( 1 + 2*e(end)*cos(vk(end)) + e(end)^2);
    %
    regime0 = sign(cos(psi(end)));
    angle   = linspace(thetaB,thetaf,100);
    %
    [  r , v , t3 , ~ , ~ , ~ , at3 , psi ,alpha,DV2, ~] = propagate_spirals_DV_plot( K1i, K2i, eei, regime0, r(end), thetaB, t0, angle);
    %
    t3 = t2(end) + t3;
    %
    data_all = [data_all;t3(2:end)'*setup.tc,r(2:end)',v(2:end)'*setup.vc,psi(2:end)',at3(2:end)'*setup.ac,alpha(2:end)',angle(2:end)'];
    %
else
    DV2 = 0;
end
%
%------------------------------------------------------------------
% COMPUTE TOTAL DV
%------------------------------------------------------------------
%
DV = abs(DV1) + abs(DV2) ;
%
psif = psi(end);
%
vf   = v(end);
%

