%------------------------------------------------------------------
% COAST ARC / KEPLERIAN ORBIT 
%------------------------------------------------------------------
%
K1  = v^2 - 2 / r * ( 1- ee1 );
K2  = r * v^2 * sin( psi);
% 
a   = (r./(2*ee1-K1*r));
e   = sqrt( K2^2+1 - 2*K2*sin(psi) );
w   = atan2( sin(theta) + K2 * cos( theta + psi ), cos(theta) - K2 * sin(theta + psi) ) + pi;
%
vk  = thetaB - w(end);
vv0 = thetaA - w(end);
%
n   = floor((thetaB-thetaA)/(2*pi));
%
r   = abs(a(end)*( 1- e(end)^2 )) ./ ( 1 + e(end)*cos(vk) );
%
v   = sqrt( 2/ r - 1/ a(end) );
%
psi = atan2 ( 1 + e(end) * cos( vk ), e(end) * sin( vk ) ) ;
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