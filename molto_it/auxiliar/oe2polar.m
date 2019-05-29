function [rc, theta, vc, psi] = oe2polar(oe, deltaT)

Mt = deltaT * sqrt(1/oe.SMA^3) + oe.M0;

%Calculate eccentric anomaly using Newton's method
E = Mt;
F = E - oe.ECC * sin(E) - Mt;
j = 0;
maxIter = 1000;
delta = 0.0000000001;
%
while (abs(F) > delta && j < maxIter)
    E = E - F / (1 - oe.ECC * cos(E));
    F = E - oe.ECC * sin(E) - Mt;
    j = j + 1;
end

if j == maxIter

    stop

end
%
%nu = 2 * atan2 (sqrt (1 + oe.ECC) * sin (E / 2), sqrt (1 - oe.ECC) * cos (E / 2)); %True anomaly
nu = atan2(sin(E)*sqrt(1-oe.ECC^2), cos(E)-oe.ECC);
%
h = sqrt(oe.SMA*(1 - oe.ECC^2));
p = oe.SMA * (1 - oe.ECC^2);
%
rc = oe.SMA * (1 - oe.ECC * cos(E)); %//Distance to central body
%
rx = rc * (cos(oe.LNODE) * cos(oe.ARGP+nu) - sin(oe.LNODE) * sin(oe.ARGP+nu) * cos(oe.INC));
ry = rc * (sin(oe.LNODE) * cos(oe.ARGP+nu) + cos(oe.LNODE) * sin(oe.ARGP+nu) * cos(oe.INC));
rz = rc * (sin(oe.INC)) * sin(oe.ARGP+nu);
%
r = [rx, ry, rz];
%

theta = atan2(ry, rx);
%
vx = rx * h * oe.ECC / (rc * p) * sin(nu) - h / rc * (cos(oe.LNODE) * sin(oe.ARGP+nu) + sin(oe.LNODE) * cos(oe.ARGP+nu) * cos(oe.INC));
vy = ry * h * oe.ECC / (rc * p) * sin(nu) - h / rc * (sin(oe.LNODE) * sin(oe.ARGP+nu) - cos(oe.LNODE) * cos(oe.ARGP+nu) * cos(oe.INC));
vz = rz * h * oe.ECC / (rc * p) * sin(nu) + h / rc * sin(oe.INC) * cos(oe.ARGP+nu);
%
v = [vx, vy, vz];
%
vc = norm(v);
%
psi = (atan2(norm(cross(r, v)), +dot(r, v)));
%


%     a = oe.SMA;
%     i = oe.INC;
%     e = oe.ECC;
%     omega = oe.LNODE;
%     w     = oe.ARGP;
%     mu = 1;
%     h = sqrt(2*a*mu*(1/(1 + e) + 1/(1 - e))^(-1));

%position and velocity vectors in perifocal reference frame
% r = (h^2/mu)*(1/(1 + e*cos(nu)))*(cos(nu)*[1;0;0] + sin(nu)*[0;1;0]); %[km]
% v = (mu/h)*(-sin(nu)*[1;0;0] + (e + cos(nu))*[0;1;0]); %[km/s]
%
% %1st Rotation matrix about the z-axis through the ascending node angle(omega)
% R3_omega = [cos(omega) sin(omega) 0 ; -sin(omega) cos(omega) 0 ; 0 0 1];
%
% % 2nd Rotation matrix about the x-axis through the inclination angle(i)
% R1_i = [1 0 0 ; 0 cos(i) sin(i) ; 0 -sin(i) cos(i)];
%
% %3rd Rotation matrix about the z-axis through the argument of the perigee(w)
% R3_w = [cos(w) sin(w) 0 ; -sin(w) cos(w) 0 ; 0 0 1];
%
% %Total transformation matrix to pass from geocentrical equatorial to
% %periphocal frame, this is why later on the transpose will be applied.
% TTm = R3_w*R1_i*R3_omega;
%
% %position and velocity vectors in geocentric equatorial reference frame,
% %applying the transpose of the transformation matrix previously computed.
% r1_vec = TTm.'*r; %[km]
% v1_vec = TTm.'*v; %[km/s]


end