function [r, v, t, a, e, w, at, psi, alpha, DV, flag] = propagate_spirals_DV_plot(K1, K2, ee, regime0, r0, theta0, t0, theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              CONTROLLED GENERALIZED LOGARITHMIC SPIRALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that computes the state vector [r(position), v(velocity), t(time)]
% and the classical orbital elements at specific angular positions for a
% given spiral defined by its control parameters K1, K2, and ee, given a
% initial position (r0,theta0)
% INPUTS ::
%        K2     :: Second control parameter
%        K1     :: First control parameter
%        r0     :: Initial Distance

%%
%        theta0 :: Initial orientation angle
%        t0     :: Initial flight time
%        ee     :: Third Control parameter
%        regime0:: Initial spiral regime (raising[+1],lowering[-1])
%        theta  :: Orientations at which the state is to be computed (scalar or vector)
% OUTPUTS ::
%        r      :: Distance at theta
%        v      :: velocity at theta
%        t      :: time at theta
%        a      :: semimajor axis      (OPTIONAL)
%        e      :: eccentricity        (OPTIONAL)
%        w      :: argument of perigee (OPTIONAL)
%        at     :: Thrust acceleration (OPTIONAL)
%        psi    :: flight path angle   (OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations are based on Roa et al. (Nov.15)
% 'Controlled generalized logarithmic spirals for low-thrust mission desing'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE 1 ! ALL VARIBLES SHOULD NON-DIMENSIONAL GIVEN THAT mu=0
% NOTE 2 ! ee = 1/2 correspond to the tangential steering case
% NOTE 3 ! We consider that when phi0 = pi/2 we are in raising regime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check for consistent values
%
flag  = 1;
PI = 3.14159265358979D0;
%
% Initial flight path angle
%
cpsi0 = regime0 * (sqrt((2 * (1 - ee) + K1 * r0)^2-K2^2) / (2 * (1 - ee) + K1 * r0));
spsi0 = K2 / (2 * (1 - ee) + K1 * r0);
%
psi0 = atan2(real(spsi0), real(cpsi0));
%
if flag == 0

elseif K1 < -abs(eps) % Elliptical spiral
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE PARTICLE NEVER SCAPE (rmax) when propagated backward or forward
    % the particle reach the origin of the central body
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rmax = (2 * (1 - ee) - K2) / (-K1);
    %
    ell = sqrt(4*(1 - ee)^2-K2^2);
    %
    thetam = theta0 + regime0 * K2 / ell * abs(acosh(rmax/r0-2*(1 - ee)/K2*(1 - rmax / r0)));
    regime = -1 * (theta > thetam) + 1 * (theta <= thetam);
    beta = ell / K2 * (theta - thetam);
    %
    % Trajectory
    %
    r = rmax * (2 * (1 - ee) + K2) ./ (2 * (1 - ee) + K2 * cosh(beta));
    v = sqrt(K1+2./r*(1 - ee));
    %
    cpsi = regime .* (sqrt((2 * (1 - ee) + K1 * r).^2-K2^2) ./ (2 * (1 - ee) + K1 * r));
    spsi = K2 ./ (2 * (1 - ee) + K1 * r);
    psi = atan2(spsi, cpsi);
    %
    % Time of flight
    %
    vm = sqrt(K2/rmax);
    kk = sqrt(-K1*rmax/(4 * (1 - ee)));
    kp = sqrt(1-kk^2);
    nn = K1 * rmax / (2 * K2);
    %
    % Elliptic integrals
    %
    v0 = v(1);
    psi0 = psi(1);
    r0 = r(1);
    %
    phi0 = real(asin(vm/v0*sqrt(2/(1 + sin(psi0)))));
    %
    [~, E01] = ELIT(kk, 90.0);
    [~, E02] = ELIT(kk, phi0*180/PI);
    %
    Pi01 = ELIT3(90.0, kk, nn);
    Pi02 = ELIT3(phi0*180/PI, kk, nn);
    %
    E0 = E02 - E01;
    Pi0 = Pi02 - Pi01;
    %
    phi = real(asin(vm./v.*sqrt(2./(1 + sin(psi)))));
    %
    [~, E03] = ELIT(kk, phi*180/PI);
    Pi03 = ELIT3(phi*180/PI, kk, nn);
    %
    E = E03 - E01;
    Pi = Pi03 - Pi01;
    %
    K4 = -regime0 .* (r0 * v0 / K1 * sqrt((1 - sin(psi0))/(1 + sin(psi0) ...
        )) + (2 * (2 * (1 - ee) * kp^2 * Pi0 - K2 * E0) * sqrt(1-ee)) / ((-K1)^(1.5) * sqrt(K2)));
    %
    t = K4 + regime .* (r .* v / K1 .* sqrt((1 - sin(psi))./(1 + sin(psi) ...
        )) + (2 * (2 * (1 - ee) * kp^2 * Pi - K2 * E) * sqrt(1-ee)) / ((-K1)^(1.5) * sqrt(K2)));
    %
    %
elseif K1 > abs(eps) % Hiperbolic  Spiral
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE PARTICLE REACHES INFINITY WITH A NON-ZERO VELOCITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if K2 < 2 * (1 - ee)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Type I Hiperbolic Spiral ( It has an asymthotic value )
        % Raising regime[ reaches infinity ], Lowering regime[ reaches origin ]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ell = sqrt(4*(1 - ee)^2-K2^2);
        %
        zeta = 2 * (1 - ee) + ell;
        thetaAs = theta0 + regime0 * K2 / ell * log(K2*(zeta - ell - K2 * sin(psi0) + ell * abs(cos(psi0)))/ ...
            (r0 * K1 * zeta * sin(psi0)));
        %
        beta = regime0 * ell / K2 * (thetaAs - theta);
        %
        r = zeta * ell^2 ./ (K1 * (sinh(beta / 2) .* (4 * zeta * (1 - ee) * sinh(beta / 2) + ...
            (zeta^2 - K2^2) * cosh(beta / 2))));
        v = sqrt(K1+2./r*(1 - ee));
        %
        cpsi = regime0 * (sqrt((2 * (1 - ee) + K1 * r).^2-K2^2) ./ (2 * (1 - ee) + K1 * r));
        spsi = K2 ./ (2 * (1 - ee) + K1 * r);
        %
        psi = atan2(spsi, cpsi);
        %
        %
        % Time of Flight
        %
        %
        kk = 1 / 2 * sqrt((2 * (1 - ee) + K2)/(1 - ee));
        nn = (2 * (1 - ee) + K2) / (2 * K2);
        %
        v0 = v(1);
        psi0 = psi(1);
        r0 = r(1);
        %
        % Ellliptic integrals
        %
        phi0 = asin(sqrt(K1*r0*sin(psi0)/(nn * K2 * (1 - sin(psi0)))));
        phi = asin(sqrt(K1*r.*sin(psi)./(nn * K2 * (1 - sin(psi)))));
        %
        % Do only if any asymptote is reached
        %
        if flag == 1

            [~, E0] = ELIT(kk, phi0*180/PI);
            Pi0 = ELIT3(phi0*180/PI, kk, nn);

            [~, E] = ELIT(kk, phi*180/PI);
            Pi = ELIT3(phi*180/PI, kk, nn);

            K4 = -regime0 * (r0 * v0 / K1 * sqrt((1 + sin(psi0))/(1 - sin(psi0))) - ...
                2 * (E0 - (1 - nn) * Pi0) * sqrt(K2*(1 - ee)) / K1^(3 / 2));

            t = K4 + regime0 * (r .* v / K1 .* sqrt((1 + sin(psi))./(1 - sin(psi))) - ...
                2 * (E - (1 - nn) * Pi) .* sqrt(K2*(1 - ee)) / K1^(3 / 2));

        end


    elseif K2 > 2 * (1 - ee)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Type II Hiperbolic Spiral
        % There is a transition between raising and lowering regime [rmin]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rmin = (K2 - 2 * (1 - ee)) / (K1);
        %
        ell = sqrt(K2^2-4*(1 - ee)^2);
        %
        thetam = theta0 - regime0 * K2 / ell * (pi / 2 + atan((2 * (1 - ee) - K2 * sin(psi0))/ ...
            (ell * abs(cos(psi0)))));
        %
        beta = ell / K2 * (theta - thetam);
        %
        r = rmin * ((2 * (1 - ee) + K2) ./ (2 * (1 - ee) + K2 * cos(beta)));
        v = sqrt(K1+2./r*(1 - ee));
        %
        regime = -1 * (theta < thetam) + 1 * (theta >= thetam);
        %
        cpsi = regime .* (sqrt((2 * (1 - ee) + K1 * r).^2-K2^2) ./ (2 * (1 - ee) + K1 * r));
        spsi = K2 ./ (2 * (1 - ee) + K1 * r);
        %
        psi = atan2(spsi, cpsi);
        %
        % Time of Flight
        %
        kk = 2 * sqrt(1-ee) / sqrt(K2+2*(1 - ee));
        nn = 2 * (1 - ee) / K2;
        v0 = v(1);
        psi0 = psi(1);
        r0 = r(1);
        %
        % Ellliptic integrals
        %
        phi0 = real(asin(sqrt(2*(1 - ee)*(1 - sin(psi0)))/(kk * sqrt(K2-2*(1 - ee)*sin(psi0)))));
        phi = real(asin(sqrt(2*(1 - ee)*(1 - sin(psi)))./(kk * sqrt(K2-2*(1 - ee)*sin(psi)))));

        if flag == 1
            %
            [F0, E0] = ELIT(kk, phi0*180/PI);
            Pi0 = ELIT3(phi0*180/PI, kk, nn);
            %
            [F, E] = ELIT(kk, phi*180/PI);
            Pi = ELIT3(phi*180/PI, kk, nn);
            %
            K4 = +regime0 * (((K2 + 2 * (1 - ee)) * K2 * E0 - K1 * rmin * (K2 * F0 + 2 * (1 - ee) * Pi0)) / ...
                (K1 * sqrt(K1*K2*(K2 + 2 * (1 - ee)))) + 2 * (1 - ee) / (K1^(3 / 2)) * asinh(sqrt( ...
                2*K1*r0*(r0 * v0^2 - K2))/(2 * sqrt(K2*r0*v0^2+(r0 * v0^2 - K2)*(1 - ee)))) + ...
                -v0 / K1^2 * sqrt(r0^2*v0^4-K2^2));
            %
            t = K4 - regime .* (((K2 + 2 .* (1 - ee)) .* K2 .* E - K1 .* rmin .* (K2 .* F + 2 .* (1 - ee) .* Pi)) ./ ...
                (K1 .* sqrt(K1.*K2.*(K2 + 2 .* (1 - ee)))) + 2 .* (1 - ee) / (K1^(3 / 2)) .* asinh(sqrt( ...
                2.*K1.*r.*(r .* v.^2 - K2))./(2 .* sqrt(K2.*r.*v.^2+(r .* v.^2 - K2).*(1 - ee)))) + ...
                -v ./ K1^2 .* sqrt(r.^2.*v.^4-K2^2));

        end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Limit case Hiperbolyc spiral
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thetaAs = theta0 - regime0 * (1 - sqrt(1+4*(1 - ee)/(K1 * r0)));
        beta = theta - thetaAs;
        %
        r = 4 * (1 - ee) ./ (K1 * beta .* (beta - regime0 * 2));
        v = sqrt(K1+2./r*(1 - ee));
        %
        cpsi = regime0 .* (sqrt((2 * (1 - ee) + K1 * r).^2-K2^2) ./ (2 * (1 - ee) + K1 * r));
        spsi = K2 ./ (2 * (1 - ee) + K1 * r);
        %
        psi = atan2(spsi, cpsi);
        %
        %
        % Time of flight
        %
        v0 = v(1);
        r0 = r(1);
        %
        H0 = v0 * sqrt(r0*(r0 * v0^2 + 2 * (1 - ee)));
        Hf = v .* sqrt(r.*(r .* v.^2 + 2 * (1 - ee)));
        t = regime0 * 1 / K1^(3 / 2) * (Hf - H0 + (1 - ee) * log((r0 * v0^2 + 1 - ee + H0)/(r * v^2 + 1 - ee + Hf)));
        %
    end

elseif abs(K1) < abs(eps) % Parabolic spiral
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE FLIGHT PATH ANGLE REMAINS CONSTANT
    % THE PARTICLE REACHES INFINITY WITH A ZERO VELOCITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    ell = sqrt(4*(1 - ee)^2-K2^2);
    %
    r = r0 * exp(regime0*ell*(theta - theta0)/K2);
    v = sqrt(K1+2./r*(1 - ee));
    psi = psi0 * ones(size(r));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the corresponding orbital elements if desired
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
e = sqrt(K2^2+1-2*K2*sin(psi));
a = (r ./ (2 * ee - K1 * r));
w = atan2(sin(theta)+K2*cos(theta+psi), cos(theta)-K2*sin(theta+psi)) + pi;
alpha = atan((1 - ee)*sin(2*psi)./(3 * ee - 1 + (1 - ee) * cos(2*psi)));
%
% Obtain the corresponding DV
%
if (abs(theta(end)-theta(1)) < 1e-12)
    % Avoid errors with zero arcs duration
    at = zeros(size(r));
    DV = 0;
    %
else
    %
    at = 1 ./ r.^2 .* sqrt(ee^2*cos(psi).^2+(1 - 2 * ee)^2*sin(psi).^2);
    integran = @(thetai, theta, at, v, psi, r) interp1(theta, at.*1./(v .* cos(psi)).*r./tan(psi), thetai, 'linear', 'extrap');
    DV = integral(@(thetai)integran(thetai, theta, at, v, psi, r), theta(1), theta(end));
    %
    if DV == 0
        flag = -1;

    end
end


end

function EL3 = ELIT3(PHI, HK, C)
%
%       =========================================================
%       Purpose: Compute the elliptic integral of the third kind
%                using Gauss-Legendre quadrature
%       Input :  Phi --- Argument ( in degrees )
%                 k  --- Modulus   ( 0 � k � 1.0 )
%                 c  --- Parameter ( 0 � c � 1.0 )
%       Output:  EL3 --- Value of the elliptic integral of the
%                        third kind
%       =========================================================
%
T = [.9931285991850949, .9639719272779138, ...
    .9122344282513259, .8391169718222188, ...
    .7463319064601508, .6360536807265150, ...
    .5108670019508271, .3737060887154195, ...
    .2277858511416451, .7652652113349734E-1];

W = [.1761400713915212E-1, .4060142980038694E-1, ...
    .6267204833410907E-1, .8327674157670475E-1, ...
    .1019301198172404, .1181945319615184, ...
    .1316886384491766, .1420961093183820, ...
    .1491729864726037, .1527533871307258];
%
LB1 = isequal(HK, 1.0E0) && (abs(PHI-90.0) <= 1.0E-8);
LB2 = isequal(C, 1.0E0) && (abs(PHI-90.0) <= 1.0E-8);
%
if (LB1 || LB2)
    EL3 = NaN;
    return
end
%
C1 = 0.87266462599716E-2 * PHI;
C2 = C1;
EL3 = 0.0E0;

for I = 1:10
    %
    C0 = C2 * T(I);
    T1 = C1 + C0;
    T2 = C1 - C0;
    %
    F1 = 1.0E0 ./ ((1.0E0 - C * sin(T1) .* sin(T1)) .* ...
        sqrt(1.0E0-HK*HK*sin(T1).*sin(T1)));
    %
    F2 = 1.0E0 ./ ((1.0E0 - C * sin(T2) .* sin(T2)) .* ...
        sqrt(1.0E0-HK*HK*sin(T2).*sin(T2)));
    %
    EL3 = EL3 + W(I) * (F1 + F2);
    %
end

EL3 = C1 .* EL3;

end

function [FE, EE] = ELIT(HK, PHI)
%
%       ==================================================
%       Purpose: Compute complete and incomplete elliptic
%                integrals F(k,phi) and E(k,phi)
%       Input  : HK  --- Modulus k ( 0 � k � 1 )
%                Phi --- Argument ( in degrees )
%       Output : FE  --- F(k,phi)
%                EE  --- E(k,phi)
%       ==================================================
%

G = 0.0;
PI = 3.14159265358979;
A0 = 1.0;
B0 = sqrt(1.0-HK*HK);
D0 = (PI / 180.0) * PHI;
R = HK * HK;
D = 0.0;

if (isequal(HK, 1.0) && isequal(PHI, 90.0))
    FE = NaN;
    EE = 1.0;
elseif isequal(HK, 1.0)
    FE = log((1.0 + sin(D0))/cos(D0));
    EE = sin(D0);
else
    FAC = 1.0;
    for N = 1:40
        A = (A0 + B0) / 2.0;
        B = sqrt(A0*B0);
        C = (A0 - B0) / 2.0;
        FAC = 2.0 * FAC;
        R = R + FAC * C * C;
        if ~isequal(PHI, 90.0)
            D = D0 + atan((B0 / A0)*tan(D0));
            G = G + C * sin(D);
            D0 = D + PI * floor(D/PI+.5);
        end
        A0 = A;
        B0 = B;
        if (C < 1.0e-7)
            break
        end
    end
    CK = PI / (2.0 * A);
    CE = PI * (2.0 - R) / (4.0 * A);

    if isequal(PHI, 90.0)
        FE = CK;
        EE = CE;
    else
        FE = D / (FAC * A);
        EE = FE * CE / CK + G;
    end

end

end
