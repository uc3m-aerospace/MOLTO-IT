function [ r , v , t , a , e , w , at , psi ,DV, flag] = propagate_spirals_DV ( K1, K2, ee, regime0, r0, theta0, t0,  theta )
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
flag = 0;
%
% Initial flight path angle
%
cpsi0   = regime0 * ( sqrt ( ( 2 * ( 1 - ee ) +  K1 * r0 )^2 - K2^2 ) / ( 2 * ( 1 - ee ) + K1 * r0 ) );
spsi0   = K2 / ( 2 * ( 1 - ee ) + K1 * r0 );
%
psi0    =  atan2( real(spsi0) , real(cpsi0) ) ;
%
if flag == 1
    
elseif K1  < -abs(eps) % Elliptical spiral
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE PARTICLE NEVER SCAPE (rmax) when propagated backward or forward
    % the particle reach the origin of the central body
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rmax    = ( 2 * ( 1- ee) - K2 ) / ( -K1 );
    %
    ell     = sqrt( 4 * ( 1 - ee )^2 - K2^2 );
    %
    thetam  = theta0 + regime0 * K2 / ell * abs( acosh ( rmax / r0 - 2 * ( 1 - ee) / K2 * ( 1 - rmax / r0) )) ;
    regime  = -1 * (theta > thetam) + 1 * (theta <= thetam) ;
    beta    = ell / K2 * ( theta - thetam );
    %
    % Trajectory
    %
    r       = rmax * ( 2 * ( 1 - ee) + K2 ) ./ ( 2 * ( 1 - ee ) + K2 * cosh ( beta ) )  ;
    v       = sqrt ( K1 + 2 ./ r * ( 1 - ee ) );
    %
    cpsi    =  regime .* ( sqrt ( ( 2 * ( 1 - ee ) +  K1 * r ).^2 - K2^2 ) ./ ( 2 * ( 1 - ee ) + K1 * r ) );
    spsi    =  K2 ./ ( 2 * ( 1 - ee ) + K1 * r );
    psi     =  atan2( spsi , cpsi ) ;
    
elseif K1 > abs(eps) % Hiperbolic  Spiral
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE PARTICLE REACHES INFINITY WITH A NON-ZERO VELOCITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if  K2 < 2 * ( 1 - ee )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Type I Hiperbolic Spiral ( It has an asymthotic value )
        % Raising regime[ reaches infinity ], Lowering regime[ reaches origin ]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ell     = sqrt( 4 * ( 1 - ee )^2 - K2^2 );
        %
        zeta    = 2 * ( 1 - ee ) + ell ;
        thetaAs = theta0 + regime0 * K2 / ell * log( K2 * ( zeta - ell - K2 * sin( psi0 ) + ell * abs( cos(psi0) ) ) / ...
            ( r0 * K1 * zeta * sin(psi0) ) );
        %
        beta    = regime0 * ell / K2 * ( thetaAs - theta );
        %
        r      = zeta * ell^2 ./  ( K1 * ( sinh( beta / 2) .* ( 4 * zeta * ( 1 - ee ) * sinh( beta / 2 ) + ...
            ( zeta^2 - K2^2 ) * cosh( beta / 2 ) ) ) )  ;
        v      = sqrt ( K1 + 2 ./ r * ( 1 - ee ) );
        %
        cpsi   = regime0 * ( sqrt ( ( 2 * ( 1 - ee ) +  K1 * r ).^2 - K2^2 ) ./ ( 2 * ( 1 - ee ) + K1 * r ) );
        spsi   = K2 ./ ( 2 * ( 1 - ee ) + K1 * r );
        %
        psi    =  atan2(spsi , cpsi ) ;
        
        
    elseif  K2 > 2 * ( 1 - ee )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Type II Hiperbolic Spiral
        % There is a transition between raising and lowering regime [rmin]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rmin = ( K2 - 2 * ( 1 - ee ) ) / ( K1 );
        %
        ell  = sqrt ( K2^2 - 4 * ( 1 - ee )^2 );
        %
        thetam   = theta0 - regime0 * K2 / ell * ( pi / 2 + atan ( ( 2 * ( 1 - ee ) - K2 * sin(psi0) ) / ...
            (ell * abs ( cos(psi0) ) ) ) );
        %
        beta = ell / K2 * ( theta - thetam );
        %
        r   = rmin * ( ( 2 * ( 1 - ee ) + K2 ) ./ ( 2 * ( 1 - ee ) + K2 * cos(beta) ) );
        v   = sqrt ( K1 + 2 ./ r * ( 1 - ee ) );
        %
        regime  = -1 * (theta < thetam) + 1 * (theta >= thetam) ;
        %
        cpsi   = regime .* ( sqrt ( ( 2 * ( 1 - ee ) +  K1 * r ).^2 - K2^2 ) ./ ( 2 * ( 1 - ee ) + K1 * r ) );
        spsi   = K2 ./ ( 2 * ( 1 - ee ) + K1 * r );
        %
        psi    =  atan2( spsi , cpsi ) ;
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Limit case Hiperbolyc spiral
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        thetaAs  = theta0 - regime0 * ( 1 - sqrt( 1 + 4 * ( 1 - ee ) / ( K1 * r0 ) ) );
        beta     = theta - thetaAs ;
        %
        r       = 4 * ( 1 - ee )./ ( K1 * beta .* ( beta - regime0 * 2 ) );
        v       = sqrt ( K1 + 2 ./ r * ( 1 - ee ));
        %
        cpsi   = regime0 .* ( sqrt ( ( 2 * ( 1 - ee ) +  K1 * r ).^2 - K2^2 ) ./ ( 2 * ( 1 - ee ) + K1 * r ) );
        spsi   = K2 ./ ( 2 * ( 1 - ee ) + K1 * r );
        %
        psi    =  atan2( spsi , cpsi ) ;
        
    end
    
elseif abs(K1)< abs(eps) % Parabolic spiral
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THE FLIGHT PATH ANGLE REMAINS CONSTANT
    % THE PARTICLE REACHES INFINITY WITH A ZERO VELOCITY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    ell     = sqrt( 4 * ( 1 - ee )^2 - K2^2 );
    %
    r       = r0 * exp( regime0 * ell * ( theta - theta0 ) / K2 ) ;
    v       = sqrt ( K1 + 2 ./ r * ( 1 - ee ));
    psi     = psi0*ones(size(r));
    %
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the corresponding orbital elements if desired
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
e  = sqrt(K2^2+1 - 2*K2*sin(psi));
a  = (r./(2*ee-K1*r));
w  = atan2( sin(theta) + K2 * cos( theta + psi ), cos(theta) - K2 * sin(theta + psi) ) + pi;
%
% Obtain the corresponding DV
%
if (abs(theta(end)-theta(1)) < 1e-11)
    % Avoid errors with zero arcs duration
    at    = 0;
    DV    = 0;
    %
else
    %
    at         = 1 ./ r.^2 .* sqrt( ee^2 * cos(psi).^2 + ( 1 - 2 * ee )^2 * sin(psi).^2 );
    integrand  = at.*1./(v.*cos(psi)).*r./tan(psi);    
    DV         = sum(integrand(1:end-1).*abs(theta(2:end)-theta(1:end-1)));
    if DV == 0
        DV = NaN;
    end
    %
end

t         = NaN;







