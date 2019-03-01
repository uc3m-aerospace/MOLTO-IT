% This function simply load the required kernels files. "path2spice" is the
% path in your computer to the MICE main folder (e.g. C:/documents/mice)
function ME = load_spice_kernels ( path2mice )


% Add the path to the SPICE (MICE) library and matlab source files
mice_src = [path2mice filesep 'src' filesep 'mice' filesep ];
mice_lib = [path2mice filesep 'lib' filesep ];
kernels  = [path2mice filesep 'kernels' filesep ];
addpath(mice_src);
addpath(mice_lib);
%
% Load the required data:
%
try
    % Leap seconds information and computation of the ephemeris time from UTC
    cspice_furnsh([kernels,'naif0010.tls'])  
    % Rotational elements for reference frames                                                                  
    cspice_furnsh([kernels,'pck00010.tpc'])  
    % Ephemerides of planets, Moon and Sun
    cspice_furnsh([kernels,'de421.bsp'])
    % Ephemerides of planets, Moon and Sun
    cspice_furnsh([kernels,'de430.bsp'])
    % Ephemerides of Jupiter
    cspice_furnsh([kernels,'jup230.bsp'])
    % Ephemerides of planets, Moon and Sun
    cspice_furnsh([kernels,'de432.bsp'])
    % Ephemerides of planets, Moon and Sun
    cspice_furnsh([kernels,'de432s.bsp'])
    % Gravitational constants of the planets
    cspice_furnsh([kernels,'gm_de431.tpc'])
    % Ephemerides of Ceres
    cspice_furnsh([kernels,'2000001.bsp'])
    %Ephemerides of Pluto
    cspice_furnsh([kernels,'2000999.bsp'])
    %Ephemerides of Comet 67P/Churyumov-Gerasimenko
    cspice_furnsh([kernels,'1000012.bsp'])
    %Ephemerides of Asteroid Didimos
    cspice_furnsh([kernels,'2065803.bsp'])
    %Ephemerides of Asteroid Ryugu
    cspice_furnsh([kernels,'2162173.bsp'])
    % If all has been loaded fine, do nothing
    ME = 'No error';
    %
catch ME % ME is an identifier containing the error message

end

return

