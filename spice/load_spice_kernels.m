% This function simply load the required kernels files. "path2spice" is the
% path in your computer to the MICE main folder (e.g. C:/documents/mice)
function ME = load_spice_kernels ( path2mice )


% Add the path to the SPICE (MICE) library and matlab source files
mice_src = [path2mice filesep 'mice/src' filesep 'mice' filesep ];
mice_lib = [path2mice filesep 'mice/lib' filesep ];
s  = what([path2mice filesep 'kernels' filesep ]);
kernels = [s.path filesep];
addpath(mice_src);
addpath(mice_lib);
addpath(kernels);
%
% Load the required data:
%
try
    % Leap seconds information and computation of the ephemeris time from UTC
    cspice_furnsh( [kernels,'naif0010.tls'] )
    % Rotational elements for reference frames
    % Ephmerides for Ceres
    cspice_furnsh( [kernels,'2000001.bsp'] )
    % Ephemerides of planets, Moon and Sun
    cspice_furnsh( [kernels,'de432s.bsp'] )
    % Gravitational constants of the planets
    cspice_furnsh( [kernels,'gm_de431.tpc'] )
    % Rotational elements for reference frames
    cspice_furnsh( [kernels,'pck00010.tpc'] )
    %
    ME = 'No error';
    %
catch ME % ME is an identifier containing the error message

end

return

