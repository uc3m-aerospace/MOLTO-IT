%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MOLTO-IT TEST CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Name of the Problem to be loaded
%
example   = 'Jupiter';
spice_dir = '/Users/davidmorante/Desktop/MicePackage/';
%
%planets = [ {'1'}   % mercury
%            {'2'}   % venus
%            {'3'}   % earth
%            {'4'}   % mars
%            {'5'}   % jupiter
%            {'6'}   % saturno
%            {'7'}   % urano
%            {'8'}   % neptuno
%            {'9'}   % pluton
          
switch(example)
    %       
    case'Ceres'
        
        problem_name  = example;
        problem_type  = 'rendezvous';
        planet_dep    = 'Earth';
        planet_arr    = '2000001'; %Ceres
        vinf0_max     =  1.6;
        planet_fb     = {'4'};
        rfb_min       = 200;
        n_fb          = [0,1];
        rev           = [0,0];
        ToF           = [100,1000];
        Initial_Date  = [{'2003 Jan 01 00:00:00'},{'2003 Dec 31 00:00:00'}];
        init_file     = [];
        output_file   = [problem_name,'.txt'];
        plot          = 0;
        useParallel   = 'yes';
        options       = [];
        maxGen        = 200;
        popsize       = 50;
        
    case'Jupiter'
        
        problem_name  = example;
        problem_type  = 'flyby';
        planet_dep    = 'Earth';
        planet_arr    = '5';
        vinf0_max     =  2;
        planet_fb     = [{'4'},{'3'},{'2'},{'4'},{'3'},{'2'}];
        rfb_min       = 200;
        n_fb          = [0,3];
        rev           = [0,0];
        ToF           = [50  50  50  50; 
                          500 500 500 1000];
        Initial_Date  = [{'2029 Jan 01 00:00:00'},{'2030 Dec 31 00:00:00'}];
        init_file     = [];
        output_file   = [problem_name,'.txt'];
        plot          = 0;
        useParallel   = 'yes';
        options       = [];
        maxGen        = 200;
        popsize       = 200;
              
end
%
% PREPARE AND SAVE INPUT STRUCTURE
%
input = load(example);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output       = molto_it(input);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













