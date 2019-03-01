%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MOLTO-IT TEST CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Name of the Problem to be loaded
%
example = 'Ceres';
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
        planet_arr    = 'Ceres';
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
        useParallel   = 'no';
        options       = [];
        maxGen        = 200;
        popsize       = 50;
        
    case'Jupiter'
        
        problem_name  = example;
        problem_type  = 'flyby';
        planet_dep    = 'Earth';
        planet_arr    = 'Jupiter';
        vinf0_max     =  2;
        planet_fb     = [{'4'},{'3'},{'2'},{'4'},{'3'},{'2'}];
        rfb_min       = 200;
        n_fb          = [0,3];
        rev           = [0,0];
        ToF           = [50,1000];
        ToF_min       = [50  50  50  50 ];
        ToF_max       = [500 500 500 1000];
        Initial_Date  = [{'2029 Jan 01 00:00:00'},{'2030 Dec 31 00:00:00'}];
        init_file     = [];
        output_file   = [problem_name,'.txt'];
        plot          = 0;
        useParallel   = 'yes';
        options       = [];
        maxGen        = 200;
        popsize       = 200;
               
     case'BepiColombo'
        
        problem_name  = example;        
        problem_type  = 'rendezvous';
        planet_dep     = 'Earth';
        planet_arr    = 'Mercury';
        vinf0_max     =  3.36;
        planet_fb     = [{'Mercury'},{'Earth'},{'Venus'}];
        rfb           = [100,7000];
        n_fb          = [4,4];
        ToF           = [100,700];
        Initial_Date  = [{'2018 Jan 01 00:00:00'},{'2018 Dec 31 00:00:00'}];
        init_file     = [];
        output_file   = [problem_name,'.txt'];
        options       = [];
        
      case'BepiColombo_Englander11'
        
        problem_name  = example;        
        problem_type  = 'rendezvous2';
        planet_dep    = 'Earth';
        planet_arr    = 'Mercury';
        vinf0_max     =  1.925;
        vinff_max     =  0.5;
        planet_fb     = [{'Earth'},{'Mercury'},{'Venus'},{'Earth'},{'Mercury'},{'Venus'},{'Earth'},{'Mercury'},{'Venus'}];
        rfb           = [100,7000];
        n_fb          = [2,4];
        ToF           = [100,1900];
        Initial_Date  = [{'2011 Jan 1 00:00:00'},{'2011 Dec 31 00:00:00'}];
        init_file     = ['BepiColombo_Englander10.txt'];
        init_file     = [];
        output_file   = [problem_name,'.txt'];
        plot          = 0;
        useParallel   = 'yes';
        options       = [];
        
    case'Pluto3'
        
        problem_name  = example;
        problem_type  = 'rendezvous';
        planet_dep    = '3';
        planet_arr    = '9';
        vinf0_max     =  8.7681;
        planet_fb     = [{'3'},{'2'},{'7'},{'8'},{'6'},{'5'}];
        %planet_fb     = [{'Jupiter'}];
        planet_fb     = [{'5'}];
        rfb           = [500,7000];
        n_fb          = [1,1];
        ToF           = [100,9000];
        Initial_Date  = [{'2028 Dec 1 00:00:00'},{'2028 Dec 31 00:00:00'}];
        init_file     = ['Pluto.txt'];
        %init_file     = [];
        output_file   = [problem_name,'.txt'];
        plot          = 0;
        useParallel   = 'yes';
        options       = [];
        
end
%
% PREPARE AND SAVE INPUT STRUCTURE
%
save(example);
input = load(example);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output       = molto_it(input);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













