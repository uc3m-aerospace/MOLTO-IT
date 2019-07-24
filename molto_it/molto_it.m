%--------------------------------------------------------------------------
function output = molto_it(input)
%--------------------------------------------------------------------------
%	MOLTO-IT Software Computation Core
%
%	This program is developed at the Universidad Carlos III de Madrid,
%   as part of a PhD program.
%
%   The software and its components are developed by David Morante González
%
%   The program is released under the MIT License
%
%   Copyright (c) 2019 David Morante González
%
%--------------------------------------------------------------------------
%
%    Function that defined the parameters for NSGA-II and call the
%    required routines to solve the interplanetary transfer
%
%--------------------------------------------------------------------------
%
% Include paths
%
output = [];
addpath(genpath('~/MOLTO-IT'))
input.spice_dir = '~/MOLTO-IT/spice';
%
% Load maximum/minimum number of flybys
%
n_fb_min   = input.n_fb(1);
n_fb_max   = input.n_fb(2);
%
% Load flyby planet list
%
planet_fb  = input.planet_fb;
%
% Determine the position of the flyby planets in the gen
%
ind.fbb(1) = 1;
ind.fbb(2) = ind.fbb(1)-1 + n_fb_max;
%
% Determine the position of the departure date in the gen
%
ind.t0     = ind.fbb(2) + 1;
%
% Determine the position of the flyby times in the gen
%
ind.ToF(1) = ind.t0 + 1;
ind.ToF(2) = ind.ToF(1) + n_fb_max;
%
% Determine the position of the revolution number in the gen
%
ind.rev(1) = ind.ToF(2) + 1;
ind.rev(2) = ind.rev(1) + n_fb_max;
%
%code_planets
%
input.ind        = ind;
input.planets    = planet_fb;
%
% Set Min/Max Values for the flight times
%
if numel(input.ToF) ==2
    ToF_min   = input.ToF(1)*ones( 1 , n_fb_max + 1 );
    ToF_max   = input.ToF(2)*ones( 1 , n_fb_max + 1 );
else
    ToF_min   = input.ToF(1,:);
    ToF_max   = input.ToF(2,:);
end
%
% Check consistency in input results
%
if ~isequal(numel(ToF_min),n_fb_max + 1) || ~isequal(numel(ToF_max),n_fb_max + 1)
    error('Check consistency in input vector')
end
%
% Set Min/Max Values for the flyby bodies
%
fbb_min   = ones(1 , n_fb_max);
fbb_max   = ((length(planet_fb)+1))*ones(1 , n_fb_max);
%
jj = 1;
%
while n_fb_min > 0
    %
    fbb_max(jj) = fbb_max(jj) -1;
    jj = jj+1;
    n_fb_min    = n_fb_min-1;
    %
end
%
% Set Min/Max Values for the initial date
%
t0_min = 0;
t0_max = 1;
%
% Set Min/Max Values for the number of revolutions
%
n_min = input.rev(1)*ones( 1 , n_fb_max + 1 );
n_max = input.rev(2)*ones( 1 , n_fb_max + 1 );
%
%
% Set Min/Max Values for the population
%
LB = [ fbb_min,  t0_min , ToF_min, n_min];
UB = [ fbb_max,  t0_max , ToF_max, n_max];
%
% Set the vartype (=1 Continuos) (=2 Discrete)
%
vartype = [2*ones(size(fbb_min)),1 ,ones(size(ToF_min)),2*ones(size(n_min))];
%
% Set the number of variables
%
nvars = numel(LB);
%
% Set the auxiliary parameter in a setup structure
%
input.planet_avai = length(planet_fb);
input.vinf0       = input.vinf0_max;
%
% Load maximum hyperbolic velocity at destination
%
if isfield(input,'vinff_max')
    input.vinff = input.vinff_max;
end
%
% Load transfer case
%
if strcmp(input.problem_type,'rendezvous')
    input.type = 1;
elseif strcmp(input.problem_type,'flyby')
    input.type = 0;
elseif strcmp(input.problem_type,'rendezvous2')
    input.type = 2;
end
%
input.available_planets = length(planet_fb);
%
% Consider Fixed Launch Date
%
if numel(input.Initial_Date) > 1
    
    input.t0_min  = input.Initial_Date(1);
    input.t0_max  = input.Initial_Date(2);
    
else
    
    input.t0_min  = input.Initial_Date(1);
    input.t0_max  = input.Initial_Date(1);
    
end
%
% Create output folder
%
if ~exist(input.output_dir, 'dir')
    mkdir(input.output_dir)
end

%
%--------------------------------------------------------------------------
% Set the NSGA-II genetic algorithm parameters
%--------------------------------------------------------------------------
%
%  Initialize options structure
%
options = nsgaopt();
%
% Set User defined population size
%
if isfield(input,'popsize')
    options.popsize    = input.popsize;
end
%
% Set User defined Maximum number of generations
%
if isfield(input,'maxGen')
    options.maxGen     = input.maxGen;
end
%
options.numObj  = 3;      % number of objectives
options.numVar  = nvars;  % number of design variables
options.numCons = 1;      % number of constraints
options.lb      = LB;     % lower bound of x
options.ub      = UB;     % upper bound of x
%
options.nameObj    = {'Time ( years )','m_p/m_0','flybys'};  % List of objectives names for the plot
options.crossover  = {'intermediate', 0.5};     % crossover operator (Ratio=1.2)--Intermediate crossover [3] creates two children from two parents: parent1 and parent2.If it lies in the range [0, 1], the children created are within the two parent. If algorithm is premature, try to set ratio larger than 1.0.
options.mutation   = {'gaussian',0.1, 0.2};     % mutation operator (scale=0.1(deviation of the random number), shrink=0.5)--> for example, shrink?[0.5, 1.0]) is usually used for local search. A large mutation range (shrink == 0) is require getting out of the local Pareto-optimal fronts
options.crossoverFraction = 0.8;                % crossover fraction of variables of an individual ( 2/numVar )-->only crossoverFraction of all variables would do crossover
options.mutationFraction  = 0.3;                % only mutaionFraction of all variables would do mutation (default 2/numvar)
options.objfun            = @fitness_nsga2;      % objective function handle
options.plotInterval      = 1;                  % interval between two calls of "plotnsga".
options.outputInterval    = 1;
options.outputfile_dir    = input.output_dir;
options.outputfile        = [input.output_dir,'/Results_extended.txt'];  %outputfile
options.useParallel       = input.useParallel;  % Parallel option.
options.vartype           = vartype;
%
if ~isempty(input.init_file)
    options.initfun =  {@initpop, input.init_file};
end
%
%--------------------------------------------------------------------------
% Load Spiece Kernels for parallel workers
%--------------------------------------------------------------------------
%
if strcmp(options.useParallel,'no')
    %
    load_spice_kernels(input.spice_dir);
else
    %
    %parpool
    delete(gcp('nocreate'));
    p= gcp;
    parfor jj  = 1:p.NumWorkers
        %
        load_spice_kernels(input.spice_dir);
        %
    end
    %
end
%
%--------------------------------------------------------------------------
% Call NSGA-II algorithm
%--------------------------------------------------------------------------
%
if input.plot == 0
    %
    % Run Genetic Algorithm
    %
    output = nsga2(options,input);
else
    %
    % Plot Result of input.plot generation
    %
    % STEP 1: obatin variables from standard file Results_extended.txt
    %
    sol = loadpopfile(options.outputfile);
    %
    % STEP 2: Take last Generation 
    %
    last_gen = sol.pops(end,:);
    %
    % STEP 3: Take desired population defined in the varibale input.plot
    %
    population = last_gen(input.plot);
    x = population.var;
    %
    % STEP 4: Propagate spiral
    %
    [~,~,outdata] = fitness_nsga2(x,input);
    %
    % STEP 5: Plot trajectory
    %
    hh = figure(1);
    hold on
    axis equal
    grid on
    xlabel('X(AU)')
    ylabel('Y(AU)')
    for ii = 1:outdata.n_fb_real +1
        
        data  = outdata.out{ii};
        r     = data(:,2);
        angle = data(:,7);
        at    = data(:,5);
        %data_all = [data_all;t2(2:end)',r(2:end)',v(2:end)',psi(2:end)',at2(2:end)',alpha(2:end)',angle(2:end)'];
        [x,y] = pol2cart(angle,r);
        %
        mark1 = plot(x(1),y(1),'black o');
        set(mark1,'MarkerFaceColor','black')
        set(mark1,'MarkerSize',7)
        %
        mark2 = plot(x(end),y(end),'black o');
        set(mark2,'MarkerFaceColor','black')
        set(mark2,'MarkerSize',7)
        %
        % Plot thrust-coast-thrust sequence
        %
        cond1 = at == 0;
        cond2 = find(cond1(1:end-1) ~= cond1(2:end))';
        idx = [1,cond2,numel(at)];
        for jj = 1:numel(cond2)+1
            if jj ==2
                plot(x(idx(jj):idx(jj+1)),y(idx(jj):idx(jj+1)),'black --')
                
                mark3 = plot(x(idx(jj)),y(idx(jj)),'black o');
                set(mark3,'MarkerFaceColor','white')
                set(mark3,'MarkerSize',7)
                
                mark4 = plot(x(idx(jj+1)),y(idx(jj+1)),'black o');
                set(mark4,'MarkerFaceColor','white')
                set(mark4,'MarkerSize',7)
            else
                plot(x(idx(jj):idx(jj+1)),y(idx(jj):idx(jj+1)),'black -')
            end
        end
    end
    
    saveas(hh,[input.output_dir,'/Trajectory',num2str(input.plot)],'epsc')
    close (hh);   
    %
    % STEP 6: Plot Thrust
    %
    hh2 = figure(2);
    hold on
    grid on
    xlabel('Time(years)')
    ylabel('Acceleration')
    t0   = 0;
    for ii = 1:outdata.n_fb_real +1
        
        data  = outdata.out{ii};
        time  = data(:,1)/(3600*24*365) + t0;
        at    = data(:,5);
        %data_all = [data_all;t2(2:end)',r(2:end)',v(2:end)',psi(2:end)',at2(2:end)',alpha(2:end)',angle(2:end)'];
        plot(time,at)
        t0    = time(end);
    end
    
        saveas(hh2,[input.output_dir,'/Accel',num2str(input.plot)],'epsc')
        close (hh2);
    
    
    
    
end
%
%