function output = genetic_problem(input,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    FUNCTION THAT DEFINES THAT CALLS THE GENETIC SOLVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load paths
%
addpath('/Users/davidmorante/Desktop/SOLVERS/NGPM')
addpath('/Users/davidmorante/Desktop/moiseevigor-elliptic-348471b')
addpath('/Users/davidmorante/Desktop/MicePackage')
%
% Define indices structure
%
n_fb_min  = input.n_fb(1);
n_fb_max  = input.n_fb(2);
planet_fb = input.planet_fb;
%
get_indices
code_planets
%
% Set Minimum and Maximum Values for Variables
%
ToF_min   = input.ToF(1)*ones( 1 , n_fb_max + 1 );
ToF_max   = input.ToF(2)*ones( 1 , n_fb_max + 1 );
%
%ToF_min  =  [ 800, 200];
%ToF_max  =  [ 1600, 1000];
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
%fbb_min = [2 1 2];
%fbb_max =fbb_min;
%
LB = [ fbb_min,  0 , ToF_min];
UB = [ fbb_max,  1 , ToF_max];
%
% Set the vartype (=1 Continuos) (=2 Discrete)
%
vartype = [2*ones(size(fbb_max)),1 ,ones(size(ToF_min))];
%
% Set the Genetic algorithm parameters
%
nvars = numel(LB);
%
% Set the auxiliary parameter in a setup structure
%
input.planet_available = length(planet_fb);
input.vinf0      = input.vinf0_max;
if isfield(input,'vinff_max')
input.vinff      = input.vinff_max;
end
input.ind        = ind;
input.planets    = planets;
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
setup = input;
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Set the NSGA-II genetic algorithm parameters
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
options = nsgaopt();                     % Initialize options structure
options.popsize = 100;                  % population size
options.maxGen  = 500;                   % max generation
options.outputfile =input.output_file;   %outputfile
options.nameObj    = {'Time ( years )','m_p'}; 
options.crossover={'intermediate', 0.5}; % crossover operator (Ratio=1.2)--Intermediate crossover [3] creates two children from two parents: parent1 and parent2.If it lies in the range [0, 1], the children created are within the two parent. If algorithm is premature, try to set ratio larger than 1.0.
options.mutation={'gaussian',0.1, 0.2};  % mutation operator (scale=0.1(deviation of the random number), shrink=0.5)--> for example, shrink?[0.5, 1.0]) is usually used for local search. A large mutation range (shrink == 0) is require getting out of the local Pareto-optimal fronts
options.crossoverFraction = 0.5;         % crossover fraction of variables of an individual ( 2/numVar )-->only crossoverFraction of all variables would do crossover
options.mutationFraction  = 0.2;         % only mutaionFraction of all variables would do mutation (default 2/numvar)
options.numObj = 2;                      % number of objectives
options.numVar = nvars;                  % number of design variables
options.numCons = 1;                     % number of constraints
options.lb = LB;                         % lower bound of x
options.ub = UB;                         % upper bound of x
options.objfun = @fitness_ga;            % objective function handle
options.plotInterval = 1;                % interval between two calls of "plotnsga". 
options.outputInterval = 1; 
options.useParallel =input.useParallel;  % Parallel option. 
options.vartype = vartype;
%
if ~isempty(input.init_file)
    options.initfun =  {@initpop, input.init_file};
end
%
%--------------------------------------------------------------------------
% Load Kernels
%--------------------------------------------------------------------------
%
if strcmp(options.useParallel,'no')
    %
    load_spice_kernels('/Users/davidmorante/Desktop/MicePackage/');
    %
else
    parfor jj  = 1:4
        %
        load_spice_kernels('/Users/davidmorante/Desktop/MicePackage/');
        %
    end
end
%

if input.plot == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = nsga2(options,setup); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
%
result = loadpopfile(input.init_file);
%
% Obtain the size of the population and the last generation
%
[last_gen,pop_size] = size(result.pops);
%
pop = result.pops(last_gen,:);
%
% Sort the population size in increasing order for Time
%
%--------------------------------------------------------------------------
% Obj2 = zeros(pop_size,1);
% 
% for i = 1:pop_size
%     
%    Obj2(i) = pop(i).obj(1);
%     
% end
% 
% [~,I] = sort(Obj2);
% 
% pop   = pop(I);
%--------------------------------------------------------------------------
path  = ['Results_',input.init_file(1:end-4)];
%
if not ( isdir(path) )
   mkdir(path)
end
%
setup.path = path;
%
fid = fopen( [path,'/GA_SUMMARY.txt'], 'wt' );
%
fprintf( fid, '---------------------------------------------------------------------');
fprintf( fid, '\n');
fprintf( fid, '---------------------------------------------------------------------');
fprintf( fid, '\n');
fprintf( fid, '                    GENETIC ALGORITHM SUMMARY                        ');
fprintf( fid, '\n');
fprintf( fid, '---------------------------------------------------------------------');
fprintf( fid, '\n');
fprintf( fid, '---------------------------------------------------------------------');
fprintf( fid, '\n');

for ii = 1:pop_size
    %
    subpath = [path,'/','TOF_',num2str(pop(ii).obj(1))];
    %
    if not ( isdir(subpath) )
        mkdir(subpath)
    end
    %
    setup.subpath = subpath;
    %
    format long
    %
    fig = figure(ii);
   % pop(ii).var = [  1.000000000000	 0.697155793443	508.410698544351	592.29550312123];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Obj,Cons,outdata] = fitness_plot(pop(ii).var,setup);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    axis equal
    grid on
    title(['Trajetory'],'Interpreter','latex','Fontsize',14)
    xlabel('X (AU)','Interpreter','latex','Fontsize',14)
    ylabel('Y (AU)','Interpreter','latex','Fontsize',14)
    saveas(fig,[subpath,'/Trajectory'],'epsc')
    saveas(fig,[subpath,'/Trajectory'],'fig')
    close(fig)
    %
    fig2 = figure(ii);
    plot(outdata.time,outdata.at)
    grid on
    title(['Trajetory'],'Interpreter','latex','Fontsize',14)
    xlabel('Time (years)','Interpreter','latex','Fontsize',14)
    ylabel('at(m/s)','Interpreter','latex','Fontsize',14)
    saveas(fig2,[subpath,'/Thrust_profile'],'epsc')
    saveas(fig2,[subpath,'/Thrust_profile'],'fig')
    close(fig2)
    %
%     fprintf( fid, '%7.5f ',Obj(1)) ;
%     %
%     fprintf( fid, '%7.5f ',Obj(2)) ;
%     %
%     fprintf( fid, outdata.date1(1,:)) ;
%     fprintf( fid, '  ') ;
%     %
%     fprintf( fid, '%7.5f ',outdata.ToF(1)/365) ;
%     %
%     fprintf( fid, outdata.date2(1,:), '  ') ;
%     fprintf( fid, '  ') ;
%     %
%     fprintf( fid, '%7.5f ',outdata.ToF(2)/365) ;
%     %
%     fprintf( fid ,outdata.date2(2,:)) ;
%     %
%     fprintf( fid, '  ') ;
%     fprintf( fid, '%7.5f ',outdata.vrel) ;
%     fprintf( fid, '\n') ;
    %
    jj = 1;
    for jj = 1: n_fb_max+1 & numel(outdata.all)>= jj
        %
        dlmwrite([subpath,'/Trajectory_Leg',num2str(jj),'.output'],outdata.all(jj).data,'delimiter','\t','precision','%.10f')
        %fid2 = fopen( [subpath,'/Trajectory_Leg',num2str(jj),'.output'], 'wt' );
        %fprintf( fid2, '%7.5f ',outdata.all(jj).data);
        %fclose(fid2);
        %
    end
    %
end

fclose (fid);

end






