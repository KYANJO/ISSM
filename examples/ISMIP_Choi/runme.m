Lx = 640000;
Ly = 80000;

steps=[1:3];
% steps=[6];

ens_id = 0;

folder = sprintf('./Models/ens_id_%d', ens_id);
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Mesh generation (Step 1)
if any(steps == 1)
    md = model();

    % Define the domain outline (x, y coordinates, closing the contour)
    domain = [
        0, 0, 0
        0, Lx, 0;      % Point 1
        Lx, Lx, Ly;    % Point 2
        Lx, 0, Ly;     % Point 3
        0, 0, 0;       % Point 4
    ];
    
    % Define the output filename
    filename = 'Domain.exp';
    
    % Open file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Unable to open %s for writing', filename);
    end
    
    % Write header
    fprintf(fid, '## Name:DomainOutline\n');
    fprintf(fid, '## Icon:0\n');
    fprintf(fid, '# Points Count  Value\n');
    fprintf(fid, '%d 1\n', size(domain, 1));
    fprintf(fid, '# X pos Y pos\n');
    
    % Write points (using columns 2 and 3 for x, y)
    for i = 1:size(domain, 1)
        fprintf(fid, '%f %f\n', domain(i, 2), domain(i, 3));
    end
    
    % Close file
    fclose(fid);

    % Use bamg for variable-resolution mesh (500 m to 10 km)
    % md = bamg(md, 'domain', './Domain.exp', 'hmax', 2000, 'splitcorners', 1);
    hvertices=[10000;500;500;10000];
    md = bamg(md, 'domain', 'Domain.exp', 'hvertices',hvertices);
    % md = bamg(md, 'domain', 'Domain.exp', 'hmax', 10000, 'hmin', 500, 'splitcorners', 0);
    % md = triangle(md,'Domain.exp',27000);
    
    % plotmodel(md,'data','mesh');
    % Save mesh
    filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    save(filename, 'md');
end

% Masks (Step 2)
if any(steps == 2)
    filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    md = loadmodel(filename);
    md = setmask(md, '', ''); % All grounded, no ice shelves
    % plotmodel(md,'data',md.mask.ocean_levelset);
    filename = fullfile(folder, 'ISMIP.SetMask.mat');
    save(filename, 'md');
end

% Parameterization (Step 2)
if any(steps == 3)
    filename = fullfile(folder, 'ISMIP.SetMask.mat');
    % filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    md = loadmodel(filename);
    md = setflowequation(md, 'SSA', 'all'); % Shelfy-stream approximation
    ParamFile = 'Mismip2.par'
    md = parameterize(md, ParamFile); % Use Mismip2.par
    
    filename = fullfile(folder, 'ISMIP.Parameterization.mat');
    save(filename, 'md');
end

% Transient Steady state and BC
if any(steps == 4)
    filename = fullfile(folder, 'ISMIP.Parameterization.mat');
    md = loadmodel(filename);

   % Find front at x ~ Lx
    % ice_front = find(md.mesh.x > (Lx - 1e3));  % within 1 km of x = Lx
    % md.stressbalance.spcvx = NaN(md.mesh.numberofvertices,1);
    % md.stressbalance.spcvy = NaN(md.mesh.numberofvertices,1);
    % md.stressbalance.spcvz = NaN*ones(md.mesh.numberofvertices,1);
    % md.stressbalance.spcvx(ice_front) = 0;
    % md.stressbalance.spcvy(ice_front) = 0;
    %  md.stressbalance.spcvz(ice_front) = 0;
    % 
    % % Stressbalance referential
    md.stressbalance.referential = NaN(md.mesh.numberofvertices, 6);
    % md.stressbalance.loadingforce = zeros(md.mesh.numberofvertices, 3);

    % Forcings
    % md.smb.mass_balance = 0.3 * ones(md.mesh.numberofvertices, 1);
    md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices, 1);
    md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);
    md.groundingline.migration = 'SubelementMigration';
    % md.calving.calvingrate = zeros(md.mesh.numberofvertices, 1);
    % md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices, 1);
    % 
    % % Thickness constraints
    md.masstransport.spcthickness = NaN(md.mesh.numberofvertices, 1);
    % 
    % % Levelset constraints
    % md.levelset.spclevelset = NaN(md.mesh.numberofvertices, 1);


    filename = fullfile(folder, 'ISMIP.BC.mat');
    save(filename, 'md');
end


% solve steady state
if any(steps == 5)
    filename = fullfile(folder, 'ISMIP.BC.mat');
    md = loadmodel(filename);

     % Transient settings
    % md.transient.ismasstransport = 1;
    % md.transient.isstressbalance = 1;
    % md.transient.isthermal = 0;
    % md.transient.isgroundingline = 1;
    % md.transient.isesa = 0;
    % md.transient.ismovingfront = 0;
   
    
    % Time stepping
     md.timestepping.start_time = 0;
    md.timestepping.time_step = 2; % years
    md.timestepping.final_time = 200;
    md.settings.output_frequency = 10;
    
    % Verbose
    md.verbose = verbose('all');

    md.cluster=generic('name',oshostname(),'np',4);

    md = solve(md, 'Transient');
    % md=solve(md,'Stressbalance');
    
    figure;
    plotmodel(md, 'data', md.geometry.thickness, 'title', 'Ice Thickness');
    figure;
    plotmodel(md, 'data', md.geometry.bed, 'title', 'Bed Topography');
    figure;
    % plotmodel(md,'data',md.results.StressbalanceSolution(end).Vel)
    plotmodel(md,'data',md.results.TransientSolution(end).Vel)

    filename = fullfile(folder, 'ISMIP.inital_simulation.mat');
    save(filename, 'md');
end

% plotting
if any(steps == 6)
    % filename = fullfile(folder, 'ISMIP.inital_simulation.mat');
    filename = fullfile(folder, 'ISMIP.BC.mat');
    md = loadmodel(filename);
    figure;
    plotmodel(md, 'data', md.geometry.thickness, 'title', 'Ice Thickness');
    figure;
    plotmodel(md, 'data', md.geometry.bed, 'title', 'Bed Topography');
    figure;
    % plotmodel(md,'data',md.results.StressbalanceSolution(end).Vel)
    plotmodel(md,'data',md.results.TransientSolution(end).Vel,'title','Velocity');

end

