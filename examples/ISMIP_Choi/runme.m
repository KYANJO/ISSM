cwd = pwd
[issmroot,~,~]=fileparts(fileparts(cwd));
newpath=fullfile(issmroot,'/src/m/dev');
addpath(newpath);
devpath;

Lx = 640000;
Ly = 80000;

steps=[1:4];
% steps=[8];
% steps=[7];

nprocs = 48;

ens_id = 0;

loadonly = 0;

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

% Parameterization (Step 2)
if any(steps == 2)
    % filename = fullfile(folder, 'ISMIP.SetMask.mat');
    filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    md = loadmodel(filename);
    md = setflowequation(md, 'SSA', 'all'); % Shelfy-stream approximation
    ParamFile = 'Mismip2.par'
    md = parameterize(md, ParamFile); % Use Mismip2.par
    
    filename = fullfile(folder, 'ISMIP.Parameterization.mat');
    save(filename, 'md');

    % write_netCDF(md, 'ISMIP_Parameterization.nc');
    % export_netCDF(md,'ISMIP_Parameterization.nc');
    % plotmodel(md, 'data', md.geometry.bed, 'title', 'Ice bed_t=0');
end

% Transient Steady state and BC
if any(steps == 3)
    filename = fullfile(folder, 'ISMIP.Parameterization.mat');
    md = loadmodel(filename);

    % % Stressbalance referential
    md.stressbalance.referential = NaN(md.mesh.numberofvertices, 6);
   
    md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices, 1);
    md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices, 1);
  
    md.masstransport.spcthickness = NaN(md.mesh.numberofvertices, 1);
    
    filename = fullfile(folder, 'ISMIP.BC.mat');
    save(filename, 'md');
end


% solve steady state
if any(steps == 4)
    filename = fullfile(folder, 'ISMIP.BC.mat');
    md = loadmodel(filename);

    md=setflowequation(md,'SSA','all');
    
    % Time stepping
    % md.timestepping=timesteppingadaptive();
    % md.timestepping.time_step_max=100;
    % md.timestepping.time_step_min=0.1;
    % md.timestepping.final_time=15000;

    md.timestepping.start_time = 0;
    md.timestepping.time_step = 2;
    md.timestepping.final_time=25000;

    md.settings.output_frequency=100;
    md.stressbalance.maxiter=100;
    md.stressbalance.abstol = NaN;
    md.stressbalance.restol = 1;
    md.settings.solver_residue_threshold=5e-5;
    
    % Verbose
    md.verbose = verbose('all');
    
    md.settings.waitonlock = 0;
    md.cluster=generic('name',oshostname(),'np',nprocs);
    md.cluster.codepath = [issmroot , '/bin'];
    md.cluster.login = 'arobel3';
    md.cluster.executionpath = [issmroot, '/execution'];
    % md.cluster.interactive = 0;
    % md.miscellaneous.name = 'MISMIP';
    % md = solve(md, 'Transient','runtimename',false,'loadonly',loadonly);
    md = solve(md, 'Transient');

    md=loadresultsfromcluster(md);

    filename = fullfile(folder, 'Transient_steadystate1.mat');
    loadonly = 1
    if loadonly == 1
        %md = solve(md, 'Transient','runtimename',false,'loadonly',loadonly);
        save(filename, 'md');
    end
end

% run without adaptive timesteping
if any(steps == 5)

    filename = fullfile(folder, 'Transient_steadystate1.mat');
    md = loadmodel(filename);

    md=setflowequation(md,'SSA','all');
    
    md = transientrestart(md);
    md.geometry.thickness =  md.results.TransientSolution(end).Thickness;
    md.geometry.surface   =  md.results.TransientSolution(end).Surface;
    md.geometry.base      =  md.results.TransientSolution(end).Base;

    % update other fields
    md.initialization.vel      = md.results.TransientSolution(end).Vel;
    md.initialization.vx       = md.results.TransientSolution(end).Vx
    md.initialization.vy       = md.results.TransientSolution(end).Vy
    md.initialization.pressure = md.results.TransientSolution(end).Pressure;
    md.smb.mass_balance        = md.results.TransientSolution(end).SmbMassBalance;
    md.mask.ocean_levelset     = md.results.TransientSolution(end).MaskOceanLevelset;


    md.timestepping.final_time=10000;
    md.settings.output_frequency=100;
    md.stressbalance.maxiter=100;

    md.stressbalance.abstol = NaN;
    md.stressbalance.restol = 1;

    md.verbose = verbose('all');

    md.settings.waitonlock = 0;
    md.cluster=generic('name',oshostname(),'np',nprocs);
    md.cluster.codepath = [issmroot , '/bin'];
    md.cluster.login = 'arobel3';
    md.cluster.executionpath = [issmroot, '/execution'];
    % md.miscellaneous.name = 'MISMIP'
    md.cluster.interactive = 0;
    loadonly = 0;
    md = solve(md, 'Transient','runtimename',false,'loadonly',loadonly);
    md=loadresultsfromcluster(md);

    loadonly = 1
    if loadonly == 1
        md = solve(md, 'Transient','runtimename',false,'loadonly',loadonly);
        filename = fullfile(folder, 'Transient_steadystate2.mat');
        save(filename, 'md');
    end
end


% plotting
if any(steps == 6)
    filename = fullfile(folder, 'ISMIP.inital_simulation.mat');
    % filename = fullfile(folder, 'ISMIP.BC.mat');
    md = loadmodel(filename);
    figure;
    plotmodel(md, 'data', md.geometry.thickness, 'title', 'Ice Thickness'); hold off;
    figure; 
    plotmodel(md, 'data', md.geometry.bed, 'title', 'Bed Topography');
    % figure;
    % plotmodel(md,'data',md.results.StressbalanceSolution(end).Vel)
    % plotmodel(md,'data',md.results.TransientSolution(end).Vel,'title','Velocity');
    for i = 1:11
        plotmodel(md,'data',md.results.TransientSolution(i).Vel)
    end

    % plot surfaces
    x = md.mesh.x;
    y = md.mesh.y;
    z = md.geometry.thickness;
    
    % Create a triangulated surface using mesh connectivity
    trisurf(md.mesh.elements, x, y, z);
    title('Ice Thickness Surface');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('Thickness (m)');
    view(3); shading interp; colorbar;
    
    z = md.geometry.bed;
    trisurf(md.mesh.elements, x, y, z);
    title('Bed Topography');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('Elevation (m)');
    view(3); shading interp; colorbar;
    
    vx = md.results.TransientSolution(end).Vx;
    vy = md.results.TransientSolution(end).Vy;
    velocity_magnitude = sqrt(vx.^2 + vy.^2);
    
    trisurf(md.mesh.elements, x, y, velocity_magnitude);
    title('Velocity Magnitude');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('|Velocity| (m/yr)');
    view(3); shading interp; colorbar;

end

if any(steps == 7)
    filename = fullfile(folder, 'ISMIP.inital_simulation.mat');
    % filename = fullfile(folder, 'ISMIP.BC.mat');
    md = loadmodel(filename);

    % Define interpolation grid
    [X, Y] = meshgrid(1:1000:700000, 1:500:80000);  % Similar to Python np.meshgrid
    
    % Interpolate model fields onto regular grid
    bed_grid = griddata(md.mesh.x, md.mesh.y, md.geometry.bed, X, Y, 'linear');
    base_grid = griddata(md.mesh.x, md.mesh.y, md.geometry.base, X, Y, 'linear');
    surface_grid = griddata(md.mesh.x, md.mesh.y, md.geometry.surface, X, Y, 'linear');
    
    % Convert to kilometers
    X_km = X / 1000;
    Y_km = Y / 1000;
    
    % Create figure
    fig = figure('Color', 'w', 'Position', [100, 100, 2000, 1200]);
    % fig = figure('Position', [100, 100, 1200, 720]);
    ax = axes('Parent', fig);
    hold on;
    
    % Plot surfaces
    surf(ax, X_km, Y_km, bed_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf(ax, X_km, Y_km, base_grid, 'EdgeColor', 'none', 'FaceColor', 'interp');
    surf_handle = surf(ax, X_km, Y_km, surface_grid, ...
        'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.6);
    
    % Set color map and color limits
    colormap(jet);
    caxis([-800, 800]);
    
    % Add colorbar
    cbar = colorbar('Location', 'eastoutside');
    cbar.Label.String = 'Elevation (m)';
    cbar.Label.FontSize = 14;
    
    % Set viewing angle
    view(ax, 30, -60);
    
    % Set aspect ratio (roughly matching Python's ax.set_box_aspect)
    % daspect([7, 0.8, 2]);
    
    % Remove pane fills (simulates transparency)
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.ZColor = 'k';
    set(ax, 'BoxStyle', 'full', 'Box', 'off');
    
    % Ticks
    yticks([0, 40, 80]);
    
    % Axis labels
    xlabel(ax, 'x (km)', 'FontSize', 14);
    ylabel(ax, 'y (km)', 'FontSize', 14);
    zlabel(ax, 'z (m)', 'FontSize', 14);
    
    % Axis label padding and alignment (not identical to Python but close)
    ax.LabelFontSizeMultiplier = 1.2;
    ax.TitleFontSizeMultiplier = 1.2;
    
    % Set font size of all tick labels
    set(ax, 'FontSize', 14);
    
    % Optional: Lighting (MATLAB has lighting tools but not pane edge controls like matplotlib)
    camlight headlight;
    lighting gouraud;
    
    hold off;

end

if any(steps ==8)
    filename = fullfile(folder, 'ISMIP.inital_simulation.mat');
    % filename = fullfile(folder, 'ISMIP.BC.mat');
    md = loadmodel(filename);

    write_netCDF(md, 'ISMIP_initial.nc');

    % Define output filename
    ncfile = fullfile(folder, 'ISMIP.initial_simulation.nc');
    
    % Define basic dimensions
    nNodes = md.mesh.numberofvertices;
    nElems = md.mesh.numberofelements;
    
    % Create and write coordinates
    nccreate(ncfile, 'x', 'Dimensions', {'node', nNodes});
    nccreate(ncfile, 'y', 'Dimensions', {'node', nNodes});
    ncwrite(ncfile, 'x', md.mesh.x);
    ncwrite(ncfile, 'y', md.mesh.y);
    
    % Create and write element connectivity
    nccreate(ncfile, 'elements', 'Dimensions', {'elem', nElems, 'nodes_per_elem', 3});
    ncwrite(ncfile, 'elements', md.mesh.elements);
    
    % Create and write geometry fields
    fields = {'thickness', 'bed', 'surface', 'base'};
    for k = 1:length(fields)
        fieldname = fields{k};
        if isfield(md.geometry, fieldname)
            nccreate(ncfile, fieldname, 'Dimensions', {'node', nNodes});
            ncwrite(ncfile, fieldname, md.geometry.(fieldname));
        end
    end
    
    % Optional: add velocity or other results
    if isfield(md.results, 'TransientSolution')
        vel = md.results.TransientSolution(end).Vel;
        nccreate(ncfile, 'velocity', 'Dimensions', {'node', nNodes, 'dim', 2});
        ncwrite(ncfile, 'velocity', vel);
    end
end

if any(steps==9)
    filename = fullfile(folder, 'ISMIP.inital_simulation.mat');
    % filename = fullfile(folder, 'ISMIP.BC.mat');
    md = loadmodel(filename);

    % Initialize GIF file
    filename = 'velocity_movie.gif';
    for i = 1:64
        plotmodel(md, 'data', md.results.TransientSolution(i).Vel); % Generate the plot
        drawnow;
        frame = getframe(gcf); % Capture frame
        im = frame2im(frame); % Convert to image
        [imind, cm] = rgb2ind(im, 256); % Convert to indexed image
        if i == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
        end
    end
end
