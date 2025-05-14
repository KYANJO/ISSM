Lx = 640000;
Ly = 80000;

steps=[1:2];

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
    md = bamg(md, 'domain', 'Domain.exp', 'hmax', 10000, 'hmin', 500, 'splitcorners', 0);
    
    % plotmodel(md,'data','mesh');
    % Save mesh
    filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    save(filename, 'md');
end

% % Masks (Step 2)
% if any(steps == 2)
%     filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
%     md = loadmodel(filename);
%     md = setmask(md, '', ''); % All grounded, no ice shelves
%     filename = fullfile(folder, 'ISMIP.SetMask.mat');
%     save(filename, 'md');
% end

% Parameterization (Step 2)
if any(steps == 2)
    % filename = fullfile(folder, 'ISMIP.SetMask.mat');
    filename = fullfile(folder, 'ISMIP.Mesh_generation.mat');
    md = loadmodel(filename);
    md = setflowequation(md, 'SSA', 'all'); % Shelfy-stream approximation
    ParamFile = 'Mismip2.par'
    md = parameterize(md, ParamFile); % Use Mismip2.par
    % Initialize friction coefficient (Weertman law, reference value)
    % md.friction = frictionweertman();
    % md.friction.coefficient = 2500 * ones(md.mesh.numberofvertices, 1);
    % md.friction.m = 3; % Weertman exponent
    % filename = fullfile(folder, 'ISMIP.Parameterization.mat');
    % save(filename, 'md');
end