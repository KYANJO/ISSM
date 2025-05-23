function runme(nprocs)

	%  read kwargs from a .mat file
	kwargs = load('model_kwargs_0.mat');

	%  access the values of the dictionary
	ParamFile = char(kwargs.ParamFile);
	Lx = double(kwargs.Lx); % length of the domain in x direction
	Ly = double(kwargs.Ly); % length of the domain in y direction
	nx = double(kwargs.nx); % number of nodes in x direction
	ny = double(kwargs.ny); % number of nodes in y direction

    %which steps to perform; steps are from 1 to 8
    %step 7 is specific to ISMIPA
    %step 8 is specific to ISMIPF
    
    % Add ISSM to MATLAB path if not already there
    % issm_dir = getenv('ISSM_DIR');
    % if ~isempty(issm_dir)
    %     addpath(genpath(issm_dir));
    %     disp(['Added ISSM directory and subdirectories from path: ', issm_dir]);
    % else
    %     error('ISSM_DIR is not set.');
    % end
    
    % steps=[1:4];
    % call the ISSM model
    steps=[1:7]; %ISMIPA
    % steps=[];
    % steps=[1:6, 8]; %ISMIPB
    
    % parameter file to be used choose between IsmipA.par or IsmipF.par
    % ParamFile='IsmipA.par'
    %ParamFile='IsmipF.par'
    
    % Run steps
    % Mesh generation #1
    if any(steps==1)
        %initialize md as a new model #help model
	    %->
	    md=model();
	    % generate a squaremesh #help squaremesh
	    % Side is 80 km long with 20 points
	    %->
	    if(ParamFile=='IsmipA.par'),
		    % md=squaremesh(md,80000,80000,20,20);
		    md=squaremesh(md,Lx,Ly,nx,ny);
	    elseif(ParamFile=='IsmipF.par'),
		    % md=squaremesh(md,100000,100000,30,30);
			md = squaremesh(md,Lx,Ly,nx,ny);
	    end
	    % plot the given mesh #plotdoc
	    %->
	    plotmodel(md,'data','mesh')
	    % save the given model
	    %->
	    save ./Models/ISMIP.Mesh_generation md;
    end
    
    % Masks #2
    if any(steps==2)
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    md = loadmodel('./Models/ISMIP.Mesh_generation');
	    % set the mask #help setmask
	    % all MISMIP nodes are grounded
	    %->
	    md=setmask(md,'','');
	    % plot the given mask #md.mask to locate the field
	    %->
	    plotmodel(md,'data',md.mask.ocean_levelset);
	    % save the given model
	    %->
	    save ./Models/ISMIP.SetMask md;
    end
    
    %Parameterization #3
    if any(steps==3)
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    md = loadmodel('./Models/ISMIP.SetMask');
	    % parametrize the model # help parameterize
	    % you will need to fil-up the parameter file defined by the
	    % ParamFile variable
	    %->
	    md=parameterize(md,ParamFile);
	    % save the given model
	    %->
	    save ./Models/ISMIP.Parameterization md;
    end
    
    %Extrusion #4
    if any(steps==4)
	    
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    md = loadmodel('./Models/ISMIP.Parameterization');
	    % vertically extrude the preceding mesh #help extrude
	    % only 5 layers exponent 1
	    %->
	    md=extrude(md,5,1);
	    % plot the 3D geometry #plotdoc
	    %->
	    plotmodel(md,'data',md.geometry.base)
	    % save the given model
	    %->
	    save ./Models/ISMIP.Extrusion md;
    end
    
    %Set the flow computing method #5
    if any(steps==5)
    
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    md = loadmodel('./Models/ISMIP.Extrusion');
	    % set the approximation for the flow computation #help setflowequation
	    % We will be using the Higher Order Model (HO)
	    %->
	    md=setflowequation(md,'HO','all');
	    % save the given model
	    %->
	    save ./Models/ISMIP.SetFlow md;
    end
    
    %Set Boundary Conditions #6
    if any(steps==6)
    
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    md = loadmodel('./Models/ISMIP.SetFlow');
	    % dirichlet boundary condition are known as SPCs
	    % ice frozen to the base, no velocity	#md.stressbalance
	    % SPCs are initialized at NaN one value per vertex
	    %->
	    md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	    %->
	    md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	    %->
	    md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);
	    % extract the nodenumbers at the base #md.mesh.vertexonbase
	    %->
	    basalnodes=find(md.mesh.vertexonbase);
	    % set the sliding to zero on the bed
	    %->
	    md.stressbalance.spcvx(basalnodes)=0.0;
	    %->
	    md.stressbalance.spcvy(basalnodes)=0.0;
	    % periodic boundaries have to be fixed on the sides
	    % Find the indices of the sides of the domain, for x and then for y
	    % for x
	    % create maxX, list of indices where x is equal to max of x (use >> help find)
	    %->
	    maxX=find(md.mesh.x==max(md.mesh.x));
	    % create minX, list of indices where x is equal to min of x
	    %->
	    minX=find(md.mesh.x==min(md.mesh.x));
	    % for y
	    % create maxY, list of indices where y is equal to max of y
	    %  but not where x is equal to max or min of x
	    % (i.e, indices in maxX and minX should be excluded from maxY and minY)
	    %->
	    maxY=find(md.mesh.y==max(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));
	    % create minY, list of indices where y is equal to max of y
	    % but not where x is equal to max or min of x
	    %->
	    minY=find(md.mesh.y==min(md.mesh.y) & md.mesh.x~=max(md.mesh.x) & md.mesh.x~=min(md.mesh.x));
	    % set the node that should be paired together, minX with maxX and minY with maxY
	    % #md.stressbalance.vertex_pairing
	    %->
	    md.stressbalance.vertex_pairing=[minX,maxX;minY,maxY];
	    if (ParamFile=='IsmipF.par')
		    % if we are dealing with IsmipF the solution is in
		    % masstransport
		    md.masstransport.vertex_pairing=md.stressbalance.vertex_pairing;
	    end
	    % save the given model
	    %->
	    save ./Models/ISMIP.BoundaryCondition md;
    end
    
    %Solving #7
    if any(steps==7)
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    md = loadmodel('./Models/ISMIP.BoundaryCondition');
	    % Set cluster #md.cluster
	    % generic parameters #help generic
	    % set only the name and number of process
	    %->
	    md.cluster=generic('name',oshostname(),'np',nprocs);
	    % Set which control message you want to see #help verbose
	    %->
	    md.verbose=verbose('convergence',true);
	    % Solve #help solve
	    % we are solving a StressBalanc
	    %->
	    md=solve(md,'Stressbalance');
	    % save the given model
	    %->
	    save ./Models/ISMIP.StressBalance md;
	    % plot the surface velocities #plotdoc
	    %->
	    plotmodel(md,'data',md.results.StressbalanceSolution.Vel)
    end
    
    %Solving #8
    if any(steps==8)
	    % load the preceding step #help loadmodel
	    % path is given by the organizer with the name of the given step
	    %->
	    md = loadmodel('./Models/ISMIP.BoundaryCondition');
	    % Set cluster #md.cluster
	    % generic parameters #help generic
	    % set only the name and number of process
	    %->
	    md.cluster=generic('name',oshostname(),'np',2);
	    % Set which control message you want to see #help verbose
	    %->
	    md.verbose=verbose('convergence',true);
	    % set the transient model to ignore the thermal model
	    % #md.transient
	    %->
	    md.transient.isthermal=0;
	    % define the timestepping scheme
	    % everything here should be provided in years #md.timestepping
	    % give the length of the time_step (4 years)
	    %->
	    md.timestepping.time_step=4;
	    % give final_time (20*4 years time_steps)
	    %->
	    md.timestepping.final_time=4*20;
	    % Solve #help solve
	    % we are solving a TransientSolution
	    %->
	    md=solve(md,'Transient');
	    % save the given model
	    %->
	    save ./Models/ISMIP.Transient md;
	    % plot the surface velocities #plotdoc
	    %->
	    plotmodel(md,'data',md.results.TransientSolution(20).Vel)
    end
end