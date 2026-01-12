%     mass_spring_2d_data
%
% the following variables are declared as global
% to be used by all the functions in this program
%
global nmx   % number of masses in the horizontal direction
global nmy   % number of masses in the vertical direction
global nm    % total number of masses
global ns    % total number of springs
global nms   % number of masses per spring
global mdof  % number of DOFs per mass
global spdof % number of DOFs per spring
global nfdof  % number of free DOFs
global dep_dofs % number of DOFs contrained by periodicity conditions
global geom  % node coordinates X and Y
global connec% connectivity matrix
global nf    % nodal freedom matrix
global load  % matrix of external loading
global load_ampl_x % amplitude of loading, X direction
global load_ampl_y % amplitude of loading, Y direction
global T     % global transformation matrix to apply periodic BCs

%% Data Initialization %%%%%%%%%%%%%%%%%%%%
%
%  the total number of masses
nm = nmx*nmy;
% the total number of springs is not known at this stage, and will be estimated
% later, while composing the connectivity matrix

% the number of masses per spring
nms = 2;
% the number of DOF per mass
mdof = 2;
% the number of DOFs per spring
spdof = nms*mdof;
% vertical distance between horizontal raws of masses
h = sqrt(3)/2*a; 
%
% Node coordinates X and Y 
%
% create an array for storing the mass corrdinates
geom = zeros(nm,2);
%
% define random shifts with max 25% of displacement from periodic position
% rand - generates univormly distributed random numbers in the interval between 0 and 1
rgeom = 0.35*rand(nmx,nmy);
%
% boundary faces must be undisturbed to enable the application of boundary
% conditions
rgeom(1,:) = 0;
rgeom(nmx,:) = 0;
rgeom(:,1) = 0;
rgeom(:,nmy) = 0;
%
% perturb the lattice nodes
mn = 0; % mass number
for i = 0:nmy-1
    for j = 0:nmx-1
        geom(mn+1,1)= (j+0.5*mod(i,2))*a + ((rand > 0.5)*2 - 1)*rgeom(j+1,i+1);
        geom(mn+1,2) = i*h + ((rand > 0.5)*2 - 1)*rgeom(j+1,i+1);
        mn = mn+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
%
% Element connectivity
%
connect = []; % empty array to store the connectivity elements, as the number of springs is unknown yet
%
% calculate the number of springs
ns = 0;
% specify horizontal strings
for i = 0:nmy-1
    for j = 1:nmx-1
        connect = [connect; nmx*i+j, nmx*i+j+1];
        ns = ns + 1;
    end
end
% specify 'right-pointed' strings
for i = 0:nmy-2
    for j = 1:nmx
        connect = [connect; nmx*i+j, nmx*(i+1)+j];
        ns = ns + 1;
    end
end
% specify 'left-pointed' strings - part 1
for i = 0:2:nmy-2
    for j = 1:nmx-1
        connect = [connect; nmx*i+j+1, nmx*(i+1)+j];
        ns = ns + 1;
    end
end
% specify 'left-pointed' strings - part 2
for i = 1:2:nmy-2
    for j = 1:nmx-1
        connect = [connect; nmx*i+j, nmx*(i+1)+j+1];
        ns = ns + 1;
    end
end
%
%% display ns
%ns
%
%% Element connectivity
%
connec = zeros(ns,2);
%
% rewrite the connectivity data to connectivity matrix (this step is redundant)
connec(:,1)=connect(:,1); 
connec(:,2)=connect(:,2); 
%
% Boundary conditions
%
% 1. Preparation
% (a) estimate max/min X and Y
Xmin = min(geom(:,1));
Xmax = max(geom(:,1));
Ymin = min(geom(:,2));
Ymax = max(geom(:,2));
eps = 1e-6;
%
% (b) determine nodes at the right and top boundares
NodesR = find((geom(:,1)>Xmax-eps & geom(:,1)<Xmax+eps));
NodesT = find((geom(:,2)>Ymax-eps & geom(:,2)<Ymax+eps));
%
% (c) create pairs of nodes at the top & bottom boundaries
npairTB = [];
for i=1:length(NodesT)
    x = geom(NodesT(i),1);
    npairTB = [npairTB; NodesT(i), find((geom(:,1)>x-eps & geom(:,1)<x+eps & geom(:,2)>Ymin-eps & geom(:,2)<Ymin+eps))];
end
%
% (d) initialize the nodal freedom matrix: 1 is for a free DOF
nf = ones(nm, mdof);
%
% 2. Prescribe clambed BC at the right boundary
%
for i=1:length(NodesR)
    nf(NodesR(i),:)=0;
end
%
% 3. Count free DOFs: nfdof
% (the free DOFs are counted (different from 0 entries of nf-matrix) and their rank assigned back into
% the matrix nf(nn, nndof)):
%
nfdof = 0;
for i = 1:nm
    for j = 1:mdof
        if nf(i,j) ~= 0
            nfdof = nfdof+1;
            nf(i,j) = nfdof;
        end
    end
end
%
% 3. Apply periodic BC at the top & bottom boundaries
%
%    These are incorporated into transformation matrix defined for the
%    global stiffness matrix
%
% (a) form the transformation matrix T
%
T = sparse(eye(nfdof,nfdof));
%
% (b) the nodes at the top boundary are redunant with dependent dofs
% (dep_dofs)
%
dep_dofs = unique([nf([NodesT],:)]);
%
% (c) for the rows with numbers in dep_dofs insert I matrix 
%
for i=1:length(npairTB)
    for j=1:mdof
        T(nf(npairTB(i,1),j),nf(npairTB(i,2),j)) = 1.;  
    end
end
%
% (d) remove columns in T matrix corresponding to Top nodes
%
T(:,dep_dofs) = [];
%
%
% Apply the loading at the middle node on the left-hand side of the
% network
% 
load = zeros(nm,mdof);
% determine nodes at the left boundary
NodesL = find((geom(:,1)>(Xmin)-eps & geom(:,1)<(Xmin)+eps));
% horizontal excitation at the middle of the left boundary
load(NodesL(fix(length(NodesL)/2)+1),1) = load_ampl_x;
load(NodesL(fix(length(NodesL)/2)+1),2) = load_ampl_y;
%

