% by A.O. Krushynska (akrushynska@gmail.com)
% last access: 12 January 2026
%
%% DYNAMIC LINEAR ANALYSIS OF A 2D RANDOM MASS-SPRING NETWORK
%
% the script generates randomly perturbed (up to 35% displacement) triangular
% networks with periodic boundary conditions at the top and bottom horizontal bounds,
% which is excited by a harmonic load of speficied frequencies in the center of the left-hand side of the lattice,
% and evaluates their dynamic response as a whole and at target nodes
%
%
%% 1) Specify a network geometry
% Define the number of nodes in horizontal and vertical directions 
% by specifying the values for 'nmx' and 'nmy' variables
%% The number in vertical direction must be uneven for a triangular lattice, to allow for periodic BCs

%
%% 2) Define the network mechanics and loading
% 2.1) specify the spring stiffness 'prop', mass of a mass 'mass'
% 2.2) specify the magnitude components 'load_ampl_x', 'load_ampl_y' and
% the excitation frequencies 'load_freqs' in Hz.
%% 
%
clc;         % Clear screen
clear all;   % Clear all variables in memory
close all;
%
% make these variables global to share them with other functions
%
global nmx   % number of masses in the horizontal direction
global nmy   % number of masses in the vertical direction
global nm    % total number of masses
global ns    % total number of springs
global nms   % number of masses per spring
global mdof  % number of DOFs per mass
global spdof % number of DOFs per spring
global nfdof % number of free DOFs
global dep_dofs % number of DOFs contrained by periodicity conditions
global geom  % node coordinates X and Y
global connec% connectivity matrix
global prop  % spring stiffness constant
global mass  % mass of the masses
global nf    % nodal freedom matrix
global load  % matrix of external loading
global load_ampl_x % amplitude of loading, X direction
global load_ampl_y % amplitude of loading, Y direction
global load_freqs  % loading frequencies
global T           % global transformation matrix to apply periodic BCs
global a
%
disp ('Executing Main.m');
%
%% INPUT DATA
% the number of masses along the horizontal and vertical direction
nmx = 21;
nmy = 29;
%
% lattice size
a = 1; % [m]
%
% spring stiffness(per unit length), mass, input amplitudes
prop = 1e6; % [N/m]
mass = 1; % [kg]
load_ampl_x = 1; % [m]
load_ampl_y = 0; % [m]
% the number of excitation frequencies
n_freq = 2;
% input frequencies
load_freqs = zeros(n_freq,1);
load_freqs(1) = 100; % [Hz] 
load_freqs(2) = 200; % [Hz] 
%
% the number of outputs, m
n_outputs = 4;
outputs = zeros(n_outputs,1); 
%
% target outputs at each frequency
t_outputs = zeros(n_freq,n_outputs);
m_2 = 1 + n_outputs/2;
t_outputs(1,1:m_2) = 0.1;
t_outputs(1,m_2+1:n_outputs) = 0.9;
t_outputs(2,1:m_2) = 0.9;
t_outputs(2,m_2+1:n_outputs) = 0.1;
%
%
% open the output file 
fID_stat = fopen('test.txt','w');
% store the network parameters
str_e = sprintf('$ network size: %u %s %u %s \r\n',  nmx, 'x', nmy, ' masses');
fprintf(fID_stat,str_e);
str_e = sprintf('$ target normalized amplitudes at outputs at f1: %1.2f  ',  t_outputs(1,:)); 
fprintf(fID_stat,str_e);
str_e = sprintf('\r\n');
fprintf(fID_stat,str_e);
str_e = sprintf('$ target normalized amplitudes at outputs at f2: %1.2f  ',  t_outputs(2,:)); 
fprintf(fID_stat,str_e);
str_e = sprintf('\r\n');
fprintf(fID_stat,str_e);
str_e = sprintf('$ excitation frequencies (in Hz): %u \r\n',  load_freqs(:));
fprintf(fID_stat,str_e);
fprintf(fID_stat,'%s \r\n',"$ number of prunned strings and the value of cost function S ");
%
% specify the number of networks 
num_networks = 220;
%
for itnum = 1:num_networks
    %
    %display the iteration number
    disp(['network number: ', num2str(itnum),' of ', num2str(num_networks)]);
    %
    % form an input file
    %% to minimize simulation time, do not save the network data in the file
    %% we generate and use them directly
    mass_spring_2d_data;
    %
    %% RANDOMLY choose the output nodes
    % array with the boundary nodes & 0, to exclude them from consideration
    bounds = [0, NodesL(:)', NodesR(:)', npairTB(:,1)', npairTB(:,2)'];
    % generate random output nodes
    outputs = int32(rand(n_outputs,1)*nm);
    % check if the selected nodes are at the boundary
    flag = 1;
    while (flag)
        is_bound = [];
        for i=1:n_outputs
            is_bound_t = find(outputs(i)==bounds);
            is_bound = [is_bound; is_bound_t(:)];
        end
        flag = max(is_bound);
        if(flag)
            outputs = int32(rand(n_outputs,1)*nm);
        end
    end
    %
    str_e = sprintf('$ outputs at f1: %u \r\n',  outputs(:));
    fprintf(fID_stat,str_e);
    str_e = sprintf('$ outputs at f2: %u \r\n',  outputs(:));
    fprintf(fID_stat,str_e);
    %%%%%%%%% end of specification of outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Initialize global stiffness, mass and force matrices for active DOFs
    matrK = sparse(nfdof,nfdof); % the global stiffness to zero
    matrM = sparse(mass*eye(nfdof,nfdof)); % the global mass to zero
    matrF = sparse(nfdof,1);
    %
    % Assemble the global stiffness matrix
    for i=1:ns
        % (a) form element stiffness & mass matrices in local coordinates
        ke = spring_ke(i);
        % (b) form tranform matrix C
        C = spring_C(i);
        % (c) transform the stiffnes element matrix to global coordinates
        kg = C*ke*C';
        % (d) retrieve the element steering vector
        g = spring_g(i);
        % (e) assemble global stiffnes matrix
        matrK = form_matrK(matrK, kg, g);
    end
    %
    % Assemble global force vector
    matrF = form_forceF(matrF);
    %%%%%%%% end of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %% find solution at the EACH frequency for the transformed matrices (transformation is needed to apply periodic BCs)
    delta = zeros(n_freq, size(T,2));
    for ifreq = 1:n_freq
        delta(ifreq,:) = ((T'*(matrK-(2*pi*load_freqs(ifreq))^2*matrM)*T)\(T'*matrF))';
    end
    % 
    % extract nodal displacements at each frequency
    nfreedofs = nfdof-length(dep_dofs);
    node_disp = zeros(n_freq, nm, mdof);
    node_amp = zeros(n_freq, nm);
    for i = 1:nm
        for j = 1:mdof
            if nf(i,j) ~=0
               if nf(i,j) <= nfreedofs 
                   % displacements at free dofs
                   for ifreq = 1:n_freq
                        node_disp(ifreq,i,j) = delta(ifreq,nf(i,j));
                   end
               else
                   % displacements at the top row of masses
                   ind = nf(i,j) - nfreedofs;
                   for ifreq=1:n_freq
                        node_disp(ifreq,i,j) = delta(ifreq,ind);
                   end
               end
            end
        end
        % evaluate the node amplitudes at EACH frequency
        for ifreq=1:n_freq
            node_amp(ifreq,i) = abs(sqrt(node_disp(ifreq,i,1)*node_disp(ifreq,i,1) + node_disp(ifreq,i,2)*node_disp(ifreq,i,2)));
        end
    end
    %
    % evaluate the normalized amplitude at EACH output at EACH frequency
    a_out = zeros(n_freq,n_outputs);
    for ifreq=1:n_freq
        max_amp_at_f = max(node_amp(ifreq,:));
        for iout = 1:n_outputs
             a_out(ifreq,iout) = node_amp(ifreq,outputs(iout))./max_amp_at_f;
        end
    end


%% plot the solution at each frequency together with the initial network configuration
fig_num = 120;
for ifreq = 1:n_freq
        fig_num = fig_num + 1;
        figure(fig_num);
        set(gca,'FontSize',16);
        % initialize plot. get a handle to graphic object
        hpl = plot(NaN, NaN); 
        axis([-1 nmx -1 nmy]);
        % draw springs
        for k = 1:ns
            mass1 = connec(k,1);
            mass2 = connec(k,2);
            x1 = geom(mass1,1); y1 = geom(mass1,2);
            x2 = geom(mass2,1); y2 = geom(mass2,2);
            x = [x1 x2];
            y = [y1 y2];
            line(x,y,'LineStyle',':','Color','k');
        end
        % draw masses
        line(geom(:,1), geom(:,2),'LineStyle','none','Marker','o','MarkerEdgeColor','k');
        % draw mass numbers
        for k=1:nm
            text(geom(k,1),geom(k,2),num2str(k));
        end
        % draw loaded masses
         for i = 1:nm
             for j = 1:mdof
                 if load(i,j) ~= 0
                   line(geom(i,1), geom(i,2),'LineStyle','none','MarkerSize',11,'Marker','o','MarkerFaceColor','k');
                 end
             end
         end
         % draw output masses
         for iout = 1:n_outputs
            line(geom(outputs(iout),1), geom(outputs(iout),2),'LineStyle','none','MarkerSize',15,'Marker','o','MarkerFaceColor','y');
         end
         % draw fixed masses
         for i = 1:nm
             for j = 1:mdof
                 if nf(i,j) == 0
                   line(geom(i,1), geom(i,2),'LineStyle','none','MarkerSize',10,'Marker','*','MarkerFaceColor','k');
                 end
             end
         end
         axis off;
        %%%%%%%%%%%%%%%%%
        %% draw the solution at a current frequency
         percentage = 0.2;
        % estimate the max amplitude
         max_amp = max(node_amp(ifreq,:));
        % draw springs
         for k = 1:ns
            mass1 = connec(k,1);
            mass2 = connec(k,2);
            x1 = geom(mass1,1)+percentage*node_disp(ifreq,mass1,1)/max_amp; 
            y1 = geom(mass1,2)+percentage*node_disp(ifreq,mass1,2)/max_amp;
            x2 = geom(mass2,1)+percentage*node_disp(ifreq,mass2,1)/max_amp; 
            y2 = geom(mass2,2)+percentage*node_disp(ifreq,mass2,2)/max_amp;
            x = [x1 x2];
            y = [y1 y2];
            line(x,y,'Color',[1 0 0]);
        end
        % draw masses
        for k=1:nm
            % define the color of a marker according to the displacement value:
            % blue = close to 0 displacement, red = close to red displacement
            colmass = node_amp(ifreq,k)./max_amp;
            line(geom(k,1)+percentage*node_disp(ifreq,k,1)./max_amp, geom(k,2)+percentage*node_disp(ifreq,k,2)./max_amp,'LineStyle','none',...
                 'Marker','o','MarkerFaceColor',[colmass 0 (1-colmass)],'MarkerEdgeColor',[colmass 0 (1-colmass)]);
        end
end

