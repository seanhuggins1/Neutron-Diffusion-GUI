% function neutron_ftcs(Params, Nspace, Ntime, tau_rel, varargin)
% Function performes a forward time-centered space integration of neutron
% diffusion in a bar of a given length
%
% Required Input:
% ===============
%
% Params (:) 3 element array of parameters [L C D] which correspond to the
%           bar length, the neutron creation rate, and the diffusion coefficient,
%           respectively
% Nspace (:) The number of spatial grid points
% Ntime (:) The number of time steps to perform the integration for
% tau_rel  (:) The relative time step in units of the critical time step
%              for scheme stability
%
% Output:
% =======
%
% No output
%
%
% Varargin:
% =========
%
% 'BCindex'	[0] set to any value that is a scalar multiple of the spatial
%               grid size and is within the domain of the bar
% 'Plot'	['3DMesh'] set to '3DMesh' to plot a 3D surface of the time
%                      evolution of neutron diffusion in a bar
%                      set to '2DColour' to plot a 2D colour plot of the time
%                      evolution of neutron diffusion in a bar
%                      set to 'snap' to plot snapshots of neutron density
%                      vs. position at every 10 time steps
%                      set to 'AvgN' to plot the average neutron density as
%                      a function of time
%
% Requires: no external m-files
% =========
%
% Example Use: neutron_ftcs([1 1 1], 61, 100, 0.5, 'plot', '2Dcolour', 'BCindex', 0.1)
% ============
%
% Author:
% =======
%
% SHuggins 10 Dec. 2018
%
function neutron_ftcs(Params, Nspace, Ntime, tau_rel, varargin)

    % If called with no arguments, echo a useage line. 
    if nargin == 0
       disp(' ')
       disp('neutron_ftcs(Params, Nspace, Ntime, tau_rel, [varargin])')
       disp(' ')
       return
    end

    %
    % Check that all varargin come in pairs.
    %
    if mod(length(varargin),2) ~= 0
      disp(' ')
      disp('Error: mis-match (odd number) of vargargin inputs')
      disp(' ')
      return
    end
    
    %Set defaults for varargins
    plotType = '3dmesh';
    BCindex = 0;
    
    for i=1:2:length(varargin)

        switch lower(varargin{i})
        
        %Check to see what type of plot the user wants
        case 'plot',
           plotType = lower(varargin{i+1});
        
        %Check to see if the user specified an index for the boundary
        %condition
        case 'bcindex',
           BCindex = varargin{i+1};

        otherwise
           disp(' ')
           disp(sprintf('WARNING: unknown varargin <%s> ignored',varargin{i}))
           disp(' ')
        end
    end
    %Grab constants from params array
    L = Params(1);  %Length
    C = Params(2);  %Neutron creation rate
    D = Params(3);  %Diffusion coefficient
    
    %Calculate the spatial grid size
    h = L/(Nspace-1);
    
    
    %Check validity of BCindex input
    if ~isnumeric(BCindex) || ~isreal(BCindex) || isinf(BCindex) || isnan(BCindex)
        disp(' ')
        disp('ERROR(neutron_ftcs): The index of the boundary condition must be a finite, real number')
        disp('Choosing default boundary condition index 0')
    end
    %Make sure that the boundary condition index is a scalar multiple of
    %the spatial grid size
    if floor(BCindex/h)~= BCindex/h
        disp(' ')
        disp('ERROR(neutron_ftcs): The index of the boundary condition must be a scalar multiple of the spatial grid size')
        disp('Choosing default boundary condition index 0')
        BCindex = 0;
    end
    %Check to see that the boundary condition index is inside the domain
    if (BCindex > L/2 || BCindex < -L/2)
        disp(' ')
        disp('ERROR(neutron_ftcs): The index of the boundary condition must be within the domain of the bar')
        disp('Choosing default boundary condition index 0')
        BCindex = 0;
    end
    
    %Check validity of time step input
    if ~isnumeric(tau_rel) || ~isreal(tau_rel) || isinf(tau_rel) || isnan(tau_rel) || tau_rel <= 0
        disp(' ')
        disp('ERROR(neutron_ftcs): The relative time step must be a finite, positive, real number greater than zero')
        disp('Choosing default relative time step 0.5')
        tau_rel = 0.5;
    end
     if tau_rel >= 1
        disp(' ')
        disp('WARNING(neutron_ftcs): A relative time step greater than or equal to unity will cause an unstable solution')
    end
    
    %Check validity of Ntime input
    if ~isnumeric(Ntime) || ~isreal(Ntime) || isinf(Ntime) || isnan(Ntime) || Ntime <= 0 || floor(Ntime)~= Ntime
        disp(' ')
        disp('ERROR(neutron_ftcs): The number of iterations must be a finite, positive integer greater than zero')
        disp('Choosing default number of iterations 100')
        Ntime = 100;
    end

    %Check validity of Nspace input
    if ~isnumeric(Nspace) || ~isreal(Nspace) || isinf(Nspace) || isnan(Nspace) || Nspace <= 0 || floor(Nspace)~= Nspace || (mod(Nspace,2) == 0)
        disp(' ')
        disp('ERROR(neutron_ftcs): The number of spatial grid points must be a finite, positive, odd integer greater than zero')
        disp('Choosing default number of spatial grid points 61')
        Nspace = 61;
    end
    
    %Check validity of Bar Length input
    if ~isnumeric(L) || ~isreal(L) || isinf(L) || isnan(L) || L <= 0
        disp(' ')
        disp('ERROR(neutron_ftcs): The bar length Params(1) must be a finite, positive, real number greater than zero')
        disp('Choosing default length 1')
        L = 1;
    end
    
    %Check validity of Bar Length input
    if ~isnumeric(L) || ~isreal(L) || isinf(L) || isnan(L) || L <= 0
        disp(' ')
        disp('ERROR(neutron_ftcs): The neutron creation rate Params(2) must be a finite, real number')
        disp('Choosing default diffusion neutron creation rate 1')
        C = 1;
    end
    
    %Check validity of Bar Length input
    if ~isnumeric(L) || ~isreal(L) || isinf(L) || isnan(L) || L <= 0
        disp(' ')
        disp('ERROR(neutron_ftcs): The diffusion coefficient Params(3) must be a finite, real number')
        disp('Choosing default diffusion coefficient 1')
        D = 1;
    end
    
    %Calculate the spatial grid size again, since L may have changed
    h = L/(Nspace-1);
       
    allowedPlotTypes = {'3dmesh','2dcolour','snap','avgn'};
    %Check validity of plot type input
    if ~any(strcmp(plotType,allowedPlotTypes)) 
        disp(' ')
        disp('ERROR(neutron_ftcs): The plot type was not recognized')
        disp('Choosing default plot type 3DMesh')
        plotType = '3dmesh';
    end
    
    
       
    %Specify our neutron density grid
    rho = zeros(Ntime,Nspace);
     
    %Specify the critical time step where the solution becomes unstable
    tau_crit = (h^2)/(2*D);
    
    %Calculate tau based on the units of tau_rel
    %tau_rel must be <= 1 for stable solutions
    tau = tau_rel*tau_crit;
    
    %Specify the critical length 
    L_crit = pi*sqrt(D/C);
    
    %Apply central delta function initial condition
    rho(1, (Nspace+1)/2 + BCindex/h) = 1/h;
    
    %Setup a constant to make things easier
    const = ((D*tau)/(h^2));
    
    %For all time
    for n = 1:(Ntime-1)
        %For space within the boundary
        for i = 2:(Nspace-1)
            %Perform forward time centred space integration
            rho(n+1,i) = rho(n,i) + C*tau*rho(n,i) + const*(rho(n,i+1) - 2*rho(n,i) + rho(n,i-1));
        end
    end
    
    %Clear the current axis
    cla;
    
    %If the user chooses to plot a 3D Mesh
    if strcmp(plotType,'3dmesh')
        %Plot a surface with the units on the axis corresponding to time
        %and physical size
        t = 1:Ntime;
        x = linspace(-L/2, L/2, Nspace);
        surf(x, t, rho(t,:))
        %Label the axes
        zlabel('neutron density per unit length');
        ylabel('time');
        xlabel('position');
        %Title the plot
        title('Neutron Diffusion 3D mesh');

    %If the user chooses to plot a 2D colour plot
    elseif strcmp(plotType,'2dcolour')
        %Plot a 2D colour plot with the units on the axis corresponding to time
        %and physical size
        t = 1:Ntime;
        x = linspace(-L/2, L/2, Nspace);
        pcolor(x, t, rho(t,:));
        %Make it pretty (I just like how this looks I apologize if it's
        %not a good data viz practice)
        shading('interp');
        %Make a colour bar for temperature
        bar = colorbar;
        %Label the colour bar and the axes
        ylabel(bar, 'neutron density per unit length')
        ylabel('time');
        xlabel('position');
        %Title the plot
        title('Neutron Diffusion 2D colour plot');
        
    %If the user chooses to plot a snapshot plot
    elseif strcmp(plotType,'snap') 
        %Make a snapshot plot of our neutron density with the units on the
        %axes corresponding to physical sizes
        x = linspace(-L/2, L/2, Nspace);
        %Plot a snapshot every 10 time steps
        t = 1:10:Ntime;
        plot(x, rho(t,:)); 
        %Label the axes
        xlabel('position');
        ylabel('neutron density per unit length');  
        %Title the plot
        title('snapshots of neutron density vs. position taken every 10 time steps');
        
        
    %If the user chooses to plot the average neutron density as a function
    %of time
    elseif strcmp(plotType,'avgn')
        %Compute the average neutron density as a function of time and plot
        %it
        n = 1:Ntime;
        %Initialize the average neutron density array
        rho_avg = zeros(length(n),1);
        for n = 1:Ntime
            rho_avg(n) = sum(rho(n,1:Nspace))/Nspace;
        end
        plot(1:Ntime, rho_avg);
        %Label the axes
        xlabel('t');
        ylabel('average neutron density');     
        %Title the plot
        title('average neutron density vs. time');
    end
    
end
