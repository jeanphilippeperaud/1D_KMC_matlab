% Jean-Philippe Peraud - December 3, 2015

%%% this script does the following:
% -loads silicon data from dataSi.txt
% -solves the LINEARIZED phonon transport problem between two walls (1D)
% with given prescribed temperatures T_l and T_r




%%%%% load the data
load 'dataSi.txt'; % Note: instead of having two TA branches in those data,
                   % we accounted for the double effect by multiplying the
                   % density of states by two for the TA mode.
% the data is organized as follows:
% - each line represents a phonon frequency
% - column 1 is the frequency
% - column 2 is the density of states
% - column 3 is the group velocity
% - column 4 is the width of the frequency cell (often a constant but not
% necessarily
% - column 5 is the relaxation time for this particular mode
% - column 6 is the polarization (1 if LA, 2 if TA). Not strictly needed
% for this code

% define reduced Planck constant and Boltzmann constant
hbar=1.054517e-34;
boltz=1.38065e-23;

% define linearization temperature (also referred to as "equilibrium" 
% temperature)
Teq = 300;

Tinit = 300; % initial temperature
Tl = 301; % left wall
Tr = 299; % right wall
L = 2e-7; % System length
N = 200000; % number of particles

NX = 20; % number of spatial cells
Ntt = 60; % number of measurement times
tmax = 3e-10; % maximum time

tt = 0:tmax/(Ntt-1):tmax; %times of measurement
Ntt = length(tt); % in some cases one may want to define tt differently than above. 
                  % Important for Ntt to be exactly the length of it
xx = L/2/NX:L/NX:L; % centroids of cells for measuring temperature
                    % not really needed but useful for plotting results
Dx = L/NX; % length of a cell
% dataSi(:,5) = 4e-11*ones(Nmodes,1);  %% <-- uncomment this for single relaxation time
SD = dataSi(:,2);  % Density of states
V = dataSi(:,3); % velocities
Dom = dataSi(:,4); % Delta frequencies
tau_inv = 1./dataSi(:,5); % relaxation rates
tau = dataSi(:,5); % relaxation times
F = dataSi(:,1); % frequencies
de_dT = (hbar*F/Teq).^2/boltz.*exp(hbar*F/boltz/Teq)./(exp(hbar*F/boltz/Teq)-1).^2; %derivative of Bose-Einstein

Nmodes = length(F); 

T  = zeros(Ntt,NX);
Qx = zeros(Ntt,NX);

% cumulative distribution functions
cumul_base = zeros(Nmodes,1);
cumul_coll = zeros(Nmodes,1);
cumul_vel  = zeros(Nmodes,1);
cumul_base(1) = SD(1)*de_dT(1)*Dom(1);
cumul_coll(1) = SD(1)*de_dT(1)*Dom(1)*tau_inv(1);
cumul_vel(1) = SD(1)*de_dT(1)*Dom(1)*V(1);
for i=2:Nmodes
    cumul_base(i) = cumul_base(i-1)+SD(i)*de_dT(i)*Dom(i);
    cumul_coll(i) = cumul_coll(i-1)+SD(i)*de_dT(i)*tau_inv(i)*Dom(i);
    cumul_vel(i) = cumul_vel(i-1) + SD(i)*de_dT(i)*V(i)*Dom(i);
end
C = cumul_base(Nmodes); % Heat capacity at Teq

%Deviational energy from the different sources
enrgInit = L*C*abs(Tinit-Teq);
enrgLeft = cumul_vel(Nmodes)*tt(Ntt)*abs(Tl-Teq)/4;
enrgRight = cumul_vel(Nmodes)*tt(Ntt)*abs(Tr-Teq)/4;

% Total deviational energy
enrgTot = enrgInit + enrgLeft + enrgRight;

% effective energy
Eeff = enrgTot/N;


% calculate thermal conductivity (optional. just to make sure parameters are realistic)
ktest = sum(SD.*Dom.*V.*V.*tau.*de_dT)/3;
tic
for i=1:N %loop over the N particles
    Ri = rand(); %draw random number to decide is particle is emitted from left or right wall
                % or from initial condition
    if Ri < enrgInit/enrgTot % this case: emitted from initial source
        x0 = L*rand();
        ind_mod = select_mode(cumul_base,Nmodes);
        R = 2*rand()-1; %this is cos(theta)
        vx = V(ind_mod)*R; %! not be confused between V and vx
        t0 = 0; % initial source => initial time of the particle is 0
        psign = sign(Tinit-Teq); % particle sign
    else
        if Ri > enrgInit/enrgTot && Ri < 1-enrgRight/enrgTot  % emitted by left wall
            x0 = 0;
            ind_mod = select_mode(cumul_vel,Nmodes);
            R = sqrt(rand()); %this is cos(theta)
            vx = V(ind_mod)*R; %! not be confused between V and vx
            psign = sign(Tl-Teq);
        else  % emitted by right wall
            x0 = L;
            ind_mod = select_mode(cumul_vel,Nmodes);
            R = -sqrt(rand()); %this is cos(theta)
            vx = V(ind_mod)*R; %! not be confused between V and vx
            psign = sign(Tr-Teq);
        end
        t0 = rand()*tt(Ntt); % emission by a wall => initial time of the particle 
                             % is anything between 0 and tmax
    end
    
    finished = false;  % as long as "false", the current particle is active
    im = 1; % index for tracking measurement times
    while tt(im)<t0
        im = im+1;
    end
    while ~finished
        Delta_t = -tau(ind_mod)*log(rand()); % time to next scattering event
        t1 = t0 + Delta_t; % time of next scattering event
        x1 = x0 + Delta_t*vx; % position of next scattering event
        
        % --- this part handles the contribution of the current particle to
        % the final estimates ------------------- %
        while (im<Ntt+1 && t0<=tt(im) && t1>tt(im))
            x_ = x0 + (tt(im)-t0)*vx; %% position at time tt(im)
            indx = floor(x_/L*NX)+1; %% index of the current spatial cell
            if (indx>0 && indx<NX+1)
                T(im,indx) = T(im,indx) + psign*Eeff/C/Dx; % temperature
                Qx(im,indx) = Qx(im,indx) + psign*Eeff*vx/Dx; % heat flux
            end
            im = im + 1;
        end
        % --------------------------------------- %
        
        % select post-collision mode
        ind_mod = select_mode(cumul_coll,Nmodes);
        
        % update particle parameters
        R = 2*rand()-1; % this is cos(theta)
        psi = 2*pi*rand();
        vx = V(ind_mod)*R;
        x0 = x1;
        t0 = t1;
        
        if (t0>tt(Ntt) || x0<0 || x0>L) % particle is terminated if it exits 
                                        % the system or if it overshoots the 
                                        % largest time of measurement
            finished = true;
        end
    end
end
toc

%% un-comment this for visualizing the solution
%movieLength = 5; % duration of movie in seconds
%figure();
%for i=1:Ntt
%    plot(xx,T(i,:) + Teq);
%    axis([0 L Teq+min(min(T)) Teq+max(max(T))]);
%    pause(movieLength/Ntt);
%end



