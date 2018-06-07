% %  dust_trajectory.m
% %   
% %  Jeffrey Walker,  2012-2014, latest major revision: January 3 2014.

function [t,q,x,y,vx,vy,RLd,ne,ni,V_time,E_xt,E_yt,B_t,f_ix,f_iy,...
    f_nx,f_ny,lambda_D,lambda_i,Kn_R0,P0,P1,Pg1,Vgrain]=...
    dust_trajectory(a,rho,r_initial,v_initial,species,ch_model,...
    profile_type,n0,Te0,Ti0,Z,P,alphm,cycles,points,filename,...
    particle_pusher)
    
% List of Dependencies (other m-files need in order for this program to run):
% improved_arctan
% profiles.m
% accumulate_charge.m
% charging_models.m
% (list of charging model modules)
%     - dust_bisection.m, dust_function_list.m
%   - oml.m
%   - oml_monoenergetic_ions.m
%   - kortshagen.m
%   - hutchinson
%   - Not written yet? - northrop.m, simple.m
% (list of particles pusher modules)
%   - iterative_pusher.m; all iterative pushers require trisolver.m
%   - sheath_iterative_pusher.m
%   - corotating_iterative_pusher.m
%   - boris_pusher.m
%   - corotating_boris_pusher.m
%   - sheath_boris_pusher.m

% explanation of inputs:
%	  species         =   atomic mass number of ion species
%	  ch_model        =   string which specifies charging model, current options
%                       (all listed in charging_models.m):
%                       'northrop','oml','simple','kortshagen','hutchinson',
%                       'oml_monoenergetic_ions', and 'constant_q'.
%	  profile_type	  =   string which sets up the inhomogenous plasma  
%                       profile of choice. IF THIS IS SOME KIND OF 
%                       CYLINDRICAL PROFILE, MAKE SURE THE STRING HAS *cyl* 
%                       SOMEWHERE IN THE LABEL YOU USE IN profiles.m
%	  n0              =   background density in m^-3
%	  Te0             =   electron thermal temperature, in eV
%	  Ti0             =   ion thermal temperature, in eV
%	  Z               =   ion charge state; generally Z=1. If the 
%                       'constant_q' profile is chosen, then Z denotes the 
%                       number of elementary charges on the dust grain.
%	  P               =   neutral gas pressure in Pascals
%   alphm           =   charge delay, where alphm<<1 means a high charge 
%                       delay; alphm<1 means the grain charges to the new 
%                       in-situ equilibrium charge more slowly
%   a   =   dust grain radius in meters. Alternatively, the user can input
%           the string 'electron', or 'ion', and dust_trajectory will
%           calculate the trajectory for an electron or ion respectively 
%           for the plasma profile specified by profiles.m
% 	rho	=   dust grain mass density in kg/m^-3; water is 1000 while  
%           melamine is 1574. Lunar Regolith is ~3000 kg/m^-3, see 
%           lunar_stratigraphy textbook for reference. Additionally, the
%           user can specify a string 'electron' or 'ion' to see resulting
%           particle trajectories for a given plasma profile specified in
%           profiles.m
%
%
%   particle_pusher     =   A string that specifies how you will push the
%                           dust grain. Allowable options: (remember to use
%                           single quotes for strings!) boris_pusher,
%                           corotating_boris_pusher, sheath_boris_pusher,
%                           iterative_pusher, sheath_iterative_pusher, and
%                           I still need to make a
%                           corotating_iterative_pusher as of Aug 29 2014
%                           (is corotating_iterative_pusher working now?)
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Explanation of outputs:
% Please Note that you can call dust_trajectory.m without specifying the
% outputs. E.g., you use:
% dust_trajectory(a,rho,r_initial,v_initial,species,ch_model,...
%   profile_type,n0,Te0,Ti0,Z,P,alphm,cycles,points,filename,...
%   particle_pusher);
%
% Instead of:
% [t,q,x,y,vx,vy,RLd,ne,ni,V_time,E_xt,E_yt,B_t,f_ix,f_iy,f_nx,f_ny,...
%   lambda_D,lambda_i,Kn_R0,P0,P1,Pg1,Vgrain]=dust_trajectory(a,rho,...
%   r_initial,v_initial,species,ch_model,profile_type,n0,Te0,Ti0,Z,P,...
%   alphm,cycles,points,filename,particle_pusher);
%
% The output variables will be saved to the target filename, so you don't
% need to keep the output variables in memory. 
% Here is the list of inputs with explanations:
%   t         = time in seconds
%   q         = charge of dust grain in coloumbs
%   x         = x position of dust grain in meters
%   y         = y position of dust grain in meters
%   vx        = velocity of dust grain in x-direction in m/s
%   vy        = velocity of dust grain in y-direction in m/s
%   RLd       = time-dependent Larmor radius of dust grain in m
%   ne        = local electron density of the plasma, sampled by grain, m^(-3)
%   ni        = local ion denisty of the plasma, sampled by the grain, m^(-3)
%   V_time    = local space potential sampled by the grain, in Volts
%   E_xt      = local x-component of electric field sampled by the grain, V/m
%   E_yt      = local y-component of electric field sampled by the grain, V/m
%   B_t       = local z-component of magnetic field sampled by the grain, in T
%   f_ix      = local x-component of ion velocity in m/s sampled by grain
%   f_iy      = local y-component of ion velocity in m/s sampled by grain
%   f_nx      = local x-component of neut velocity in m/s sampled by grain
%   f_ny      = local y-component of neut velocity in m/s sampled by grain
%   lambda_D  = local Debye length of plasma sampled by grain, in m
%   lambda_i  = local ion-neutral ch.exch. mean free path sampled by grain, in m
%   Kn_R0     = time-dependent Knudsen capture radius of grain, dimensionless
%   P0        = time-dep probability ion has no collisions in grain sheath
%   P1        = time-dep probability ion has 1 collision in grain sheath
%   Pg1       = time-dep probability ion has many collisions in grain sheath
%   Vgrain    = time-dependent difference of the surface potential of the 
%                     dust grain and local space potential, normalized by Te; 
%                     Vgrain = Vsurf - Vspace = q/(C Te) (Vspace=V_time here) 
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% some fundamental constants:
qe=1.6e-19;
me=9.1e-31;
mp=1.67e-27;
eps0=8.854e-12;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %	Error checking

% % do some error checking on the ion species specified by user
if rem(species,1)~=0 || species<=0 || species >300
	exception='species must be an integer value and less than 300';
	error(exception)
end
% % Now that we've checked that the species chosen is marginally valid, 
% % compute ion mass and neutral atom mass
mi = species*mp;
mr=me/mi;
% assumption that neutral gas atoms have the same mass as ions.
m_neut=mi;
% % do some error checking on the r_initial and v_initial arrays to make 
% % sure they are the right size, then load the arrays into the x0, y0, 
% % v_0x, and v_0y constants.
[m,n]=size(r_initial);
if m~=1 || n~=2
    tempchar = 'Initial position array is in the wrong dimensions.';
    tempchar = strcat(tempchar,' Please use r_initial=[xi yi].');
	exception = tempchar;
    clear tempchar;
    error(exception)
end
x0=r_initial(1);
y0=r_initial(2);

[m,n]=size(v_initial);
if m~=1 || n~=2
    tempchar='Initial velocity array is in the wrong dimensions.';
    tempchar=strcat(tempchar,' Please use v_initial=[vxi vyi].'); 
	exception=tempchar;
    clear tempchar;
	error(exception)
end
v_0x=v_initial(1);
v_0y=v_initial(2);

if points<0 || cycles<=0
	tempchar='Choose positive definite values for the number of';
    tempchar=strcat(tempchar,' gyrocycles and points/gyrocycle');
    exception=tempchar;
    clear tempchar;
	error(exception)
end
% check to make sure background plasma density, electron temperature, ion
% temperature, and adjustable charge rate parameter are all greater than
% zero.
if n0<=0 || Te0<=0 ||Ti0<=0 ||alphm<=0
	tempchar='Choose positive definite values for plasma density,';
    tempchar=strcat(tempchar,' ion/electron temperatures (in eV), and');  
    tempchar=strcat(tempchar,' adjustable charge rate parameter alpha,m.');
    exception=tempchar;
    clear tempchar;
	error(exception)
end
% % now that we now Ti0 is okay, set Tn=Ti0.
% May 2014: Maybe this should be set as an input??
Tn=Ti0;
% determine if a is input as a number (the usual state of affairs) or as a 
% string. If 'electron' or 'ion' are input for a, the code will run using
% the electron or ion mass as the dust grain mass. Hence, it can be used to
% produce electron/ion trajectories for a given profile.
if ischar(a)==0 
    if a<=0 || a>0.01
        tempchar='Choose positive definite values for the dust grain';
        tempchar=strcat(tempchar,' radius. Also, choose a<0.01 m');
        exception=tempchar;
        clear tempchar;
        error(exception)    
    end
    % if the user has chosen a floating value for a, or a numerical value
    % for a, and has not specified 'constant_q' as the charge model, then 
    % the ion charge number Z must be restricted to positive definite 
    % values, and probably less than 118, since that is the atomic number 
    % of the heaviest known element.
%     if (Z<=0 || Z>118 || floor(Z)~=Z) && ...
%             strcmp(ch_model,'constant_q')==0
%         tempchar='You must chose an integer number for the ion charge';
%         tempchar=strcat(tempchar,' state, which is positive definite,');
%         tempchar=strcat(tempchar,' greater than zero, and less than 118.');
%         exception=tempchar;
%         clear tempchar;
%         error(exception)   
%     % can I replace this end with an else?
%     end 
    % If we have gotten this far, then 'constant_q' has been chosen for the
    % charge model, and we must ensure that an integer value has been 
    % chosen for the charge state.
    if floor(Z)~=Z
        tempchar='You must chose an integer number for the number of';
        tempchar=strcat(tempchar,' elementary charges on the dust grain.');
        exception=tempchar;
        clear tempchar;
        error(exception);
    end
    % Since we have safely passed all of the tests, set the dust
    % temperature to the neutral temperature:
    Td=Tn;
else
	% user has input either 'ion' or 'electron' for the dust grain size
  	% a, so set the "dust grain mass" to the ion or electron mass. The
 	% code will continue along, calculating the trajectory for a test
  	% elementary particle (ion or electron).
  	if strcmp(a,'ion')==1 || strcmp(a,'electron')==1
       	if strcmp(a,'ion')
            % the mass we should use for the "dust mass" is the mass of an 
            % ion.
          	md=mi;
            % May as well set "dust temperature" to the ion temperature
            % here:
            Td=Ti0;
            % may as well set rho='ion' just in case user made a mistake
            rho='ion';
            if Z<=0 || Z>118 || floor(Z)~=Z
               tempchar='You must choose 118>Z>=1 for the ion charge';
               tempchar=strcat(tempchar,' state.');
               exception=tempchar;
               clear tempchar;
               error(exception)
            end
        end
        if strcmp(a,'electron')
            % the mass we should use for the "dust mass" is the mass of an
            % electron.
          	md=me;
            % May as well set "dust temperature" to the electron 
            % temperature here:
            Td=Ti0;
            % may as well set rho='electron' just in case user made a 
            % mistake
            rho='electron';
            if Z~=-1
               tempchar='The charge state for an electron is -1';
               tempchar=strcat(tempchar,' elementary charge.');
               exception=tempchar;
               clear tempchar;
               error(exception)
            end
        end
  	% user has input a string value for a, but the string is not 'ion' 
   	% or 'electron'.
    else
      	tempchar='If you are inputting a string for dust grain radius';
       	tempchar=strcat(tempchar,' a, then the only permissible');
       	tempchar=strcat(tempchar,' options are ion or electron');
       	exception=tempchar;
       	clear tempchar;
       	error(exception);
    end
end
% determine the dust grian mass density rho has been input as a string. If
% not, then ischar(rho)==0.
if ischar(rho)==0
    if rho<=0 || rho>3e7
        tempchar='Choose positive definite values for the dust grain';
        tempchar=strcat(tempchar,' density. Also, choose rho<30,000,000');
        tempchar=strcat(tempchar,'kg/m^3');
        exception= tempchar;
        clear tempchar;
        error(exception)
    end
    % Since we have passed the test above, we can now calculate the mass of
    % the dust grain, in kg m^-3
    md=rho*(4*pi*a^3/3);	
% ischar(rho)==1, which means that the user must specify either ion or
% electron.
else
    if strcmp(a,'ion')==1 || strcmp(a,'electron')==1
       	if strcmp(a,'ion')
            % set the ion radius equal to the bohr radius in meters for 
            % now; obviously this needs to be changed based on whatever 
            % species is chosen.
            a=5.29e-11;
            % use mi/(4/3*pi*a^3) for ion density.
            rho=mi/(4/3*pi*a.^3);
        end
        if strcmp(a,'electron')
            % set the electron radius equal to the classical electron
            % radius in meters:
            a=2.282e-15;
            % use me/(4/3*pi*a^3) for electron density.
            rho=me/(4/3*pi*a.^3);
        end
  	% user has input a string value for rho, but the string is not 'ion' 
   	% or 'electron'.
    else
      	tempchar='If you are inputting a string for dust grain density';
       	tempchar=strcat(tempchar,' rho, then the only permissible');
       	tempchar=strcat(tempchar,' options are ion or electron');
       	exception=tempchar;
       	clear tempchar;
       	error(exception);
    end
end
% Make sure the user does not input a negative value for the neutral gas
% pressure!
if P<0
	tempchar='Choose positive definite values for the neutral gas';
    tempchar=strcat(tempchar,' pressure (P>=0).');
    exception=tempchar;
    clear tempchar;
	error(exception)
end
% This segment of code may no longer be necessary.	
% if Z<=0 || Z>200
% 	tempchar='Choose a positive definite value for the ion charge state.';
%     tempchar=strcat(tempchar,' Also, Z is restricted to <200.');
%     exception=tempchar;
%     clear tempchar;
% 	error(exception)
% end
% % if a filename is not specified, then the output will be written to 
% % test.mat
if length(filename)==0 || ischar(filename)==0
	filename='test';
    tempchar='A valid filename was not chosen, so data will be written';
    tempchar=strcat(tempchar,' to test.mat');
    disp(tempchar);
    clear tempchar;
end
% % Check to see which particle mover the user has chosen.
if strcmp(particle_pusher,'boris_pusher')==0 && ...
        strcmp(particle_pusher,'iterative_pusher')==0 && ...
        strcmp(particle_pusher,'corotating_boris_pusher')==0 && ...
        strcmp(particle_pusher,'corotating_iterative_pusher')==0 && ...
        strcmp(particle_pusher,'sheath_boris_pusher')==0 && ...
        strcmp(particle_pusher,'sheath_iterative_pusher')==0
	% % current options are the boris pusher, and the iterative pusher. 
	% % If you want anything else, you'll have to build it!
	tempchar='Choose a valid way of time-advancing the dust grain. Your';
    tempchar=strcat(tempchar,' options are: boris_pusher,');
    tempchar=strcat(tempchar,' iterative_pusher, corotating_boris_pusher');
    tempchar=strcat(tempchar,', corotating_iterative_pusher');
    tempchar=strcat(tempchar,', sheath_iterative_pusher,');
    tempchar=strcat(tempchar,', or sheath_boris_pusher.');
    tempchar=strcat(tempchar,' Make sure to put single quotes around');
    tempchar=strcat(tempchar,' your choice of time-advancement pusher,');
    tempchar=strcat(tempchar,' i.e., it must be input as a string.');
    exception=tempchar;
    clear tempchar;
	error(exception)
end
	
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % A few quick definitions and constants:

% % Initial perp-temperature of the dust; vi should be comparable to this.
%Td=Ti0;   	% Td~Ti~Tn; unless Td~100 eV (this is measured in some 
            % experiments)
            
% compute the expected thermal speed of the test particle/grain, assuming a
% kinetic temperature of Td (determined during the error checking above).
vd=sqrt(2*qe*Td/md)         % the resulting velocity, based on Tdust.
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % Initial conditions for the plasma parameters and other derived
% % quantities; initialize some parameters and DEFINE THE TIMESTEP.

% % Make sure that there is at least one point per gyrocycle. Maybe this 
% % can be put up at the beginning with the other error checks?
if points<1
	points=1;
end
H=1/points;     % what fraction of a gyroperiod or bounce period do we 
                % increment by?
% % compute the number of timesteps needed to execute the number of 
% % gyro-cycles if we are at the equilibrium charge. Also depends on how 
% % many points we want per gyro-cycle.
nsteps=round(cycles*points);	% I added a rounding feature in case the 
                                % user selects a fraction of a gyrocycle.
                                
% % it is occaisionally useful to short circuit the main loop and only get 
% % one timestep; if that is the case and nsteps=0, then only compute 2 
% % timesteps.
if nsteps==0 || nsteps==1
	nsteps=2;
	disp('Only computing one time step.')
end

% estimate equilibrium charge Q0 and use it to define Newton timestep. 
% First, call profiles to get the right plasma conditions. Use t=0 in
% profile
[V_space,Ex,Ey,B,vix,viy,vex,vey,vnx,vny,gx,gy,ni0,ne0,alph,Ti0,Te0,...
    nneut,lambda_i0,lambda_D0,corot_period]=...
    profiles(Ti0,Te0,n0,0,x0,y0,profile_type,P,species);
% % Compute the capacitance after profiles!
C0=4*pi*eps0*a*(1+a/lambda_D0);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % some extra stuff added in here.
% Tr=Te0/Ti0;
% mr=me/mi;
% M=[0 0];
% eta=ne0/ni0;
% Kna=1e5;
% KnD=lambda_D0/a;
% R_le=sqrt(pi/4)*me*sqrt(2*qe*Te0/me)/qe/B;
% e_mag=a/R_le;
% mag_ratio=1;
% NDe=4/3*pi*n0.*lambda_D0.^3;
% % initialize counters and stuff
% Z=0;
% cnt=1;
% tchg=0;
% % find in situ equilibrium charge:
% while   Z<=0
%  	dZdt=dimensionless_charger(ch_model,Z,Tr,mr,M,eta,Kna,KnD,alph,...
%             e_mag/mag_ratio);
%         dt = 1/3*KnD/(1+1/KnD)/(1+Tr/eta)/NDe/abs(dZdt);
%         dZ=1/3*KnD/(1+1/KnD)/(1+Tr/eta)/NDe*sign(dZdt);
%         Zarr(cnt)=Z;
%         tarr(cnt)=tchg;
%        	if cnt>2&&Zarr(cnt)==Zarr(cnt-2) 
%             % break out of the loop when the charge begins oscillating back
%             % and forth between two values.
%         	break
%         end
%         cnt=cnt+1;
%         tchg=tchg+dt;
%         Z=Z+dZ; 
%         %disp(num2str(Z))
% %         drawnow;
% %         figure(1);clf;plot(tarr,Zarr)
%     
% end
% clear Zarr;clear tarr;
% Zeq=Z
% Q0_different=Zeq*C0*Te0/qe;
% disp(strcat('Number of charges using different method:',num2str(Q0_different)));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% % Compute the initial grain speed with respect to electrons and ions; 
% % needed when the dust grain is moving in a given charge model. If the
% % simulation is done in a sheath, then the following workaround is used.
if strcmp(particle_pusher,'sheath_boris_pusher') ||...
        strcmp(particle_pusher,'sheath_iterative_pusher')
    % use vey to store information about the flow velocity of electrons in
    % the phi direction, i.e., vey=(r d/dt phi)_e
    phi=improved_arctan(x0,y0);
    we = sqrt((v_0x-(vey)*sin(phi)).^2+(v_0y+(vey)*cos(phi)).^2);
    % use vex to store information about the flow velocity of ions in the z
    % direction
    wi=sqrt((v_0x-vix).^2+(v_0y-viy).^2 + (vex).^2);
    w=[we wi];
else
    % the regular, non-convoluted method which works for everything else
    we=sqrt((v_0x-vex).^2+(v_0y-vey).^2);
    wi=sqrt((v_0x-vix).^2+(v_0y-viy).^2);
    w=[we wi];
end
% % first argument of charging_models is qflag variable; if it is 1,  
% % then you are calculating equilibrium charge, but if it is 0 then do not  
% % calculate equilibrium charge; base the timestep off of this.
% % Get predicted equilibrium charge, starting with a grain charge of 0 for 
% % simplicity.
[Itoti,Q0,Kn_R0i,P0i,P1i,Pg1i]=...
    charging_models(1,ch_model,a,alph,Te0,Ti0,ne0,ni0,B,Z,C0,0,...
    lambda_D0,lambda_i0,w,species);	
% % display the equilibrium charge
number_of_charges=(Q0/qe);
disp(strcat('number of charges:',num2str(number_of_charges)));
% % The initial gyro frequency:
w_cdi=abs(Q0*B/md)



% % August 2013: may want to base the timestep off of the dust-neutral 
% % collision frequency if that is a larger quantity than the dust 
% % gyro-frequency. Additionally, if a gravitational field is present, one
% % also needs to pick a relevant timescale for gravity as well.

cn=sqrt(8*qe*Tn/pi/m_neut);
vxdriftn=v_0x-vnx;
vydriftn=v_0y-vny;
vn=sqrt(vxdriftn.^2+vydriftn.^2);
% % Using the Knudsen number, decide to use hyd. or kinetic regime for 
% % dust-neutral drag. Figure out whether to use Epstein drag force, or if 
% % it should go as velocity^2.
if vn<cn
	
	delta=1.26;		% coefficient for melamine??
	% % dust-neutral collision frequency for epstein drag force
	nu_dn=delta*nneut*cn*m_neut/a/rho
else
    delta=1.26;
	% Not sure about what to do if w>cn; how do I characterize the neutral 
	% dust collision frequency if the drag force is dependent on v^2 
    % instead of v?
	nu_dn=delta*nneut*cn*m_neut/a/rho	
end

% % determine which is the larger quantity: neutral-dust collision 
% % frequency, or gyro frequency.
% % Base the time step off of whichever quantity is largest. 
% % But first: check to see if both frequecies are zero!
if w_cdi==0 && nu_dn==0
	% % This section of code is for error checking; if B=0, also check to 
    % % see if Ex or Ey is zero.
	disp('Warning! No background magnetic field or neutral gas!') 
	% % right below this comment is where the electric field should be 
	% % checked. it may be possible to define some frequency which relies 
    % % on |E|.
	
	% % attempt to find a bounce period in order to specify a time step.
	dtNwt=0.001;	% Just a guess right now
end
% % Now check which is the largest freqency, and choose the smallest 
% % timescale for the timestep. If all quantities are zero, you will bypass 
% % this loop and dtNwt will be fixed according to the if statement above.
if w_cdi>nu_dn
    if corot_period == 0
        dtNwt=H*2*pi/w_cdi
    else
        if w_cdi>2*pi/corot_period
            dtNwt=H*2*pi/w_cdi 
        else
            dtNwt=H*corot_period
        end
    end
else 
    if corot_period == 0
        % for this case, check to see if B=0; this next set of statements
        % is used for small oscillations in a sheath. This requires that
        % ion flow must be non-zero. Dust charge must also be non-zero. See
        % notebook #7, page 53.
        if B==0 && w(2)~=0 && Q0~=0
            tempchar1='Using small oscillation frequency.';
            tempchar2=' Make sure grain starts slightly below the sheath';
            tempchar3=' boundary.';
            disp(strcat(tempchar1,tempchar2,tempchar3));
            %dtNwt=H*2*pi/nu_dn/50
            % use profiles again to find out what is the potential
            % difference between the plasma and the lower electrode:
            [V0,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~]=...
                profiles(Ti0,Te0,n0,0,0,0,profile_type,P,species);
            % small oscillation frequency around a stable levitation
            % height:
            w_small=sqrt(abs(Q0/md))*sqrt(qe*ni0*abs(viy)/...
                eps0*sqrt(mi/2/qe))/(-1.5*sqrt(qe*ni0*abs(viy)*...
                sqrt(mi/2/qe))*y0+(-V0)^(3/4))^(1/3)
            dtNwt=H*2*pi/w_small;
        else
            dtNwt=H*2*pi/nu_dn
        end
    else
        if nu_dn>2*pi/corot_period
            dtNwt=H*2*pi/nu_dn
        else
            dtNwt=H*corot_period
        end
    end
end

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % NOW initialize or declare arrays, after nsteps has been defined based  
% % on the user inputs cycles*points, which leads to arrays of size nsteps:
t=zeros(1,nsteps);
q=zeros(size(t));
Itot=zeros(size(t));
x=zeros(size(t));
y=zeros(size(t));
vx=zeros(size(t));
vy=zeros(size(t));
C=zeros(size(t));
% % profile variables:
ne=zeros(size(t));
ni=zeros(size(t));
Te=zeros(size(t));
Ti=zeros(size(t));
Tn=zeros(size(t));
E_xt=zeros(size(t));
E_yt=zeros(size(t));
B_t=zeros(size(t));
v_ex=zeros(size(t));
v_ey=zeros(size(t));
v_ix=zeros(size(t));
v_iy=zeros(size(t));
v_nx=zeros(size(t));
v_ny=zeros(size(t));
f_ix=zeros(size(t));
f_iy=zeros(size(t));
f_nx=zeros(size(t));
f_ny=zeros(size(t));
lambda_D=zeros(size(t));
lambda_i=zeros(size(t));
V_time=(size(t));
Vgrain=(size(t));
Kn_R0=(size(t));
P0=(size(t));
P1=(size(t));
Pg1=(size(t));
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % Initial conditions for the plasma parameters and other derived
% % quantities.
% % Intialize arrays with their starting values.
t(1)=0;
cntarr(1)=1;
x(1)=x0;
% % set x(1)=v_0y*md/Q0/B to make sure the dust grain goes through "both 
% % sides" of the inhomogeneity.
%x(1)=-v_0y*md/Q0/B;
y(1)=y0;
vx(1)=v_0x;
vy(1)=v_0y;
% Make q(1)=equilibrium charge if desired, to prevent transient
q(1)=Q0; 
% otherwise, use q(1)=0 or q(1)=some other value, but this will produce a 
% transient. Sometimes, this is particularly useful, other times not as
% much. q(1)=0 is quite useful for the case of an initially uncharged grain
% becoming charged after falling into plasma, such as in the case of an ice
% crystal being launched from Enceladus.
%q(1)=0;

% initialize capacitance and other quantities.
C(1)=C0;
Itot(1)=Itoti;	 
E_xt(1)=Ex;
E_yt(1)=Ey;
B_t(1)=B;
V_time(1)=V_space;
Vgrain(1)=q(1)/C0/Te0;
ni(1)=ni0;
ne(1)=ne0;
Ti(1)=Ti0;
Te(1)=Te0;
%Tn(1)=Ti0;
Tn(:)=Ti0;
v_ex(1)=vex;
v_ey(1)=vey;
v_ix(1)=vix;
v_iy(1)=viy;
v_nx(1)=vnx;
v_ny(1)=vny;
lambda_D(1)=lambda_D0;
lambda_i(1)=lambda_i0;
Kn_R0(1)=Kn_R0i;
P0(1)=P0i;
P1(1)=P1i;
Pg1(1)=Pg1i;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % advance the particle backward in time by half a timestep. this gets the 
% % velocities at t=-dtNwt/2, which is necessary for any leapfrog style 
% % method.
switch particle_pusher
	case 'iterative_pusher'
		%dtNwt=H*2*pi/abs(Q0);  
		%dtNwt=dtNwt*md/abs(B)
		[xtemp,ytemp,vx(1),vy(1),w]=...
            iterative_pusher(-dtNwt/2,a,rho,q(1),x(1),y(1),vx(1),vy(1),...
            species,E_xt(1),E_yt(1),B_t(1),gx,gy,ne(1),ni(1),nneut,...
            vnx,vny,vex,vey,vix,viy,Te(1),Ti(1),lambda_D(1),ch_model);
	
    case 'sheath_iterative_pusher'
		%dtNwt=H*2*pi/abs(Q0);  
		%dtNwt=dtNwt*md/abs(B)
		[xtemp,ytemp,vx(1),vy(1),w]=...
            sheath_iterative_pusher(-dtNwt/2,a,rho,q(1),x(1),y(1),vx(1),...
            vy(1),species,E_xt(1),E_yt(1),B_t(1),gx,gy,ne(1),ni(1),...
            nneut,vnx,vny,vex,vey,vix,viy,Te(1),Ti(1),lambda_D(1),...
            ch_model);
        
    case 'corotating_iterative_pusher'
		[xtemp,ytemp,vx(1),vy(1),w]=...
            corotating_iterative_pusher(-dtNwt/2,a,rho,q(1),x(1),y(1),vx(1),...
            vy(1),species,E_xt(1),E_yt(1),B_t(1),gx,gy,ne(1),ni(1),...
            nneut,vnx,vny,vex,vey,vix,viy,Te(1),Ti(1),Tn(1),lambda_D(1),...
            ch_model,corot_period,nu_dn);
        
    case 'boris_pusher'
		[xtemp,ytemp,vx(1),vy(1),w]=...
            boris_pusher(-dtNwt/2,md,q(1),x(1),y(1),vx(1),vy(1),...
            E_xt(1),E_yt(1),B_t(1),nu_dn,gx,gy,vex,vey,vix,viy,vnx,vny);
        
    case 'corotating_boris_pusher'
        [xtemp,ytemp,vx(1),vy(1),w]=...
            corotating_boris_pusher(-dtNwt/2,md,q(1),x(1),y(1),...
            vx(1),vy(1),E_xt(1),E_yt(1),B_t(1),nu_dn,gx,gy,vex,vey,...
            vix,viy,vnx,vny,corot_period);
        
    case 'sheath_boris_pusher'
       	[xtemp,ytemp,vx(1),vy(1),w]=...
         	sheath_boris_pusher(-dtNwt/2,md,q(1),x(1),y(1),...
          	vx(1),vy(1),E_xt(1),E_yt(1),B_t(1),nu_dn,...
           	gx,gy,vex,vey,vix,viy,vnx,vny);
end 
% % xtemp and ytemp are not needed when we go back just a half step.
% % you will get the wrong answer if you try to get x(1) and y(1) from the 
% % command written above.
clear xtemp; clear ytemp;

% % Now, all the initial values have been assigned, time to start the main
% % loop.
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % The Main loop.
t_acc=0;    % the accumulated time since a charge update starts at zero.
for k=2:nsteps
	%tic;
	
   % % display progress.
    if ~mod(k-2,nsteps/20)
		disp([num2str((k-2)/nsteps*100,'%.2f') '%'])
    end
   
    % diagnostic stuff
%     r_sat = 60268000;
%     drawnow;
%     subplot(1,2,1);
%     plot(x,y,'.b');hold on;
%     phi_sat = linspace(0,2*pi,1e2);
%     x_sat = r_sat*cos(phi_sat);
%     y_sat = r_sat*sin(phi_sat);
%     xlim([-10*r_sat 10*r_sat]);
%     ylim([-10*r_sat 10*r_sat]);
%     plot(x_sat,y_sat,'--r');
%     subplot(1,2,2);
%     plot(t,q/qe,'.g');hold on;
    %pause
    %q(k-1)
    
    %%~~~~#1 and #2:
    % % Use iterative leapfrog method, sheath iterative  
    % % leapfrog method, Boris method, corotating Boris method, 
    % % sheath Boris method to push the dust grain.
    % %   	
    switch particle_pusher
    	case 'iterative_pusher'
                                                   
            [x(k),y(k),vx(k),vy(k),w]=...
                iterative_pusher(dtNwt,a,rho,q(k-1),x(k-1),y(k-1),...
                vx(k-1),vy(k-1),species,E_xt(k-1),E_yt(k-1),B_t(k-1),...
                gx,gy,ne(k-1),ni(k-1),nneut,vnx,vny,vex,vey,vix,viy,...
                Te(k-1),Ti(k-1),lambda_D(k-1),ch_model);
            
            % % if using the iterative_pusher:
    		t(k)=t(k-1)+dtNwt;
    	%%~~~~#3    
    	% % final step in leapfrog method: calculate qd(t), E(t), and ion 
    	% % drag force. These quantities are computed at the full timestep 
    	% % Use profiles to calculate E(t) and E(space), and update plasma 
    	% % conditions. (they are spatially dependent.) NOTE: should 
    	% % profiles be called before or after chargeup?? 
        % % 7/16/2013 -->BEFORE.
            [V_time(k),E_xt(k),E_yt(k),B_t(k),vix,viy,vex,vey,vnx,vny,...
                gx,gy,ni(k),ne(k),alph,Ti(k),Te(k),nneut,lambda_i(k),...
                lambda_D(k),corot_period]=...
                profiles(Ti0,Te0,n0,t(k),x(k),y(k),profile_type,P,species);
        %figure(1);drawnow;plot(t,q/1.6e-19);%plot(x,y);%
   	
            % Compute the capacitance after profiles!
            C(k)=4*pi*eps0*a*(1+a/lambda_D(k));	
         	
            [q(k),Itot(k),Kn_R0(k),P0(k),P1(k),Pg1(k),t_acc]=...
                accumulate_charge(0,ch_model,a,alph,Te(k),Ti(k),...
                ne(k),ni(k),B_t(k),Z,C(k),q(k-1),dtNwt,alphm,...
                lambda_D(k),lambda_i(k),w,t_acc,species);
            
            clear cntarr;
            % Vgrain is the difference between the grain surface potential 
            % and the local space potential, normalized to the electron
            % temperature (Te is in volts.)
    		Vgrain(k)=q(k)/C(k)/Te(k);	
            
          	% store the electron, ion, and neutral flows for each time
            % step.
            v_ex(k)=vex;
            v_ey(k)=vey;
            v_ix(k)=vix;
            v_iy(k)=viy;
            v_nx(k)=vnx;
            v_ny(k)=vny;
    
    	% % put the 1/2 velocity update timestep here if you want to 
    	% % calculate vx and vy at the same timesteps as the positions. 
    	% % The time axis would no longer be staggered. Essentially, 
    	% % you are rotating and  accelerating half a timestep
    
    	case 'sheath_iterative_pusher'
                                                   
            [x(k),y(k),vx(k),vy(k),w]=...
                sheath_iterative_pusher(dtNwt,a,rho,q(k-1),x(k-1),y(k-1),...
                vx(k-1),vy(k-1),species,E_xt(k-1),E_yt(k-1),B_t(k-1),...
                gx,gy,ne(k-1),ni(k-1),nneut,vnx,vny,vex,vey,vix,viy,...
                Te(k-1),Ti(k-1),lambda_D(k-1),ch_model);
            
            % % if using the iterative_pusher:
    		t(k)=t(k-1)+dtNwt;
    	%%~~~~#3    
    	% % final step in leapfrog method: calculate qd(t), E(t), and ion 
    	% % drag force. These quantities are computed at the full timestep 
    	% % Use profiles to calculate E(t) and E(space), and update plasma 
    	% % conditions. (they are spatially dependent.) NOTE: should 
    	% % profiles be called before or after chargeup?? 
        % % 7/16/2013 -->BEFORE.
            [V_time(k),E_xt(k),E_yt(k),B_t(k),vix,viy,vex,vey,vnx,vny,...
                gx,gy,ni(k),ne(k),alph,Ti(k),Te(k),nneut,lambda_i(k),...
                lambda_D(k),corot_period]=...
                profiles(Ti0,Te0,n0,t(k),x(k),y(k),profile_type,P,species);
        %figure(1);drawnow;plot(t,q/1.6e-19);%plot(x,y);%
   	
            % Compute the capacitance after profiles!
            C(k)=4*pi*eps0*a*(1+a/lambda_D(k));	
         	
            [q(k),Itot(k),Kn_R0(k),P0(k),P1(k),Pg1(k),t_acc]=...
                accumulate_charge(0,ch_model,a,alph,Te(k),Ti(k),...
                ne(k),ni(k),B_t(k),Z,C(k),q(k-1),dtNwt,alphm,...
                lambda_D(k),lambda_i(k),w,t_acc,species);
            
            clear cntarr;
            % Vgrain is the difference between the grain surface potential 
            % and the local space potential, normalized to the electron
            % temperature (Te is in volts.)
    		Vgrain(k)=q(k)/C(k)/Te(k);	
            
          	% store the electron, ion, and neutral flows for each time
            % step.
            v_ex(k)=vex;
            v_ey(k)=vey;
            v_ix(k)=vix;
            v_iy(k)=viy;
            v_nx(k)=vnx;
            v_ny(k)=vny;
    
    	% % put the 1/2 velocity update timestep here if you want to 
    	% % calculate vx and vy at the same timesteps as the positions. 
    	% % The time axis would no longer be staggered. Essentially, 
    	% % you are rotating and  accelerating half a timestep
    		
    	case 'corotating_iterative_pusher'
                                                   
            [x(k),y(k),vx(k),vy(k),w]=...
                corotating_iterative_pusher(dtNwt,a,rho,q(k-1),x(k-1),...
                y(k-1),vx(k-1),vy(k-1),species,E_xt(k-1),E_yt(k-1),...
                B_t(k-1),gx,gy,ne(k-1),ni(k-1),nneut,vnx,vny,vex,vey,...
                vix,viy,Te(k-1),Ti(k-1),Tn(k-1),lambda_D(k-1),ch_model,...
                corot_period,nu_dn);
            
            % % if using the iterative_pusher:
    		t(k)=t(k-1)+dtNwt;
    	%%~~~~#3    
    	% % final step in leapfrog method: calculate qd(t), E(t), and ion 
    	% % drag force. These quantities are computed at the full timestep 
    	% % Use profiles to calculate E(t) and E(space), and update plasma 
    	% % conditions. (they are spatially dependent.) NOTE: should 
    	% % profiles be called before or after chargeup?? 
        % % 7/16/2013 -->BEFORE.
            [V_time(k),E_xt(k),E_yt(k),B_t(k),vix,viy,vex,vey,vnx,vny,...
                gx,gy,ni(k),ne(k),alph,Ti(k),Te(k),nneut,lambda_i(k),...
                lambda_D(k),corot_period]=...
                profiles(Ti0,Te0,n0,t(k),x(k),y(k),profile_type,P,species);
        %figure(1);drawnow;plot(t,q/1.6e-19);%plot(x,y);%
   	
            % Compute the capacitance after profiles!
            C(k)=4*pi*eps0*a*(1+a/lambda_D(k));	
         	
            [q(k),Itot(k),Kn_R0(k),P0(k),P1(k),Pg1(k),t_acc]=...
                accumulate_charge(0,ch_model,a,alph,Te(k),Ti(k),...
                ne(k),ni(k),B_t(k),Z,C(k),q(k-1),dtNwt,alphm,...
                lambda_D(k),lambda_i(k),w,t_acc,species);
            
            clear cntarr;
            % Vgrain is the difference between the grain surface potential 
            % and the local space potential, normalized to the electron
            % temperature (Te is in volts.)
    		Vgrain(k)=q(k)/C(k)/Te(k);	
            
          	% store the electron, ion, and neutral flows for each time
            % step.
            v_ex(k)=vex;
            v_ey(k)=vey;
            v_ix(k)=vix;
            v_iy(k)=viy;
            v_nx(k)=vnx;
            v_ny(k)=vny;
    
    	% % put the 1/2 velocity update timestep here if you want to 
    	% % calculate vx and vy at the same timesteps as the positions. 
    	% % The time axis would no longer be staggered. Essentially, 
    	% % you are rotating and  accelerating half a timestep
        
    	case 'boris_pusher'
            [x(k),y(k),vx(k),vy(k),w]=...
                boris_pusher(dtNwt,md,q(k-1),x(k-1),y(k-1),...
                vx(k-1),vy(k-1),E_xt(k-1),E_yt(k-1),B_t(k-1),nu_dn,...
                gx,gy,vex,vey,vix,viy,vnx,vny);
            % % update the time array.
            t(k)=t(k-1)+dtNwt;
         
    	%%~~~~#3    
    	% % final step in leapfrog method: calculate qd(t), E(t), and ion 
    	% % drag force. These quantities are computed at the full timestep 
    	% % (they are spatially dependent.) Use profiles to calculate E(t) 
    	% % and E(space), and update plasma conditions. NOTE: should 
        % % profiles be called before or after chargeup?? 
    	% % 7/16/2013 -->BEFORE.
            [V_time(k),E_xt(k),E_yt(k),B_t(k),vix,viy,vex,vey,vnx,vny,...
                gx,gy,ni(k),ne(k),alph,Ti(k),Te(k),nneut,lambda_i(k),...
                lambda_D(k),corot_period]=...
                profiles(Ti0,Te0,n0,t(k),x(k),y(k),profile_type,P,species);
        
         	% Compute the capacitance after profiles!
            C(k)=4*pi*eps0*a*(1+a/lambda_D(k));		
     
            [q(k),Itot(k),Kn_R0(k),P0(k),P1(k),Pg1(k),t_acc]=...
                accumulate_charge(0,ch_model,a,alph,Te(k),Ti(k),...
                ne(k),ni(k),B_t(k),Z,C(k),q(k-1),dtNwt,alphm,...
                lambda_D(k),lambda_i(k),w,t_acc,species);
            
            clear cntarr;
            
            % Vgrain is the difference between the grain surface potential 
            % and the local space potential, normalized to the electron
            % temperature (Te is in volts.)
    		Vgrain(k)=q(k)/C(k)/Te(k);	
            
            % store the electron, ion, and neutral flows for each time
            % step.
            v_ex(k)=vex;
            v_ey(k)=vey;
            v_ix(k)=vix;
            v_iy(k)=viy;
            v_nx(k)=vnx;
            v_ny(k)=vny;
            
    	% % For ion and neutral drag, you have to assume vions-vdust>>vthi,
        % % and v_neut-vdust>>vthn in order for Boris method to work.
    	
        % % If you want to compute the drag forces as functions of time,
        % % put them here.
    	% % I think the ion drag should actually use q(k) instead of 
    	% % q(k-1), because q should be evaluated at a spatial location - 
        % % 5/23/2013
        
            %[f_ix(k),f_iy(k)]=ion_drag(q(k),a,Te(k),Ti(k),ni(k),C,lambda_D(k),lambda_i(k),vix,viy,0,0,ch_model,Kn_R0(k),P0(k),P1(k),Pg1(k),Z,B,gx,gy,species);
            %[f_nx(k),f_ny(k),nu_dn]=neutral_drag(a,rho,nneut,0,0,vx(k-1),vy(k-1),Tn,gx,gy,species);
   	
    	% % put the 1/2 velocity update timestep here if you want to calculate
    	% % vx and vy at the same timesteps as the positions. The time axis
    	% % would no longer be staggered. Essentially, you are rotating and 
    	% % accelerating half a timestep
 
    	
        
        % use corotating boris pusher when you are dealing with dust and a
        % planet or moon, or other celestial body.
        case 'corotating_boris_pusher'
            [x(k),y(k),vx(k),vy(k),w]=...
                corotating_boris_pusher(dtNwt,md,q(k-1),x(k-1),y(k-1),...
                vx(k-1),vy(k-1),E_xt(k-1),E_yt(k-1),B_t(k-1),nu_dn,...
                gx,gy,vex,vey,vix,viy,vnx,vny,corot_period);

            % % update the time array.
            t(k)=t(k-1)+dtNwt;
    	%%~~~~#3    
    	% % final step in leapfrog method: calculate qd(t), E(t), and ion 
    	% % drag force. These quantities are computed at the full timestep 
    	% % (they are spatially dependent.) Use profiles to calculate E(t) 
    	% % and E(space), and update plasma conditions. NOTE: should 
        % % profiles be called before or after chargeup?? 
    	% % 7/16/2013 -->BEFORE.
            [V_time(k),E_xt(k),E_yt(k),B_t(k),vix,viy,vex,vey,vnx,vny,...
                gx,gy,ni(k),ne(k),alph,Ti(k),Te(k),nneut,lambda_i(k),...
                lambda_D(k),corot_period]=...
                profiles(Ti0,Te0,n0,t(k),x(k),y(k),profile_type,P,species);

            % Compute the capacitance after profiles!
            C(k)=4*pi*eps0*a*(1+a/lambda_D(k));
    	
     
            [q(k),Itot(k),Kn_R0(k),P0(k),P1(k),Pg1(k),t_acc]=...
                accumulate_charge(0,ch_model,a,alph,Te(k),Ti(k),...
                ne(k),ni(k),B_t(k),Z,C(k),q(k-1),dtNwt,alphm,...
                lambda_D(k),lambda_i(k),w,t_acc,species);
            
            clear cntarr;
            
            % Vgrain is the difference between the grain surface potential 
            % and the local space potential, normalized to the electron
            % temperature (Te is in volts.)
    		Vgrain(k)=q(k)/C(k)/Te(k);	
            
            % store the electron, ion, and neutral flows for each time
            % step.
            v_ex(k)=vex;
            v_ey(k)=vey;
            v_ix(k)=vix;
            v_iy(k)=viy;
            v_nx(k)=vnx;
            v_ny(k)=vny;
            
    	% % For ion and neutral drag, you have to assume vions-vdust>>vthi,
        % % and v_neut-vdust>>vthn in order for Boris method to work.
    	
        % % If you want to compute the drag forces as functions of time,
        % % put them here.
    	% % I think the ion drag should actually use q(k) instead of 
    	% % q(k-1), because q should be evaluated at a spatial location - 
        % % 5/23/2013
        
            %[f_ix(k),f_iy(k)]=ion_drag(q(k),a,Te(k),Ti(k),ni(k),C,lambda_D(k),lambda_i(k),vix,viy,0,0,ch_model,Kn_R0(k),P0(k),P1(k),Pg1(k),Z,B,gx,gy,species);
            %[f_nx(k),f_ny(k),nu_dn]=neutral_drag(a,rho,nneut,0,0,vx(k-1),vy(k-1),Tn,gx,gy,species);
   	
    	% % put the 1/2 velocity update timestep here if you want to 
    	% % calculate vx and vy at the same timesteps as the positions. 
    	% % The time axis would no longer be staggered. Essentially, you  
    	% % are rotating and accelerating half a timestep
        
        % January 2014: use sheath_boris_pusher if you want work with dust
        % grains levitated in a planar sheath.
    	case 'sheath_boris_pusher'
            [x(k),y(k),vx(k),vy(k),w]=...
                sheath_boris_pusher(dtNwt,md,q(k-1),x(k-1),y(k-1),...
                vx(k-1),vy(k-1),E_xt(k-1),E_yt(k-1),B_t(k-1),nu_dn,...
                gx,gy,vex,vey,vix,viy,vnx,vny);
            % % update the time array.
            t(k)=t(k-1)+dtNwt;
         
    	%%~~~~#3    
    	% % final step in leapfrog method: calculate qd(t), E(t), and ion 
    	% % drag force. These quantities are computed at the full timestep 
    	% % (they are spatially dependent.) Use profiles to calculate E(t) 
    	% % and E(space), and update plasma conditions. NOTE: should 
        % % profiles be called before or after chargeup?? 
    	% % 7/16/2013 -->BEFORE.
            [V_time(k),E_xt(k),E_yt(k),B_t(k),vix,viy,vex,vey,vnx,vny,...
                gx,gy,ni(k),ne(k),alph,Ti(k),Te(k),nneut,lambda_i(k),...
                lambda_D(k),corot_period]=...
                profiles(Ti0,Te0,n0,t(k),x(k),y(k),profile_type,P,species);
        
            % Compute the capacitance after profiles!
            C(k)=4*pi*eps0*a*(1+a/lambda_D(k));
     
            [q(k),Itot(k),Kn_R0(k),P0(k),P1(k),Pg1(k),t_acc]=...
                accumulate_charge(0,ch_model,a,alph,Te(k),Ti(k),...
                ne(k),ni(k),B_t(k),Z,C(k),q(k-1),dtNwt,alphm,...
                lambda_D(k),lambda_i(k),w,t_acc,species);
           
            clear cntarr;
            
            % Vgrain is the difference between the grain surface potential 
            % and the local space potential, normalized to the electron
            % temperature (Te is in volts.)
    		Vgrain(k)=q(k)/C(k)/Te(k);	
            
            % store the electron, ion, and neutral flows for each time
            % step.
            v_ex(k)=vex;
            v_ey(k)=vey;
            v_ix(k)=vix;
            v_iy(k)=viy;
            v_nx(k)=vnx;
            v_ny(k)=vny;
            
    	% % For ion and neutral drag, you have to assume vions-vdust>>vthi,
        % % and v_neut-vdust>>vthn in order for Boris method to work.
    	
        % % If you want to compute the drag forces as functions of time,
        % % put them here.
    	% % I think the ion drag should actually use q(k) instead of 
    	% % q(k-1), because q should be evaluated at a spatial location - 
        % % 5/23/2013
        
            %[f_ix(k),f_iy(k)]=ion_drag(q(k),a,Te(k),Ti(k),ni(k),C,lambda_D(k),lambda_i(k),vix,viy,0,0,ch_model,Kn_R0(k),P0(k),P1(k),Pg1(k),Z,B,gx,gy,species);
            %[f_nx(k),f_ny(k),nu_dn]=neutral_drag(a,rho,nneut,0,0,vx(k-1),vy(k-1),Tn,gx,gy,species);
   	
    	% % put the 1/2 velocity update timestep here if you want to calculate
    	% % vx and vy at the same timesteps as the positions. The time axis
    	% % would no longer be staggered. Essentially, you are rotating and 
    	% % accelerating half a timestep
           
    end

 
    % % if you want to see the trajectory evolve in real time, here is the
    % % code below, just uncomment it:
    %figure(1);drawnow;plot(x,y)
    
    % % July 2013: I had originally intended to allow for a variable
    % % timestep even in the main loop that reflects the variable 
    % % gyroperiod; I have abandoned this however because Boris and other 
    % % leapfrog algorithms require a fixed timestep. I have left this 
    % % vestigal code below.
    
    % % update the newton timestep to reflect the possibility that the 
    % % gyrofrequency is varying in time. should I include a varying newton 
    % % timestep? This means that:
    %dtNwt=H*2*pi*md/abs(q(k))/B;
	%toc;
	%disp(k);
	%pause;
    
end
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % figure out whether we have used a cylindrical or non-cylindrical   
% % profile, and order the densities, fields, etc. accordingly.
%cyl_type=strfind(profile_type,'cyl');
%
%if length(cyl_type)==0
	% % put the E,ni,ne profiles in order, where s will become the ordered 
    % % set of x-positions (low to high x values)
	%[s,i]=sort(x); 
	%E_sx = E_xt(i);
	%E_sy=E_yt(i);
	%%V_s = V_time(i);
	%ni_s=ni_time(i);
	%ne_s=ne_time(i);
    %
    %% % optional: produce a quick diagnostic figure
    %h1=figure;clf;
    % % This plot (2/27/2013) shows the density profiles and on top of  
    % % this, shows the motion of the dust grain.
    %%subplot(2,1,1);
    %%fill([0 0.02 0.02 0],[0 0 0.2 0.2],'y')
    %hold on;
    %plot(x,y,'-b','LineWidth',2);
    %set(gca,'fontsize',16);
    %xlabel('x position');
    %ylabel('y position');
    %title('Slab Profile')
    %%  subplot(2,1,2);
    %%  plot(s,ni_s,'-r','LineWidth',1);hold on;plot(s,ne_s,'--g','LineWidth',2);
    %%  set(gca,'fontsize',16);
    %%  axis square;
    %%  xlabel('x position (m)')
    %%  ylabel('density (m^{-3})')
    %%  legend('Ion Density','Electron Density')
    %
    %% % here is a diagnostic plot to show the electric field and density profiles.
    %% h2=figure;clf
    %% subplot(1,2,2);
    %% %set(gca,'fontsize',16);
    %% [AX,H1,H2]=plotyy(s,E_s,s,V_s,'plot');
    %% %set(gca,'fontsize',16);
    %% set(get(AX(1),'Ylabel'),'String','Electric Field (V/m)') 
    %% set(get(AX(2),'Ylabel'),'String','V_{space} (V)') 
    %% xlabel('x position (m)')
    %% % xlabel('x position (m)')
    %% % ylabel('Electric Field (V/m)')
    %% subplot(1,2,1);
    %
    %% % This plot (4/18/2013) is for showing the x and y componenets of the 
    %% % ion drag force.
    %%h3=figure;
    %%plot(t,f_x,'-c')
    %%hold on;
    %%plot(t,f_y,':m')
    %%set(gca,'fontsize',16);
    %%xlabel('time');
    %%ylabel('Force (N)');
    %%legend('Ion Drag Force, x-direction','Ion Drag Force, y-direction')
%
%else
	%%% if cyl_type gives us any value other than an empty array, then we have
	%%% a cylindrical profile of some type. organize things according.
	%r=sqrt(x.^2+y.^2);
	%[s,i]=sort(r);
	%E_sx = E_xt(i);
	%E_sy=E_yt(i);
	%E_sr=sqrt(E_sx.^2+E_sy.^2);
	%%V_s = V_time(i);
	%ni_s=ni_time(i);
	%ne_s=ne_time(i);
	%s=[-s,s];
	%E_sr=[-E_sr,E_sr];
	%ni_s=[ni_s,ni_s(end:-1:1)];
	%ne_s=[ne_s,ne_s(end:-1:1)];
    %
    %% % optional: produce a quick diagnostic figure
    %figure(1)
    %hold on;
    %plot(x,y,'-b','LineWidth',2);
    %set(gca,'fontsize',16);
    %xlabel('x position');
    %ylabel('y position');
    %title('Cylindrical Profile')
    %R=0.225;
    %xleft=-1.1*R;
    %xright=1.1*R;
    %ylower=-1.1*R;
    %yupper=1.1*R;
%
    %phi_circ=linspace(0,2*pi,1e4);
    %xcirc=R*cos(phi_circ);
    %ycirc=R*sin(phi_circ);
    %plot(xcirc,ycirc,'--r', 'LineWidth',2);
    %axis equal;
%end
% % save some memory by trashing old time series:
%clear E_time;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%
% % KE needs to be reworked so that it's calculated at the full step;
% % currently this is the KE at the half step.
%KE=md*(vx.^2+vy.^2)/2;
% % PE of the dust grain in the field; plug into PE the desired spatial  
% % dependence of the electric field
%PE = q.*E_time;

% % RLd also needs to be reworked for the same reasons as KE
% zero_B = find(B_t==0)
% zero_q = find(q==0)
% finite_q_or_B = find(B_t ~=0)
% if length(zero_q)~=0 || length(zero_B)~=0
%   RLd(zero_q)=inf;
%   RLd(zero_B)=inf;
% 

% Maybe get rid of this definition of larmor radius, because it's wrong.
% if B==0
% 	RLd=inf;	% Infinite gyro-radius if there is no magnetic field.
% else
%     % make sure you take the absolute value of the gyro-frequency!
%     RLd=md.*sqrt(vx.^2+vy.^2)./abs(q)./abs(B_t);
% end
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %
% % Use the following command to save the output variables and the 
% % parameters. Save all variables; almost all quantities here should be 
% % saved. It is faster and easier to clear the unnecessary ones.
clear Q0;clear i;clear k;clear nsteps;
save(strcat(filename,'.mat'));
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end
