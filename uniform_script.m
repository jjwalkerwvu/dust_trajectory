% example_dust_trajectory_script.m

% This is just a script function that runs dust_trajectory.m and 
% gyrophaser.m. Just replace all of the variables or strings with your
% desired inputs, press f5 in matlab and this script will run the
% simulation, then it will run gyrophaser.m, which appends the target 
% filename with the gyro-averaged quantities. After this file is created,
% you can analyze the data in gyro-phase or in time.

% Here are the variables:

% dust grain radius in meters
a=8e-7;	% to make sure a~1.6 microns in diameter.
% density of the dust grain in kg/m^3
rho=1e3;	% density of water
% initial position of the grain, with the units given in meters
r_initial=[0 0];
% initial velocity of the grain, with the units given in meters/second
v_initial=[0.003 0];
% mass number of ions in the plasma
species=40; % argon mass number = 40
% choose a charge model, can be 'oml', 'hutchinson', 'kortshagen', or 
% 'phgk'.
ch_model='oml';
% choose a profile type; look at profiles.m for examples. must be entered
% as a string.
profile_type='uniform';
% baseline density, in units of ions/electrons m^{-3}
n0=1e14;
% baseline electron temperature in units of eV
Te0=1.6;
% baseline ion temperature in units of eV
Ti0=1/40;
% ionization of plasma
Z=1;    % Z=1 means singly ionized ions
% neutral gas pressure in mTorr
P_mtorr=0;
% convert pressure in mTorr to Pa
P=P_mtorr/7.5;	
% ADJUSTABLE CHARGING RATE PARAMETER: alphm=1 means there is no
% restriction; alphm<1 implies that the grain charges more slowly than
% normal.
alphm=1;
% approximate number of gryocycles
cycles=10;
% approximate number of points per gyrocycle
points=2e3; % 2000 points/cycle is a reasonable number to work with; 
            % a smaller number will make simulations run quicker but with
            % less accuracy, while a larger number will make simulations
            % run slower but with greater accuracy.
            
% pick a filename for your output; do not need to include .mat extension 
% because the code does this for you.
filename='test';

% choose a method of advancing the trajectory of the dust grain; your
% options are 'boris_pusher', 'iterative_pusher',
% 'corotating_boris_pusher', 'sheath_boris_pusher', 
% still need to add: 'corotating_iterative_pusher'! ~ April 9 2014
particle_pusher='boris_pusher';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% alternative: concatenate all of the relevant information into the data
% file. comment these lines out if you don't want to use them.

% append the dust size information to the file, in units of um.
a_string = num2str(a/1e-6);
a_string(find(a_string=='.'))='_';
% append the density information to the file (in units of m^-3).
n_string = num2str(n0);
indices = find(n_string=='0');
n_string = strcat('_n_',n_string(1),n_string(2:indices(1)-1),...
    'e',num2str(length(n_string)-1));

p_string=num2str(P_mtorr); 
p_string(find(p_string=='.'))='_';
p_string=strcat('_P_',p_string,'mTorr');

filename=strcat(ch_model,'_',profile_type,'_',particle_pusher,'_a',...
    a_string,'um',n_string,p_string);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Now that all the inputs are given, run the simulation.
% The '...' ellipses are used to make the function call more readable.

% Please Note that you can call dust_trajectory.m without specifying the
% outputs. E.g., you use:
% dust_trajectory(a,rho,r_initial,v_initial,species,ch_model,...
%   profile_type,n0,Te0,Ti0,Z,P,alphm,cycles,points,filename,...
%   particle_pusher);
%
% Instead of:
% [t,q,x,y,vx,vy,RLd,ne,ni,V_time,E_xt,E_yt,B_t,f_ix,f_iy,f_nx,f_ny,...
%     lambda_D,lambda_i,Kn_R0,P0,P1,Pg1,Vgrain]=...
%     dust_trajectory(a,rho,r_initial,v_initial,species,ch_model,...
%     profile_type,n0,Te0,Ti0,Z,P,alphm,cycles,points,filename,...
%     particle_pusher);
%
% The output variables will be saved to the target filename, so you don't
% need to keep the output variables in memory. Thus, you can use the line
% below:
dust_trajectory(a,rho,r_initial,v_initial,species,ch_model,...
    profile_type,n0,Te0,Ti0,Z,P,alphm,cycles,points,filename,...
    particle_pusher);

load(strcat(filename,'.mat'));
plot(x*1e3,y*1e3);
set(gca,'fontsize',16)
xlabel('x (mm)')
ylabel('y (mm)')
% make the graph background transparent
set(gcf, 'Color', [1,1,1]);
xlim([-1.1*max(x)*1e3 1.1*max(x)*1e3]);
axis square;
%export_fig('uniform_gyromotion.png')

% Now run gyro_phaser.m. Pick a starting angle in degrees. 
% If vx>0 and vy=0, then use 270 as a starting angle.
% If vy>0 and vx=, then use 0 as a starting angle, etc.
%phi_start=270;
%gyrophaser(filename,phi_start);