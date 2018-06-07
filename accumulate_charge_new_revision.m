%% Last modified by Jeff Walker, August 23, 2012, modified Feb 21 2013
%% updated to have density etc. calculations done here instead of in dqdt_models
%% to save computation time.
function [q,Kn_R0,P0,P1,Pg1,t_acc]=accumulate_charge(qflag,ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,q,tmax,alph_m,lambda_D,lambda_i,w,t_acc,species); 



% % save the initial charge in qi.
qi=q;

% % explanation of inputs:
% %	qflag 	= whether or not to evaluate equilibrium charge for a given 
% %             model; qflag=1 means do calculate q_eq, qflag=0 means do 
% %             not calculate q_eq.
% % ch_model = which charge model you are using, e.g., 'oml', 'kortshagen',
% %             etc.
% % a       = grain size in meters
% % alph    = coefficient of UV illumination (I like to use 0.25)
% % Te      = electron temperature in eV
% % Ti      = ion temperature in eV
% % ne      = local electron density in m^(-3)
% % ni      = local ion density in m^(-3)
% % B       = strength of the magnetic field in Tesla
% % Z       = ionization state of the plasma; I usually just use Z=1.
% % C       = capacitance of the dust grain (4*pi*eps0*a)
% % q       = the current charge on the dust grain. 
% % tmax    = A comparison for how long accumulate_charge should charge up
% %             the grain. Generally, the Newton timestep will be used here.
% % alph_m  = the charge delay parameter
% % lambda_D= the local Debye length, in meters.
% % lambda_i= the local mean free path for ion-neutral charge exchange
% %             collisions.

% % I've decided to get rid of global vars; they are commented if you feel 
% % like using them again.
%global qe;
%global me;
%global mp;
%global eps0;
qe=1.6e-19;
me=9.1e-31;
mp=1.67e-27;
eps0=8.854e-12;
mi=species*mp;


eta=ne/ni;
Tau=Te/Ti;
mr=me/mi;
vthe=sqrt(2*qe*Te/me);	% local electron thermal speed, m/s
vthi=sqrt(mr/Tau)*vthe;	% local ion (proton) thermal speed, m/s

% this expression below is for the ion collision frequency is from David  
% Pace's website. Have to divide ni by 1e6 to get it in the right units(?)
% Additionally, this is applicable for ion-ion collisions, not ion-neutral
% collisions.
%nu_i=(4.8e-8)*((ni/1e6)*c_log*Z^4)/(sqrt(species)*Ti^(3/2)); 


Ze=1;	% number of electrons collected per adaptive timestep
cnt=1;  % initialize loop counter
tchg=0; % initialize the adaptive timestep; for when dt<tmax

[Itot(cnt),q0,Kn_R0,P0,P1,Pg1]=charging_models(qflag,ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,q,lambda_D,lambda_i,w,species);
% % This is the line of code that I have used in the past to make plots for
% % posters, and the JPP paper.
dt=Ze*qe/abs(Itot(cnt));

% % Maybe the line below is better?? Charge is added/subtracted 1 electron
% % at a time.
dt=Ze*qe/alph_m/abs(Itot(cnt));

% % Alright, here's how it works: If the time for the grain to charge by an
% % increment/decrement of 1 electron is greater than the newton timestep,
% % update the charge by a fractional amount!
if dt>tmax
    %dq=alph_m*tmax*Itot(cnt); 
    %q=q+dq;
    % if dtchg>dtnwt, then start timing how long it's been since the last
    % charge update
    t_acc=t_acc+tmax;
    % % bunchabs: ensure that the adaptive timestep works even if
    % dtchg>dtNwt!!!
    %t_rem=dt-tmax;
    if t_acc>dt
        %dq=alph_m*tmax*Itot(cnt); 
        % written in such a way so that the grain is updated by +/- 1
        % electron!
        dq=alph_m*dt*Itot(cnt);
        %dq=dt*Itot(cnt);
        % q(t) gets "lumped" into a time that is some integer multiple of
        % dtNwt, but the time since the last charge update is maintained in
        % the t_acc variable.
        q=q+dq; 
        % find out how far into THIS newton timestep we went before the
        % particle was given +/- 1 electron
        %[Itot(cnt),q0,Kn_R0,P0,P1,Pg1]=charging_models(qflag,ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,q,lambda_D,lambda_i,w);
        %dt_update=Ze*qe/abs(Itot(cnt));
        % % figure out what the UPDATED current should have been when the
        % % charge should have been updated. this means we will get the
        % % correct value of the current at the charge update.
        t_rem=rem(t_acc,dt);
        % reset the time since the last update; which is 
        t_acc=tmax-t_rem;
        
    end
   
% % If the time for the grain to charge up or down by one electron is less
% % than the newton timestep, use the adaptive timestep to charge the
% % grain.
else

    while tchg<tmax     % iterate to find equilibrium charge; 
                        % tmax is the newton timestep
        %disp('tchg<dtNwt')
        % % q0 is an output of charging models; maybe it is not needed if 
        % % using accumulate_charge function
        [Itot(cnt),q0,Kn_R0,P0,P1,Pg1]=charging_models(qflag,ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,q,lambda_D,lambda_i,w,species);
	
        
        % % if there is no current to the dust grain, this will be dividing 
        % % by zero!
        dt=Ze*qe/abs(Itot(cnt));	
    % % Maybe the line below is better??
        %dt=Ze*qe/alph_m/abs(Itot(cnt));
    % % the next line is the charge delay
        dq=alph_m*dt*Itot(cnt);
	% % if not using an alpha:
	%dq=dt*Itot(cnt);
%	if abs(dq/qe)<1
%		if dt>tmax
%			'!'
%		else
%			dt=qe/abs(alph_m*dqdt(cnt));
%			dq=qe;
%		end
%	end
        q=q+dq;
        tchg=tchg+dt;
        % % How do I fix things to include t_acc? August 2013
        %tchg=t_acc+tchg+dt;
        
        qarr(cnt)=q;
        % optional: make a time array
        t(cnt)=tchg;
	
        if cnt>2&&qarr(cnt)==qarr(cnt-2)
            % % because we have broken out of the loop early, we need to
            % % find out how much time was left in the newton timestep.
            %t_rem=tmax-tchg;
            %t_acc=t_rem;    % the grain reached its final charge value during 
                            % the newton timestep at t_rem, so it's been at
                            % this value for time = t_acc.
            break
        end
        cnt=cnt+1;
	
    % various diagnostic commands, not usually necessary
    %disp(cnt)
    %drawnow;
    %subplot(2,1,2);plot(t,qarr/qe,'.');grid on;
    end
    % % after exiting the while loop, tchg>tmax. Find out what was tchg one
    % % charging increment before this.
    %t_rem=tchg-dt;  % dt has not changed since tchg was calcualted; subtract 
                    % off dt from tchg to find out what tchg was before we
                    % exited the while loop.
    %t_acc=tmax-t_rem;   % We've been at this charge for time = t_acc
end
%cnt
% % figure out what is the change in charge since the last time step.
delta_q=(q-qi)/qe; 
%pause

% close all;
% drawnow;
% subplot(2,1,1);plot(t,dqdt/qe,'.');grid on
% subplot(2,1,2);plot(t,qarr/qe,'.');grid on;pause

