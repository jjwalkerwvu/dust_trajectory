% % Last modified by Jeff Walker, August 23, 2012, modified Feb 21 2013
% % updated to have density etc. calculations done here instead of in 
% % dqdt_models to save computation time.
function [q,Itot,Kn_R0,P0,P1,Pg1,t_acc]=accumulate_charge(qflag,...
    ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,q,tmax,alph_m,lambda_D,lambda_i,...
    w,t_acc,species) 

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
% %             the grain. Generally, the Newton timestep will be used here,
% %             but it depends.
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

Ze=1;	% number of electrons collected per adaptive timestep
cnt=1;  % initialize loop counter
tchg=0; % initialize the adaptive timestep; for when dt<tmax

[Itot(cnt),q0,Kn_R0,P0,P1,Pg1]=charging_models(qflag,ch_model,a,alph,...
    Te,Ti,ne,ni,B,Z,C,q,lambda_D,lambda_i,w,species);
% % This is the line of code that I have used in the past to make plots for
% % posters, and the JPP paper.
dt=Ze*qe/abs(Itot(cnt));

% % Maybe the line below is better?? Charge is added/subtracted 1 electron
% % at a time.
dt=Ze*qe/alph_m/abs(Itot(cnt));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
        [Itot(cnt),q0,Kn_R0,P0,P1,Pg1]=charging_models(qflag,ch_model,a,...
            alph,Te,Ti,ne,ni,B,Z,C,q,lambda_D,lambda_i,w,species);
	
        
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
            %t_acc=t_rem;       % the grain reached its final charge value  
                                % during the newton timestep at t_rem, so 
                                % it's been at this value for time = t_acc.
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
    %t_rem=tchg-dt;     % dt has not changed since tchg was calcualted; subtract 
                        % off dt from tchg to find out what tchg was before we
                        % exited the while loop.
    %t_acc=tmax-t_rem;   % We've been at this charge for time = t_acc
end
%cnt
% % figure out what is the change in charge since the last time step.
delta_q=(q-qi)/qe; 
% % wait, what is that line above used for???
%pause
Itot=Itot(end);

% close all;
% drawnow;
% subplot(2,1,1);plot(t,dqdt/qe,'.');grid on
% subplot(2,1,2);plot(t,qarr/qe,'.');grid on;pause

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% stuff I've tacked on for instantaneous charging:
% FIND EQUILIBRIUM SURFACE POTENTIAL.
% figure out dimensionless plasma parameters

% e_mag=a./(sqrt(pi/4)*me*sqrt(2*qe*Te/me)/qe/B);
% Kna=lambda_i/a;
% Tr=Te/Ti;
% Me=w(1)/vthe;
% % could have problems with this later down the road, just sayin. If you use
% Mi=w(2)/vthi;
% % initialize charging loop parameters
% cnt=1;
% tchg=0;
% Z=0;
% % if ions are flowing, we cant use the usual tn_fact=(1+Tr/eta), so prepare
% % for this! The thought here is that the ion temperature 
% if Mi==0
%     tn_fact=(1+Tr/eta);
%     KnD=lambda_D/a;
%     NDe=4/3*pi*ne*lambda_D.^3;
% else
%     % if the string input is mono-energetic, use the definition below for
%     % tn_fact.
%     if strcmp(ch_model,'oml_monoenergetic_ions')==1 || ...
%             strcmp(ch_model,'kortshagen_monoenergetic_ions')==1 || ...
%             strcmp(ch_model,'hutchinson_monoenergetic_ions')==1
%         % This expression uses Mi in terms of the bohm speed, which is 
%         % correct for a mono-energetic population of ions. If it is 
%         % flow-shifted, reconfigure for Mi in terms of the bohm speed!
%         Mi=w(2)/sqrt(e*Te/mi);
%         tn_fact=(1+1/Mi.^2/eta);
%     else
%         % if it's not mono-energetic, use this definition!
%         tn_fact=(1+1/Mi.^2*Tr/2/eta);
%     end
%     KnD=sqrt(eps0*(Te*Te)/(ni*Te+ne*Te)/qe)/a;
%     NDe=4/3*pi*ne.*sqrt(eps0*(Te*Te)/(ni*Te+ne*Te)/qe).^3;
% end
% M=[Me,Mi];
% % currently, no UV is assumed for region 1 (x<0).
% while   Z<=0        % Note that this is a bogus statement; the point is to 
%                     % run the loop until Z repeats itself, that's when we 
%                     % have reached equilibrium surface potential.          
%    	dZdt=dimensionless_charger(ch_model,Z,Tr,mr,M,eta,Kna,KnD,0,e_mag);
%     dt=1/3*KnD/(1+1/KnD)/(tn_fact)/NDe/abs(dZdt);
%    	dZ=1/3*KnD/(1+1/KnD)/(tn_fact)/NDe*sign(dZdt);
%   
%    	Zarr(cnt)=Z;
%    	tarr(cnt)=tchg;
%    	tchg=tchg+dt;
%   	Z=Z+dZ; 
%   	if cnt>2&&Zarr(cnt)==Zarr(cnt-2)
%      	% break out of the loop when the charge begins oscillating back
%        	% and forth between two values.
%       	break
%   	end
%   	cnt=cnt+1;
%         
% %   	drawnow;
% %   	figure(1);plot(tarr,Zarr)
%     
% end
% clear Zarr;clear tarr;
% q=C*Z*Te;
% % itot:
% Itot=C*(sqrt(ne*qe.^2/me/eps0)/2/pi)*dZdt;
% % output garbage to tacc for now:
% t_acc=0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
