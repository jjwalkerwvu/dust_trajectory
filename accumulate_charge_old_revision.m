%% Last modified by Jeff Walker, August 23, 2012, modified Feb 21 2013
%% updated to have density etc. calculations done here instead of in dqdt_models
%% to save computation time.

% % Last modified August 28 2013: I made a better algorithm than this, but I have since
% % reverted to this older code for the JPP paper.
% % See: accumulate_charge_new_revision.m for the "better" code, which still needs some
% % tuning.


function [q,Itot,Kn_R0,P0,P1,Pg1,t_acc]=accumulate_charge(qflag,ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,q,tmax,alph_m,lambda_D,lambda_i,w,t_acc,species); 
%[q]=accumulate_charge(a,alph,alph_m,species,mi,Ti,vthi,Rli,Te,vthe,Rle,ni,ne,B,Z,c_log,C,q,tmax,er,lambda_D,lambda_i,ch_model,V_space);
	
% % explanation of inputs:
% %	qflag 	= whether or not to evaluate equilibrium charge for a given model;
% %		    qflag=1 means do calculate q_eq, qflag=0 means do not calculate q_eq.

% % I've decided to get rid of global vars; they are commented if you feel like using them again.
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
vthi=sqrt(mr/Tau)*vthe;		% local ion (proton) thermal speed, m/s


%% the block of code below used to be in dqdt_models.m, but I've placed it
%% here in the hope that this will reduce computations by only having to
%% calculate it once per coarse timestep, or the value of dt found in
%% dust_gyro_models.m

%nu_i=(4.8e-8)*((ni/1e6)*c_log*Z^4)/(sqrt(species)*Ti^(3/2)); % this expression for the ion collision frequency is from David Pace's website. Have to divide ni by 1e6 to get it in the right units(?)


Ze=1;	% number of electrons collected per adaptive timestep
%
cnt=1;
tchg=0;
while tchg<tmax	% iterate to find equilibrium charge; tmax is the newton timestep
	[Itot(cnt),q0,Kn_R0,P0,P1,Pg1]=charging_models(qflag,ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,q,lambda_D,lambda_i,w,species);
	%% q0 is an output of charging models; maybe it is not needed if using accumulate_charge function
	
	dt=Ze*qe/abs(Itot(cnt));	%% if there is no current to the dust grain, this will be dividing by zero!
   	%% the next line is the charge delay
	dq=alph_m*dt*Itot(cnt);
	%% if not using an alpha:
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
	t(cnt)=tchg;
	qarr(cnt)=q;
	%
	if cnt>2&&qarr(cnt)==qarr(cnt-2)
		break
	end
	cnt=cnt+1;
	%disp(cnt);

end

% just put some stuff down here.
t_acc=0;
I_tot=0;
%cnt

% close all;
% drawnow;
% subplot(2,1,1);plot(t,dqdt/qe,'.');grid on
% subplot(2,1,2);plot(t,qarr/qe,'.');grid on;pause

