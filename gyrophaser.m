% % This program produces gyrophase-style plots, instead of "configuration 
% % space" plots.
% % This corresponds to Jeff Idea #1, see my notebook #5, page 39
% % 
% % Please note: this code only runs well for re-tracing cycloidal orbits.
% % This means that if you had a high electric field, you need to have a
% % much higher initial velocity or else you will get garbage results.

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % this function should be called from a script where the global variables 
% % listed below are defined within the script

function [phase,q_polar,I_polar,phase_polar,xgc,ygc,rgc,phi_gc,...
    vxgc,vygc,vrgc,vphi_gc,xc,yc,tgc]=gyrophaser(data_file,phi_start); 

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Brief Explanation of inputs and how this program works
% %	Gyrophaser.m will find the guiding center drifts of a dust grain. 
% %	By default, it will use slab geometry to compute vxgc and vygc, but if 
% %	a cylindrical profile has been chosen, then rgc and phi_gc will be 
% %	calculated. In the slab case, vygc generally corresponds to a 
% %	diamagnetic-like drift, while vxgc corresponds to a 
% %	gyrophase-like drift, whereas in the cylindrical case vrgc corresponds 
% % to a gyrophase-like drift and vphi_gc corresponds to a diamagnetic-like 
% % drift.

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

% % ERROR CHECKING: IF phi_start is not given, use a default of 
% % phi_start=90. Apparently this does not actually work, figure out how to 
% % do the error checking
if isempty(phi_start);
    phi_start=360;
    disp('phi_start=360 degrees');
end
% Older code below VVVVVVVVVVVVVVVVVVVVVVVVVVV
% % set the starting angle for the gyroavering analysis in degrees. This 
% % can be phase(1), or the initial gyro-angle, but it does not have to be.
% %  However, phi_start MUST BE POSITIVE, AND IT MUST BE LESS THAN 2*PI. 
%phi_start=90;
%phi_start=270;  % Chosen for the JPP paper!!!
% Older code above ^^^^^^^^^^^^^^^^^^^^^^^^^^^
% % June 3 2013 NOTE:
% % gyrophaser does not give the correct drift magnitudes for abrupt 
% % inhomogeneities; this will require further refinement.
% Older code above ^^^^^^^^^^^^^^^^^^^^^^^^^^^

% % use this function to put things into a gyrophase plot, and also get 
% % drifts. MAKE SURE YOU'VE LOADED THE DATA FILE!
load(strcat(data_file,'.mat'));
% This line below is garbage!
%RLd=md.*sqrt(vx.^2+vy.^2)./abs(q)./abs(B_t);
v_perp=sqrt(vx.^2+vy.^2);

% % prepare to scroll through all of the data; these should all be arrays 
% % the size of t.
nsteps=length(t);
phase=zeros(size(t));
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %	Options for computing gyrocenter:
% %	

% % average_method: get xgc,ygc by averaging x(t), y(t) over a gyrocycle to 
% % get the gyrocenters
average_method=0;

% % larmor_average: use the time-average Larmor radius in x,y configuration 
% % space to get xgc,ygc.
if average_method==1
	larmor_average=0;
else
	larmor_average=1;
end

% % Initialize some counters and other variables:
% % initialize a general counter for the loop. This counter is separate  
% % from the cycle_counter, and is just used to update the phase angle 
% % properly.
counter=1;
% % set phi_prev=0
phi_prev=0;
% % initialize the first element of the array to 0; this will get changed 
% % later.
new_cycle(1)=1;


% % start the cycle counter at 1; everytime the phase angle passes 
% % phi_start this counter will increment by one.
cycle_count=1;

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% figure out if we have a cylindrical profile or not.
if isempty(strfind(profile_type,'cyl'))~=1
% if the value of the above statement is true, we have a cyl. profile
    ar=zeros(size(t));
    aphi=zeros(size(t));
    % make arrays for radial and azimuthal grain velocities:
    r=sqrt(x.^2+y.^2);
    vr=(x.*vx+y.*vy)./r;
    vphi=(x.*vy-y.*vx)./r;
    
    % find the phi-coordinate in cylindrical geometry; this is NOT the
    % gyro-phase angle, it's the angle in configuration space!
    phi_config(1)=improved_arctan(x(1),y(1));
    phi_last=phi_config(1);
    for i=2:length(x);
        phi_temp=improved_arctan(x(i),y(i));
        dphi=phi_temp-phi_last;
        if abs(dphi)>pi/2;
            % if this statement is true, then we have switched from 0 radians
            % to ~<2*pi radians.
            dphi=phi_temp-2*pi;
        end
        phi_config(i)=phi_config(i-1)+dphi;
        phi_last=phi_temp;
    
    %     % this is for the arc length:
    %     s(i-1)=(v_perp(i)+v_perp(i-1))*(t(i)-t(i-1));
        if i<length(x)         
            ar(i)=(vr(i+1)-vr(i-1))/2/(t(2)-t(1));
            aphi(i)=(vphi(i+1)-vphi(i-1))/2/(t(2)-t(1));
        end
    end;
    clear phi_last;clear phi_temp;
    
    % these lines are needed, kinda bs but tacked on anyway to keep
    % acceleration arrays the same size as velocity arrays.
    ar(1)=ar(2);
    ar(end)=ar(end-1);
    aphi(1)=aphi(2);
    aphi(end)=aphi(end-1);
    
    % now calculate larmor radius.
    RLd=v_perp.^3./(vr.*aphi-vphi.*ar);
    
    for k=1:nsteps
        % The point of all these logical statements is to ensure that we are 
        % correctly solving for the gyro-phase angle in configuration space.
	
        % In quadrant 1.) of the standard circle, but exactly at 0 degrees in 
        % cofig. space
        if vphi(k)>=0 && vr(k)==0
            phi=0;
        end
        % In quadrant 2.) of the standard circle, but exactly at 90 degrees in 
        % cofig. space
        if vphi(k)==0 && vr(k)<=0
            phi=pi/2;
        end
        % In quadrant 3.) of the standard circle, but exactly at 180 degrees in 
        % cofig. space
        if vphi(k)<=0 && vr(k)==0
            phi=pi;
        end
        % In quadrant 4.) of the standard circle, but exactly at 270 degrees in 
        % cofig. space
        if vphi(k)==0 && vr(k)>=0
            phi=1.5*pi;
        end

        % In quadrant 1.) of the standard circle.
        if vphi(k)>0 && vr(k)<0	
            phi=atan(abs(vr(k)/vphi(k)));
            % update the "revolutions" counter, which is good for the first 2
            % iterations of the loop.
            if (k==1 || k==2) && phi_prev>=pi/2
                counter=counter+1; 
                %new_cycle(counter-1)=k;
            end
            % update revolutions counter for ccw rotation
            if k>=3 && strcmp(direction,'ccw') && phi_prev>=pi/2
			counter=counter+1; 
			%new_cycle(counter-1)=k;
            end      
        end
        % In quadrant 2.) of the standard circle.
        if vphi(k)<0 && vr(k)<0
            phi=atan(abs(vphi(k)/vr(k)))+pi/2;
			
        end
        % In quadrant 3.) of the standard circle.
        if vphi(k)<0 && vr(k)>0
            phi=atan(abs(vr(k)/vphi(k)))+pi;	 		
        end
        % In quadrant 4.) of the standard circle.
        if vphi(k)>0 && vr(k)>0
            phi=atan(abs(vphi(k)/vr(k)))+1.5*pi;
            % have to add this line to deal with 'clockwise' rotation
            if k>=3 && (strcmp(direction,'cw') && phi_prev<=pi/2)
                counter=counter+1;
            end
        end    
        % to take care of the first element of the phase array in this loop.
        if k==1
            phase(1)=phi;
        end
    
        % at k=2, if the current angle phi is greater than the previous
        % calculated angle phi_prev, OR if the current angle phi 
        if k==2 && (phi>phi_prev || (phi<pi/2 && phi_prev> 3*pi/2))
            direction ='ccw';
        end
        % the cw case:
        if k==2 && (phi<phi_prev || (phi>3*pi/2 && phi_prev<pi/2))
            direction ='cw';  
        end
    
        % after k=2 case has arisen, we can now specify how to update the phase
        % angle.
        % the ccw case:
        if k>=2 && strcmp(direction,'ccw')
            phase(k)=phi+2*pi*(counter-1);
        end
        % the cw case:
        if k>=2 && strcmp(direction,'cw') 
            if (k==2 && phi>3*pi/2 && phi_prev<pi/2)
                phase(1)=phi_prev+2*pi;
            end
            phase(k)=phi-2*pi*(counter-1);
        end
    
    
    % August 27, 2013: 
    % This line below has become quite complicated compared to how it 
    % started out. You need to make sure that the computed angle 
    % 0<phi<2*pi is greater than the designated starting angle in order to  
    % update which gyro-cycle you are on. Also, either the previously 
    % calculated angle phi is less than the starting angle, or phi is less 
    % than the previously calculated value of phi.
    % JANUARY 23, 2014: 
    % The statement must be ammended for grains that gyrate in a clockwise
    % direction!!
     
    % The "counter-clockwise" gyrating statement
        if k>=2 && phi>phi_start*pi/180 && ...
            (phi_prev<=phi_start*pi/180 || phi<phi_prev) && ...
            strcmp(direction,'ccw')

            new_cycle(cycle_count)=k; 
            cycle_count=cycle_count+1;
        end
        % The "clockwise" gyrating statement
        if k>=2 && (phi<phi_start*pi/180 && ...
            (phi_prev>=phi_start*pi/180 || phi>phi_prev)) && ...
            strcmp(direction,'cw')
        
            new_cycle(cycle_count)=k; 
            cycle_count=cycle_count+1;  
        end  
            
        % Need to set phi_prev equal to the current value of phi.
        phi_prev=phi;    
    end

% initialize the arrays which will become the instantaneous guiding center
% position; these arrays are filled up with xc, yc arrays during each cycle 
% iteration

    rc=[];
    phi_c=[];
    rc_temp=[];
    phi_c_temp=[];

    % the part you need to subtract off:
    %rs=sqrt(4*r.^2+(sin(phi_config)).^2-RLd.^2.*(sin(phase)).^2);

    xc=[];
    yc=[];
    xc_temp=[];
    yc_temp=[];
    for i=2:cycle_count-1
    % Need to get the correct indices for each gyrocycle, which will be
    % used in the lines below. This algorithm currently forces the dust 
    % grain to start at zero degrees.
  
        cycle_indices=new_cycle(i-1):new_cycle(i);
        tgc(i)=t(new_cycle(i))-t(new_cycle(i-1));
        
    	xc=x(cycle_indices)-...
         	abs(RLd(cycle_indices)).*cos(phase(cycle_indices));
       
     	%xgc(i)=trapz(t(cycle_indices),xc)/tgc(i);
    	
    	% % Do the same thing as above, but for y??
      	yc=y(cycle_indices)-...
         	abs(RLd(cycle_indices)).*sin(phase(cycle_indices));

     	%ygc(i)=trapz(t(cycle_indices),yc)/tgc(i);
    
    % % HOW DO I FIX THIS FOR CYL. GEOMETRY???
%         rc=r(cycle_indices)-...
%             RLd(cycle_indices).*(sin(phase(cycle_indices)).*phi_config(cycle_indices)+...
%             cos(phase(cycle_indices)).*phi_config(cycle_indices));
%         
%         phi_c=0;
%         
%         rgc(i)=trapz(t(cycle_indices),rc)/tgc(i);
    	
    	% % Do the same thing as above, but for phi; need to use phi_config
    	% % here; 
%         phi_c=r(cycle_indices).*phi_config(cycle_indices)-...
%            abs(RLd(cycle_indices)).*sin(phase(cycle_indices));
%         phi_c=phi_config(cycle_indices)-...
%         	abs(RLd(cycle_indices)).*sin(phase(cycle_indices));
%     phi_c=r(cycle_indices).*phi_config(cycle_indices)-...
%         abs(RLd(cycle_indices)).*sin(phase(cycle_indices));
    
%     phi_c=r(cycle_indices).*phi_config(cycle_indices)-...
%         RLd(cycle_indices).*sin(cycle_indices);
%   	phi_gc(i)=trapz(t(cycle_indices),phi_c)/tgc(i);
       
%         for k=1:length(xc)
%         	
%             [V_space,Ex,Ey,B,vix,viy,vex,vey,vnx,vny,gx,gy,n_i,n_e,alph,...
%                 T_i,T_e,nneut,l_i,l_D,corot_period]=...
%                 profiles(Ti0,Te0,n0,t(k),xc(k),yc(k),...
%                 profile_type,P,species);
%             % capacitance of grain with adjustment for debye shielding
%             Cap=4*pi*eps0*a*(1+a/l_D);
%             w_e=sqrt((vx(k)-vex).^2+(vy(k)-vey).^2);
%             w_i=sqrt((vx(k)-vix).^2+(vy(k)-viy).^2);
%             w_vec=[w_e w_i];
%             [Itot_eq,qeq,KnR_eq,P_0,P_1,P_g1]=...
%                 charging_models(1,ch_model,a,alph,T_e,T_i,n_e,n_i,B,Z,Cap,0,...
%                 l_D,l_i,w_vec,species);
%             q_polar(cycle_indices(k))=q(cycle_indices(k))/qeq;
%             I_polar(cycle_indices(k))=Itot(cycle_indices(k))/Itot_eq;
%             phase_polar(cycle_indices(k))=phase(cycle_indices(k));
%             % % while we are in this loop, determine what is phi_c, or the
%             % % angle of the instantaneous guiding center for cylindrical
%             % % coordinates.
%             phi_c(k) = improved_arctan(xc(k),yc(k));
%             % % OPTIONAL LINE BELOW:
%             % % COMPUTE THE GYRO-AVERAGED VALUE OF dq(x)/dx!
%             % % --how to do this??
%         end
       
        %rc_temp=[rc_temp,rc];
        %phi_c_temp=[phi_c_temp,phi_c];
        xc_temp=[xc_temp,xc];
        yc_temp=[yc_temp,yc];
        
    
        
       
        % ideally, put some code here in case neither of the above cases are 
        % selected.
   
        % % for cyl. coordinates:
        %vrgc(i-1)=(rgc(i)-rgc(i-1))/(tgc(i));
        vrgc(i-1)=trapz(t(cycle_indices),vr(cycle_indices));
        % % the phi component of the guiding center drift is tricky, because 
        % % rgc might be changing so the usual vphi=r*dphi/dt may give us wrong 
        % % answers; just use an average value of r for now.
        %vphi_gc(i-1)=(rgc(i)+rgc(i-1))*(phi_gc(i)-phi_gc(i-1))/(tgc(i))/2;
        vphi_gc(i-1)=trapz(t(cycle_indices),vphi(cycle_indices));
    
	
    end
    
    % arrays I either don't need or have not figured out for the cyl.
    % profile.
 	q_polar=0;
    I_polar=0;
    phase_polar=0;
    %xgc=0;
    %ygc=0;
    %rgc=0;
    phi_gc=0;
    vxgc=0;
    vygc=0;
    vrgc=0;
    vphi_gc=0;
    %xc=0;
    %yc=0;
    tgc=0;
    %rc=rc_temp;
    xc=xc_temp;
    yc=yc_temp;

else
% If the value of the if statement was false, we don't have a cyl. profile,
% but rather a slab profile. proceed as normal.
    ax=zeros(size(t));
    ay=zeros(size(t));
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % 	In the loop below, we must do the following:
% %	get instantaneous gyrocenter by using v_perp method
% January 2014: I think that the code below only works when phi is
% increasing, i.e., the grain gyrates in a counter-clockwise direction.
    for k=1:nsteps
	% The point of all these logical statements is to ensure that we are 
	% correctly solving for the gyro-phase angle in configuration space.
	
	% In quadrant 1.) of the standard circle, but exactly at 0 degrees in 
    % cofig. space
        if vy(k)>=0 && vx(k)==0
            phi=0;
        end
	% In quadrant 2.) of the standard circle, but exactly at 90 degrees in 
    % cofig. space
        if vy(k)==0 && vx(k)<=0
            phi=pi/2;
        end
	% In quadrant 3.) of the standard circle, but exactly at 180 degrees in 
    % cofig. space
        if vy(k)<=0 && vx(k)==0
            phi=pi;
        end
	% In quadrant 4.) of the standard circle, but exactly at 270 degrees in 
    % cofig. space
        if vy(k)==0 && vx(k)>=0
            phi=1.5*pi;
        end

	% In quadrant 1.) of the standard circle.
        if vy(k)>0 && vx(k)<0	
            phi=atan(abs(vx(k)/vy(k)));
        % update the "revolutions" counter, which is good for the first 2
        % iterations of the loop.
            if (k==1 || k==2) && phi_prev>=pi/2
                counter=counter+1; 
                %new_cycle(counter-1)=k;
            end
            % update revolutions counter for ccw rotation
            if k>=3 && strcmp(direction,'ccw') && phi_prev>=pi/2
                counter=counter+1; 
                %new_cycle(counter-1)=k;
            end      
        end
        % In quadrant 2.) of the standard circle.
        if vy(k)<0 && vx(k)<0
            phi=atan(abs(vy(k)/vx(k)))+pi/2;
			
        end
        % In quadrant 3.) of the standard circle.
        if vy(k)<0 && vx(k)>0
            phi=atan(abs(vx(k)/vy(k)))+pi;	 		
        end
        % In quadrant 4.) of the standard circle.
        if vy(k)>0 && vx(k)>0
            phi=atan(abs(vy(k)/vx(k)))+1.5*pi;
            % have to add this line to deal with 'clockwise' rotation
            if k>=3 && (strcmp(direction,'cw') && phi_prev<=pi/2)
                counter=counter+1;
            end
        end

        %phase(k)=phi+2*pi*(counter-1);
    
        % to take care of the first element of the phase array in this loop.
        if k==1
            phase(1)=phi;
        end
    
    % at k=2, if the current angle phi is greater than the previous
    % calculated angle phi_prev, OR if the current angle phi 
        if k==2 && (phi>phi_prev || (phi<pi/2 && phi_prev> 3*pi/2))
            direction ='ccw';
        end
        % the cw case:
        if k==2 && (phi<phi_prev || (phi>3*pi/2 && phi_prev<pi/2))
            direction ='cw';  
        end
    
        % after k=2 case has arisen, we can now specify how to update the 
        % phase angle.
        % the ccw case:
        if k>=2 && strcmp(direction,'ccw')
            phase(k)=phi+2*pi*(counter-1);
        end
        % the cw case:
        if k>=2 && strcmp(direction,'cw') 
            if (k==2 && phi>3*pi/2 && phi_prev<pi/2)
                phase(1)=phi_prev+2*pi;
            end
            phase(k)=phi-2*pi*(counter-1);
        end
    
    
    % August 27, 2013: 
    % This line below has become quite complicated compared to how it 
    % started out. You need to make sure that the computed angle 
    % 0<phi<2*pi is greater than the designated starting angle in order to  
    % update which gyro-cycle you are on. Also, either the previously 
    % calculated angle phi is less than the starting angle, or phi is less 
    % than the previously calculated value of phi.
    % JANUARY 23, 2014: 
    % The statement must be ammended for grains that gyrate in a clockwise
    % direction!!
     
        % The "counter-clockwise" gyrating statement
        if k>=2 && phi>phi_start*pi/180 && ...
            (phi_prev<=phi_start*pi/180 || phi<phi_prev) && ...
            strcmp(direction,'ccw')

            new_cycle(cycle_count)=k; 
            cycle_count=cycle_count+1;
        end
        % The "clockwise" gyrating statement
        if k>=2 && (phi<phi_start*pi/180 && ...
            (phi_prev>=phi_start*pi/180 || phi>phi_prev)) && ...
            strcmp(direction,'cw')
        
            new_cycle(cycle_count)=k; 
            cycle_count=cycle_count+1;  
        end  
            
        % Need to set phi_prev equal to the current value of phi.
        phi_prev=phi;  
        
        % last little part here is needed for calculating RLd.
        if k>=2 && k<=nsteps-1
            ax(k)=(vx(k+1)-vx(k-1))/2/(t(2)-t(1));
            ay(k)=(vy(k+1)-vy(k-1))/2/(t(2)-t(1));
        end
    end
    % for the clockwise direction, I need to subtract off 180 degrees to 
    % get the correct phase.
    if strcmp(direction,'cw')
        phase=phase-pi;
    end
    
    % these lines are needed, kinda bs but tacked on anyway to keep
    % acceleration arrays the same size as velocity arrays.
    ax(1)=ax(2);
    ax(end)=ax(end-1);
    ay(1)=ay(2);
    ay(end)=ay(end-1);
    
    % now calculate larmor radius.
    RLd=v_perp.^3./(vx.*ay-vy.*ax);
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % % make a plot to check on the phase, if desired.
    %figure(1);clf;
    %plot(phase,'.')
    % % Not sure if the following can go into the loop above, but the idea   
    % % here is to calculate the average guiding center over a gyroperiod.
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % % set the initial xgc position to the initial x-position
    % below must be fixed for arbitrary starting gyro-average angle
    % xgc(1)=x(new_cycle(1))-RLd(new_cycle(1))*cos(new_cycle(1));	
    % %^^^^ Is that right???
    % %% set the initial ygc position:
    % ygc(1)=y(new_cycle(1))+RLd(new_cycle(1))*sin(new_cycle(1));
    xgc(1)=x(new_cycle(1));     %%<<< Is this right???
    % % set the initial ygc position:
    ygc(1)=y(new_cycle(1))-RLd(new_cycle(1))*sign(q(new_cycle(1)));

    % % Can now set the initial rgc position:
    rgc(1)=sqrt(xgc(1).^2+ygc(1).^2);
    phi(1)=improved_arctan(xgc(1),ygc(1));

    % % obviously, tgc(1) should be then zero:
    tgc(1)=0;
% % KEEP IN MIND: the guiding centers listed above are the INITIAL guiding 
% % center coordinates BEFORE any motion has happened. The real guiding 
% % centers for the first gyrocycle corresponds to xgc(2) and ygc(2); the 
% % above is merely necessary to set up so that the guiding center 
% % velocities can be calculated, starting even with the first gyrocycle. 
% % You can remove these.

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% initialize the arrays which will become the instantaneous guiding center
% position; these arrays are filled up with xc, yc arrays during each cycle 
% iteration
    xc_temp=[];
    yc_temp=[];
    Bc_temp=[];
    rc_temp=[];
    phi_c_temp=[];

    % initialize arrays to save some computational time:
    tgc=zeros(1,cycle_count-1);
    xgc=zeros(size(tgc));
    ygc=zeros(size(tgc));
    rgc=zeros(size(tgc));
    phi_gc=zeros(size(tgc));
    % find out what xgc(1) and ygc(1) should be! This is at
    % gyro-phase=phi_start.
    xgc(1)=x(new_cycle(1))-...
        abs(RLd(new_cycle(1))).*cos(phase(new_cycle(1)));
    ygc(1)=y((new_cycle(1)))-...
        abs(RLd((new_cycle(1)))).*sin(phase((new_cycle(1))));

    for i=2:cycle_count-1
    % Need to get the correct indices for each gyrocycle, which will be
    % used in the lines below. This algorithm currently forces the dust 
    % grain to start at zero degrees.
  
        cycle_indices=new_cycle(i-1):new_cycle(i);
        tgc(i)=t(new_cycle(i))-t(new_cycle(i-1));
    
        if average_method==1
        % UNFORTUNATELY, I THINK THE AVERAGE METHOD JUST DOES NOT WORK FOR
        % THE PURPOSE OF FINDING THE GYROCENTER. USE LARMOR_AVERAGE 
        % INSTEAD.
        % -JJW June 2013
    	% below is the simpler, more intuitive, averaging over x during a 
        % cycle method:
            xgc(i)=trapz(t(cycle_indices),x(cycle_indices))/tgc(i);
    
    	% % below is the simpler, averaging over y during a cycle method:
            ygc(i)=trapz(t(cycle_indices),y(cycle_indices))/tgc(i);
  
    	% % for cyl. coordinates
            rgc(i)=sqrt(xgc(i).^2+ygc(i).^2);
            phi_gc(i)=improved_arctan(xgc(i),ygc(i));
        end
    
        if larmor_average==1
    	% % first, here's the xgc(t) method, only valid if vx(0)=negative:
            xc=x(cycle_indices)-...
                abs(RLd(cycle_indices)).*cos(phase(cycle_indices));
       
            xgc(i)=trapz(t(cycle_indices),xc)/tgc(i);
    	
    	% % Do the same thing as above, but for y??
            yc=y(cycle_indices)-...
                abs(RLd(cycle_indices)).*sin(phase(cycle_indices));

            ygc(i)=trapz(t(cycle_indices),yc)/tgc(i);
        
        % compute the instantaneous guiding center position in cylindrical 
        % coordinates. phi_c gets computed later in the loop below.
            rc=sqrt(xc.^2+yc.^2);
            rgc(i)=trapz(t(cycle_indices),rc)/tgc(i);
        % initialize the phi_c array to save time:
            phi_c = zeros(1,length(xc));
        
            for k=1:length(xc)
        	
                [V_space,Ex,Ey,B,vix,viy,vex,vey,vnx,vny,gx,gy,n_i,n_e,...
                    alph,T_i,T_e,nneut,l_i,l_D,corot_period]=...
                    profiles(Ti0,Te0,n0,t(k),xc(k),yc(k),...
                    profile_type,P,species);
            % capacitance of grain with adjustment for debye shielding
                Cap=4*pi*eps0*a*(1+a/l_D);
                w_e=sqrt((vx(k)-vex).^2+(vy(k)-vey).^2);
                w_i=sqrt((vx(k)-vix).^2+(vy(k)-viy).^2);
                w_vec=[w_e w_i];
                [Itot_eq,qeq,KnR_eq,P_0,P_1,P_g1]=...
                    charging_models(1,ch_model,a,alph,T_e,T_i,n_e,n_i,...
                    B,Z,Cap,0,l_D,l_i,w_vec,species);
                Bc(cycle_indices(k))=B;
                q_polar(cycle_indices(k))=q(cycle_indices(k))/qeq;
                I_polar(cycle_indices(k))=Itot(cycle_indices(k))/Itot_eq;
                phase_polar(cycle_indices(k))=phase(cycle_indices(k));
            % % while we are in this loop, determine what is phi_c, or the
            % % angle of the instantaneous guiding center for cylindrical
            % % coordinates.
                %phi_c(k) = improved_arctan(xc(k),yc(k));
            % % OPTIONAL LINE BELOW:
            % % COMPUTE THE GYRO-AVERAGED VALUE OF dq(x)/dx!
            % % --how to do this??
            end
            xc_temp=[xc_temp,xc];
            yc_temp=[yc_temp,yc];
            Bc_temp=[Bc_temp,Bc];
            rc_temp=[rc_temp,rc];
            phi_c_temp=[phi_c_temp,phi_c];
        % fix this for cyl. coordinates!!! May 28 2013
    	% for cyl. coordinates
    	%rgc(i)=sqrt(xgc(i).^2+ygc(i).^2);
    	%phi_gc(i)=improved_arctan(xgc(i),ygc(i));
        
        % % now that phi_c has been calculated, compute the guiding center
        % % averaged angle of the guiding center position for cylindrical
        % % geometry.
            phi_gc(i) = trapz(t(cycle_indices),phi_c)/tgc(i);
        end

    % ideally, put some code here in case neither of the above cases are 
    % selected.
   
        % % compute the guiding center drifts
        vxgc(i-1)=(xgc(i)-xgc(i-1))/(tgc(i));
        vygc(i-1)=(ygc(i)-ygc(i-1))/(tgc(i));
        % % compute gyro-averaged magnetic field, for grad-B drifts
        Bgc(i-1)=trapz(t(cycle_indices),Bc(cycle_indices))/tgc(i);
    
    % % Alternative method, which lines up with my theory for abrupt
    % % inhomogeneity:
    %vxgc(i-1)=trapz(t(cycle_indices),vx(cycle_indices))/(tgc(i));
    %vygc(i-1)=trapz(t(cycle_indices),vy(cycle_indices))/(tgc(i));
    
        % % for cyl. coordinates:
        vrgc(i-1)=(rgc(i)-rgc(i-1))/(tgc(i));
	% % the phi component of the guiding center drift is tricky, because 
	% % rgc might be changing so the usual vphi=r*dphi/dt may give us wrong 
	% % answers; just use an average value of r for now.
        vphi_gc(i-1)=(rgc(i)+rgc(i-1))*(phi_gc(i)-phi_gc(i-1))/(tgc(i))/2;
	
    %disp(i)
    end
% WARNING: I THINK MAYBE THE FIRST ELEMENT OF THESE QUANTITIES SHOULD NOT
% BE DROPPED??? 4/15/2014
% % The first element of vxgc, vygc, xgc, ygc was JUST used to help us 
% % calculate some things; discard these now.
%xgc(1)=[];ygc(1)=[];vxgc(1)=[];vygc(1)=[];
% % likewise, discard first elements of the cyl. arrays.
%rgc(1)=[];phi_gc(1)=[];vrgc(1)=[];vphi_gc(1)=[];
% % turn xc_temp and yc_temp into the finished arrays xc and yc:
    xc=xc_temp;
    yc=yc_temp;
    Bc=Bc_temp;
    rc=rc_temp;
    phi_c=phi_c_temp;
    % done with the temporary arrays.
    clear xc_temp; clear yc_temp; clear Bc_temp;
    clear rc_temp; clear phi_c_temp;
% % Need to fix all of the *_polar arrays, because they are initialized to
% % zero for elements 1:new_cycle(1)-1.
    q_polar(1:new_cycle(1)-1)=[];
    %I_polar(1:new_cycle(1)-1)=[];
    phase_polar(1:new_cycle(1)-1)=[];

    % % Optional plotting commands to make sure everything is working 
    % % correctly.
    figure(2);clf;
    plot(x,y,'-b','linewidth',2);
    hold on;
    plot(xc,yc,'-k');
    plot(xgc,ygc,'sk','markerfacecolor','k','markersize',12);
    axis square;
    hold off;
    set(gcf, 'Color', [1,1,1]);
    % figure(3);clf;
    % plot(vxgc,'b');
    % hold on;
    % plot(vygc,'--r');
    % hold off;
    % set(gcf, 'Color', [1,1,1]);

    % % Comment this out when necessary.
    %polar(phase_polar,q_polar);

end
% % It takes a lot of time now for gyrophaser.m to run, so make sure the 
% % file, including the original filename gets saved!
save(strcat(filename,'_gyrophased.mat'))
% % don't need the old file anymore.
delete(strcat(data_file,'.mat'));

end