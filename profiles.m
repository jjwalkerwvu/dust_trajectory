% % profiles.m

function [V,E_x,E_y,B,vi_x,vi_y,ve_x,ve_y,vn_x,vn_y,g_x,g_y,ni,ne,alph,...
    Ti,Te,nneut,lambda_i,lambda_D,corot_period]=...
    profiles(Ti0,Te0,n0,t,x,y,profile_type,P,species)
% Perhaps add statements for Z, eta, tau, and ion/electron larmor radii here so 
% that they do not have to be computed elsewhere? 


% explanation of inputs:
% Ti0           = initial or equilibrium ion temperature in eV
% Te0           = initial or equilibrium electron temperature in eV
% n0            = equilibrium plasma density in m^(-3)
% t             = time, in seconds; this is only used if the profile has a 
%                 time-dependent inhomogeneity, like a pulsating UV light 
%                 (see 'uv_time' profile)
% x             = x-coordinate, in m
% y             = y-coordinate, in m
% profile_type  = a string specifying which profile you want to use; scroll 
%                 through this .m file to find one you want to use, or make your 
%                 own and put it in here
% P             = equilibrium neutral gas pressure, in pascals
% species       = the atomic mass number for the ions and neutrals in the plasma
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Explanation of outputs:
% V             = local space potential, in Volts
% E_x           = local electric field in x-direction, in V/m
% E_y           = local electric field in y-direction, in V/m
% B             = local magnetic field in z-direction, in T
% vi_x          = local ion flow velocity in x-direction, in m/s
% vi_y          = local ion flow velocity in y-direction, in m/s
% ve_x          = local electron flow velocity in x-direction, in m/s; for 
%                 situations where the grain is assumed to be levitating in the 
%                 sheath edge, ve_x is reserved for vi_z
% ve_y          = local electron flow velocity in y-direction, in m/s
% vn_x          = local neutral flow velocity in x-direction, in m/s
% vn_y          = local neutral flow velocity in y-direction in m/s
% g_x           = local acceleration due to gravity in x-direction, m/s^2
% g_y           = local acceleration due to gravity in y-direction, m/s^2
% ni            = local ion density in m^(-3)
% ne            = local ion density in m^(-3)
% alph          = local coefficient of UV illumination; f_uv/ne/vthe
% Ti            = local ion temperature in eV
% Te            = local electron temperature in eV
% nneut         = local neutral density in m^(-3)
% lambda_i      = local ion-neutral charge exchange mean free path, in m
% lambda_D      = local Debye length, in m
% corot_period  = time for planetary body to rotate once; this is the period of
%                 the plasma frame that co-rotates with the planetary body

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % I've decided to get rid of global vars; they are commented if you feel 
% % like using them again.
%global qe;
%global me;
%global mp;
%global eps0;
%global mi
qe=1.6e-19;
me=9.1e-31;
mp=1.67e-27;
eps0=8.854e-12;
mi=species*mp;
mr=me/mi;
%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%	All of the information below should go into each profile!

% % assume Tn = Ti, may not always be true. I will eventually need Tn0 as 
% % an input, and Tn as an output of profiles.m
nneut=P/qe/Ti0;
% % bit of code I am trying out to change pressure with respect to time.
% if t>=3
%     Pmtorr=.05;
%     nneut=(Pmtorr/7.5)/qe/Ti0;
% end

Tau=Te0/Ti0;

vthe=sqrt(2*qe*Te0/me);	% local electron thermal speed, m/s
vthi=sqrt(mr/Tau)*vthe;		% local ion (proton) thermal speed, m/s

% % calculate the mean free path; assuming this is the mean free path for 
% % neutrals. To save myself the headache of rewriting much of this code, I 
% % will follow Lampe 2003 pop, and set lambda_i = lambda_mfp = 
% % lambda_neut. 5.8e-19 is the cross-section for resonant charge exchange 
% % that I got out of Smirnov's book on ionized gases. This value is for a  
% % 1 ev beam of Ar ions incident on Ar neutrals, but according to Tsendin   
% % this energy dependence is weak so I'm using this value. This will have
% %  to be changed for different species! IT SHOULD BE: 
% % vthn=sqrt(2*qe*(Tn0=Ti0)/(mi=mn)), CORRECT???
% The definition below is from Tsendin's textbook:
nu=32*(5.8e-19)*nneut*sqrt(qe*Ti0/mi)/3/sqrt(pi);
if nu==0
	lambda_i=inf;
else
	lambda_i=vthi/nu;	% lambda_i will be infinite if nu=0!
    
    % This is an expression from Pascal Chabert's textbook, I think it is
	% not as correct as the Tsendin expression above.
    %lambda_i=1/5.8e-19/nneut;
end
% % error checking code below doesn't work, don't know why.
% if ch_model=='kortshagen'&&P==0
% 	exception = 'YOU MUST USE A POSITIVE DEFINITE VALUE FOR THE PRESSURE 
%   WHEN USING KORTSHAGEN CHARGE MODEL!!!';
% 	error(exception)
% end



% % if no profile is specified, then use a uniform profile
% if length(profile_type)==0;
%    disp('no profile was selected, so a uniform profile has been chosen');
%    profile_type='uniform';
% end

% % 	All the code above needs to go into each profile!
% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
switch profile_type
	case 'uniform'
        alph0=0;
		V=0;
		E_x=-100;
		E_y=0;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
	% % uniform plasma, but with an electric field
    case 'non_uniform_n0'
        alph0=0;
		V=0;
		E_x=0;
		E_y=0;
        % "L" marks where the density reaches (1-exp(-1))*n0: 
        L=0.01;
		ni=n0*sech(x/L);
		ne=ni;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
	% % uniform plasma, but with an electric field
	case 'uniform_E'
        alph0=0;
        L=1;	%% scale length of plasma
		V=0;
		%E_x=-Te0/2/L;
        E_x=-100;
		E_y=0;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'ExB_theory_simulation'
        alph0=0;
        L=1;	%% scale length of plasma
		V=0;
        E_x=-100;
		E_y=-100;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'ExB_theory_simulation_ion_drag'
        alph0=0;
        L=1;	%% scale length of plasma
		V=0;
        E_x=-0;
		E_y=-0;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=-200;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		B=1;
        % corot_period is unused for this profile.
        corot_period=0;
  	case 'ExB_theory_simulation_drag'
        alph0=0;
        L=1;	%% scale length of plasma
		V=0;
        E_x=-100;
        %E_x=0;
		E_y=-0;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=-100;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		B=4;
        % corot_period is unused for this profile.
        corot_period=0;     
    case 'grad_B_slab'
        alph0=0; 
		V=0;
        E_x=-0;
		E_y=-0;
        B0=4;
        % length scale:
        L=1000;	
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=-100;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		%B=4;
        B=B0-x/L*B0;
        % corot_period is unused for this profile.
        corot_period=0;
	% % uniform plasma, but with electric field only on for -1<E<1:
	case 'uniform_E_step'
        alph0=0;
		V=0;
		E_x=-10*(x>-.1)*(x<.1);
		E_y=0;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		Ti=Ti0;Te=Te0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'density_step'
        alph0=0;
        V=0;
        E_x=0;
        E_y=0;
        if x<0;
            ne=2*n0;
            ni=n0;
       	else
            ne=n0;
            ni=2*n0;
            
        end
        
        alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		Ti=Ti0;Te=Te0;
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % corot_period is unused for this profile.
        corot_period=0;
    % % for an E*L*sech(x/L) plasma potential:
    % % Also bear in mind that if you input E=positive value, you get a 
    % % positive E-field peak at x=-L, and a negative peak at x=L.
	case 'E_tanh_sech'
        alph0=0;
		E=100; 	%% max strength of E-field in V/m
		L=0.1;
		ne = n0*exp(E*L*sech(x/L)/Te);
		ni = n0*exp(-E*L*sech(x/L)/Ti);
        % % multiply by 2 here for ions and electrons??
        % % sept 2012 note to question above: yes.
		E_q = -(2)*(E*tanh(x/L)*sech(x/L)); 	  
        alph=alph0;                                        
        g_x=0;g_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    % % for an infinitely big electric field in -x direction:
    case 'E_lin_inf'
% 		if x<0.5&&x>0
% 			ne=n0*exp(-E*x/Te)*(x<0.5)*(x>0);
% 		else
% 			ne=n0;
% 		end
% 		ni=n0;
        alph0=0;
        E=100; %% max strength of E-field in V/m
        E_x = -E;
        E_y=0;
        V = E*x;
        ne=n0;
        ni=n0;
        alph=alph0;
        g_x=0;g_y=0; 
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    % % 2013/03/01: This profile is not a valid model.
    case 'V_step'
    % % I get a dirac delta function when x = 0; not sure how to adjust
    % % the electric field. I set it to zero for now.
        alph0=0;
        E=-100; %% max strength of E-field in V/m
        E_x = E;
        E_y=0;
        V = (x<0)*(-50);    %% 100 is just an arbitrary number in volts.
        ne=n0;
        ni=n0;
        alph=alph0;
        g_x=0;g_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'E_tanh'
        alph0=0;
		n1=n0;
		n2=n0*3;
		L=0.001;
		ne=n1+(n2-n1)/2*tanh(x/L);
		ni=ne;
		E_q=0;
        alph=alph0;
        g_x=0;g_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'B_step'
         alph0=0;
         E=0;
         E_x=E;
         E_y=E;
         V=0;
         ne=n0;
         ni=n0;
         alph=alph0;
         g_x=0;g_y=0;
         vi_x=0;vi_y=0;
         ve_x=0;ve_y=0;
         vn_x=0;vn_y=0;
         Ti=Ti0;Te=Te0;
         lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
         % Make the step function right at x=0
         %B=4*(x>0)+5*(x<=0);
         % should I reverse the step function, so that the greater negative  
         % charge (q1) is for x<=0, and the lesser negative charge (q2) is 
         % for x>0?
         B=4.05*(x<=0)+5*(x>0);	% B=4.05 T is chosen, with 1.6 micron dia. 
                                % grains in mind to demonstrate stationary 
                                % guiding center.
         % corot_period is unused for this profile.
        corot_period=0;
    case 'graded_uv'
        alph0=0.5;
        ni=n0;
        ne=n0;
        %% scale length:
        L=.075;
        if x<L&&x>-L
            	alph=(alph0*x+alph0/2)*(x<L)*(x>-L);
        elseif x>L
           		alph=alph0;
        else
            	alph=0;
        end
     
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'uv_step'
        %alph0=0.5;
        alph0=0.25;
        ni=n0;
        ne=n0;
        alph=(x>0)*alph0;
        V=0;
        E_x=-4;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
     
     case 'uv_step_sheath'
        %alph0=0.5;
        % boundary region, where it switches from UV present/absent
        %xb=-0.05;
        % default value should be xb=0.
        xb=0;
        Ti=Ti0;Te=Te0;
        ni=n0;
        ne=n0;
        % from Dove et. al.
        f_uv=1e17;
        % the following expression ensures that we have 20 times as much uv
        % flux as Dove et. al.
        alph0=20*f_uv/ne/vthe;
        %alph0=0.25;
        alph=(x>xb)*alph0;
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_y=0;
        % assume the bohm speed for ions at sheath for now, check this 
        % later. remember that ve_x is vi_z!!!
        ve_x=sqrt(qe*Te0/mi);
        vn_x=0;vn_y=0;
        
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        %lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % because ions flow at the bohm speed:
        lambda_D=sqrt((eps0*Te)/qe/(ne+ni));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;   
        
        
    % Just like 'uv_step' above, except now the UV gets modulated in time.  
    case 'uv_time'
        % figure out what frequency you want to modulate the UV. I suggest
        % using a frequency that is close to the dust cyclotron frequency
        % for the "shadowed" condition. This means you have to input all of
        % the necessary parameters before you pick a frequency. Likewise,
        % you only need to compute this timescale once.
        
        % NOTE: PICK A REASONABLE FLUX OF PHOTOELECTRONS, AND REDO ALL OF
        % THE STUFF BELOW FOR AN ECR PLASMA!!
        % UV profile, need to calculate appropriate alph0 again, taking
        % into account distance of the grain from the sun. calculated in a
        % seemingly weird way due to how the photo-current is computed in
        % the various charge models.
        solar_distance = 9.5;   % distance from the sun in AU
        effic = 1;  % conductors have efficiencies of ~1, oxides have ~0.1
        % this uv photon flux is for solar radiation and regolith
        f_UV = 2.8e13*effic/solar_distance.^2;
        % coefficient for solar radiation and regolith; consider redo-ing 
        % for water
        alph0 = 0.25*sqrt(4*pi)*f_UV/n0/vthe;
        % or, just input a value.
        alph0=0.25; % alph0=1/4 implies that you get a flux of photo
                    % electrons that is equal to the electron current that
                    % would be present if the grain were at the local space
                    % potential.
        
        % assume singly ionized plasma for now.
        Z=1;
        ni=n0;
        ne=n0;     
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % strongest magnetic field along axis for an ECR is 0.0875 T.
        B=0.0875;
        if t==0
            % Input what size of dust you are interested in, to find out
            % what frequency to modulate the UV.
            dust_size=15*1e-9;  % I've gone with 15 nm radius grains.
            C=4*pi*eps0*dust_size;
            alph=0;
            qflag=1;
            ch_model='oml';     % oml for now, but you are free to pick a 
                                % different charge model
            % assume for now that the grain moves at the neutral flow 
            % speed; might want to think about this assumption.
            w_grain=[sqrt(ve_x^2+ve_y^2) sqrt(vi_x^2+vi_y^2)];
            [Itot,q,Kn_R0,P0,P1,Pg1]=...
                charging_models(qflag,ch_model,dust_size,alph,Te,Ti,...
                ne,ni,B,Z,C,0,lambda_D,lambda_i,w_grain,species);
            % frequency at which to modulate UV:
            rho=1e3;
            md=4/3*pi*rho*dust_size^3;
            T_UV=abs(2*pi*md/q/B);
        end
        % YOU HAVE TO SHAPE THE PULSE, OTHERWISE THE GRAIN GOES AROUND IN A
        % CIRCLE!
        T1=0.1476;
        T2=.3095;
        %T_UV=0.1476;
        T_UV = T1/2+T2/2;
        % integer multiples of T_UV: 
        n_whole = floor(t/T_UV);
        
        %alph=alph0*(((n_whole+.5)*T_UV<t) && (n_whole+1)*T_UV>t);
        % I think this bottom line works:
        alph=alph0*((n_whole*T_UV+T1/2<t) && (n_whole+1)*T_UV>t);
        % corot_period is unused for this profile.
        corot_period=0;    
    case 'uv_gradient'
        alph0=0.25;
        ni=n0;
        ne=n0;
        alph=(x>0)*alph0;
        % % use region where q=q2=q1/2 as the scale length.
        %L=.0015;
        % % alternatively, use a smaller gradient:
        L=0.00015
        % alph increases from 0 at x=0 to a max value at x=L, and stays at 
        % this max value for x>L.
        alph=alph0*(x>=L) + alph0*x/L*(x>=0)*(x<=L);
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % corot_period is unused for this profile.
        corot_period=0;
        
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %	The child langmuir sheath provides the parameters for a dust grain 
% %	suspended in a simple DC discharge sheath. An electrode is placed at 
% %	the bottom, with an applied voltage relative to the plasma potential at   
% % the center of the discharge (space potential is taken as 0 in the 	
% %	middle of the discharge. This requires cold ions, or Te>>Ti. The sheath 
% % is oriented in the y-direction, with the planar surface starting at 
% % y=0.
 	case 'cl_sheath'
        % length of the plasma, in meters:
        L=.05;
        % ionization constant for Argon gas, in units of m^3 s^-1:
        Ki=5e-14;
        alph0=0;
 		V0=-10;     % % This is the voltage of the Sheath electrode 
                    % % relative to the plasma's space potential far away 
                    % % from the electrode
        % temp profile??
        Ti=Ti0;Te=Te0;
 		% % specify the bohm speed, either here or in the main program
 		ub=sqrt(qe*Te0/mi);
 		ns=n0*exp(-1/2);
 		if V0>Te0;
        % % if the above condition is false, CL condition does not hold.
        % % Note that V0 will probably be negative relative to Te0.
            %return
            exception='V0>Te, CL condition=false';
            error(exception);
 		end
 		%ji=qe*n0*sqrt(qe*Te0/mi);
 		ji=qe*ns*ub;
        % % recheck this!
 		%V=((4/3)*(2*x*(ji/eps0)^(1/2)*(2*qe/mi)^1/4-(-V0^3/4))^(4/3);
        % % calculate the sheath boundary. Might be necessary to include
        % % the -mi/2/qe*(ub).^2 term, required for the space potential to 
        % % be continuous.
 		s=2*((-V0)^(3/4))*sqrt(eps0/ji)*((2*qe/mi)^(1/4))/3;	
        % plasma density at the sheath boundary:
        ns=n0*exp(-1/2);
 		% % y=0 is the planar electrode which is well defined. Set the 
 		% % sheath boundary (s) above this. See Chabert's 2011 textbook.
 		% % What follows in this if statement is the plasma solution.
 		if y>=s;
 			% % for y>sheath edge location, the plasma is quasineutral. One
 			% % still needs to obtain the n0(y), however.

            % % Uncomment the lines below to get rid of complications:
            V=0;
            % electric field in the presheath is half an electron
            % temperature. This field must be present to accelerate ions to
            % the bohm speed at the plasma sheath boundary. I have this
            % field going in the negative direction. The potential drop is
            % Te/2, so you have to divide by the length of the presheath
            % (the size of the plamsa)
            %E_y=-Te/L;
            % No ion flow in x direction, although a cylindrical glow
            % discharge might have a radial ion flow.
            vi_x=0;
            % consider using the Tonks Langmuir solution!
            % first, compute the sheath length:
            %ys=2/3*ub/nneut/Ki;
            % or, just use:
            ys=L/2-s;
            ynorm=2/3*(L/2-y)/ys;
            % ion flow is given by the commented equation below, but must
            % be rewritten so that matlab can handle it.
            %vi_y=2*ub*cos(1/3*(4*pi-atan(sqrt(4/9/ynorm/ynorm-1))));
            % the rewritten version.
            vi_y=2*ub*cos(4*pi/3)*cos(atan(sqrt(4/9/ynorm/ynorm-1))/3)-...
            	2*ub*sin(4*pi/3)*sin(atan(sqrt(4/9/ynorm/ynorm-1))/3);
            %vi_y=0;
            V=-mi/2/qe*vi_y.^2+mi/2/qe*ub^2;
            % See my notebook #7, page 54 for the derivation. Still needs
            % to be debugged! May 13, 2014
%             E_y=(-2/3/ys)*mi/qe*vi_y*(8*ub/27/sqrt(4/9/ynorm.^2-1)/...
%                 (1+(4/9/ynorm.^2-1).^2)/ynorm.^3)*...
%                 (-sqrt(3)/2*cos(1/3*atan(sqrt(4/9/ynorm.^2-1)))+...
%                 +.5*sin(1/3*atan(sqrt(4/9/ynorm.^2-1))));
         	E_y=(-2/3/ys)*mi/qe*vi_y*(2*ub/3/ynorm/sqrt(4/9/ynorm.^2-1))...
                *(-sqrt(3)/2*cos(1/3*atan(sqrt(4/9/ynorm.^2-1)))+...
                +.5*sin(1/3*atan(sqrt(4/9/ynorm.^2-1))));
            %E_y=0;
            ni=n0*exp((V-mi/2/qe*ub^2)/Te);
            ne=ni;
            E_x=0;
            %ne=n0;ni=n0;
            % use boundary values? Or find the plasma solution, otherwise
            % the densities are discontinuous.
            %ne=ns;ni=ns;
            ve_x=0;ve_y=0;
            vn_x=0;vn_y=0;
 		% below the boundary.	
 		else
 			% I added the extra -mi/2/qe*(ub).^2 term so that the 
 			% potentials match on either side of the boundary! Check this 
            % in the future.
 			V=-(-3/2*sqrt(ji/eps0)*(2*qe/mi)^(-1/4)*y+...
                (-V0)^(3/4)).^(4/3);
 			E_y=(-2*sqrt(ji/eps0)*(2*qe/mi)^(-1/4))*...
                (-3/2*sqrt(ji/eps0)*(2*qe/mi)^(-1/4)*y+(-V0)^(3/4)).^(1/3);
           	% I take away the term -mi/2/qe*ub.^2 so that density is 
            % continuous.
            ne=ns*exp(V/Te0);
            %ni=ji*sqrt(-mi/2/qe/V)/qe 	% huge problems when V~0 at sheath 
                                        % edge!                             
            ni=ns*(1-2*qe*(V)/mi/ub^2).^(-1/2);
            % specify the ion flow at this spatial location, assuming the 
            % sheath is located at y=s. Also, ions should be flowing
            % downward:
            vi_y=-ns*ub/ni;
            E_x=0;
            vi_x=0;
            % % Uncomment to get rid of ion streaming
            %vi_y=0;
            % electron streaming?
            ve_x=0;ve_y=0;
            % should there be neutrals streaming??
            vn_x=0;vn_y=0;
 		end
 	
 		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
 		% % specify acceleration due to gravity; can put this in x or y 
 		% % direction but generally g_y=9.8 m/s/s on Earth experiments
 		g_x=0;
 		g_y=-9.8;
        
        % % try no magnetic field first.
 		B=0;
        %B=4;
        % UV illumination??
        alph=alph0;
        % corot_period is unused for this profile.
        corot_period=0;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %	cl_sheath for RF plasma. Still need to add in the self bias, April 10,
% % 2014. This is applicable for RF frequencies that are above the ion
% % plasma frequency, but below the electron plasma frequency.
 	case 'cl_sheath_rf'
        % the rf amplitude:
        Vrf=50;
        % the floating potential of a planar sheath when a sinusoidal RF
        % voltage is applied:
        Vself=Te0*(.5*log(2*pi*me/mi)-log(besseli(0,Vrf/Te0)));
        % alternatively, for a square wave instead of sinusoidal:
        %Vself=Te0*(.5*log(2*pi*me/mi)-log(cosh(Vrf/Te0)));
        
        % ionization constant for Argon gas, in units of m^3 s^-1:
        Ki=5e-14;
        alph0=0;
 		V0=-50;     % % This is the voltage of the Sheath electrode 
                    % % relative to the plasma's space potential far away 
                    % % from the electrode
        % temp profile??
        Ti=Ti0;Te=Te0;
 		% % specify the bohm speed, either here or in the main program
 		ub=sqrt(qe*Te0/mi);
 		ns=n0*exp(-1/2);
 		if V0>Te0;
        % % if the above condition is false, CL condition does not hold.
        % % Note that V0 will probably be negative relative to Te0.
            %return
            exception='V0>Te, CL condition=false';
            error(exception);
 		end
 		%ji=qe*n0*sqrt(qe*Te0/mi);
 		ji=qe*ns*ub;
        % % recheck this!
 		%V=((4/3)*(2*x*(ji/eps0)^(1/2)*(2*qe/mi)^1/4-(-V0^3/4))^(4/3);
        % % calculate the sheath boundary.
 		s=2*((-V0)^(3/4))*sqrt(eps0/ji)*((2*qe/mi)^(1/4))/3;	
        % plasma density at the sheath boundary:
        ns=n0*exp(-1/2);
 		% % y=0 is the planar electrode which is well defined. Set the 
 		% % sheath boundary (s) above this. See Chabert's 2011 textbook
 		if y>=s;
 			% % for y>sheath edge location, the plasma is quasineutral. One
 			% % still needs to obtain the n0(y), however.
            
            % the stuff below is really not appropriate for a CL sheath; it
            % is more relevant to the plasma solution.
%             % Specify the ion flow at this spatial location, assuming the 
%             % sheath is located at y=s:
%             vi_y=2*ub*cos(1/3*(4*pi-sqrt(atan(4/9/(2*y/s)^2-1))));
%  			V=-mi/2/qe*vi_y.^2;
%  			E_y=-2*mi/qe*vi_y*ub;
%             ni=n0*exp(V/Te);
%             ne=ni;
 
            % % Uncomment the lines below to get rid of complications:
            V=0;
            % electric field in the presheath is half an electron
            % temperature. This field must be present to accelerate ions to
            % the bohm speed at the plasma sheath boundary. I have this
            % field going in the negative direction.
            E_y=-Te/2;
            vi_x=0;
            % ions are at the bohm speed everywhere in the presheath, and
            % stream 
            vi_y=-ub;  
            E_x=0;
            ne=n0;ni=n0;
            ve_x=0;ve_y=0;
            vn_x=0;vn_y=0;
 		% below the boundary.	
 		else
 			V=-(-3/2*sqrt(ji/eps0)*(2*qe/mi)^(-1/4)*y+(-V0)^(3/4)).^(4/3);
 			E_y=(-2*sqrt(ji/eps0)*(2*qe/mi)^(-1/4))*...
                (-3/2*sqrt(ji/eps0)*(2*qe/mi)^(-1/4)*y+(-V0)^(3/4)).^(1/3);
            ne=ns*exp(V/Te0);
            %ni=ji*sqrt(-mi/2/qe/V)/qe 	% huge problems when V~0 at sheath 
                                        % edge!
            ni=ns*(1-2*qe*V/mi/ub^2).^(-1/2);
            % specify the ion flow at this spatial location, assuming the 
            % sheath is located at y=s. Also, ions should be flowing
            % downward:
            vi_y=-ns*ub/ni;
            E_x=0;
            vi_x=0;
            % % Uncomment to get rid of ion streaming
            %vi_y=0;
            % electron streaming?
            ve_x=0;ve_y=0;
            % should there be neutrals streaming??
            vn_x=0;vn_y=0;
 		end
 	
 		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
 		% % specify acceleration due to gravity; can put this in x or y 
 		% % direction but generally g_y=9.8 m/s/s on Earth experiments
 		g_x=0;
 		g_y=-9.8;
        
        % % try no magnetic field first.
 		B=0;
        B=4;
        % UV illumination??
        alph=alph0;
        % corot_period is unused for this profile.
        corot_period=0;   
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %	The child langmuir sheath provides the parameters for a dust grain 
% %	suspended in a simple DC discharge sheath. An electrode is placed at 
% %	the bottom, with an applied voltage relative to the plasma potential at   
% % the center of the discharge (space potential is taken as 0 in the 	
% %	middle of the discharge. This requires cold ions, or Te>>Ti. The sheath 
% % is oriented in the y-direction, with the planar surface starting at 
% % y=0.
 	case 'cl_sheath_collisional'
        % length of the plasma, in meters:
        L=.05;
        % ionization constant for Argon gas, in units of m^3 s^-1:
        Ki=5e-14;
        alph0=0;
 		V0=-10;     % % This is the voltage of the Sheath electrode 
                    % % relative to the plasma's space potential far away 
                    % % from the electrode
        % temp profile??
        Ti=Ti0;Te=Te0;
 		% % specify the bohm speed, either here or in the main program
 		ub=sqrt(qe*Te0/mi);
 		ns=n0*exp(-1/2);
 		if V0>Te0;
        % % if the above condition is false, CL condition does not hold.
        % % Note that V0 will probably be negative relative to Te0.
            %return
            exception='V0>Te, CL condition=false';
            error(exception);
 		end
 		%ji=qe*n0*sqrt(qe*Te0/mi);
 		ji=qe*ns*ub;
        % % recheck this!
 		%V=((4/3)*(2*x*(ji/eps0)^(1/2)*(2*qe/mi)^1/4-(-V0^3/4))^(4/3);
        % % calculate the sheath boundary. Might be necessary to include
        % % the -mi/2/qe*(ub).^2 term, required for the space potential to 
        % % be continuous.
 		%s=2*((-V0)^(3/4))*sqrt(eps0/ji)*((2*qe/mi)^(1/4))/3;	
        s=(-5/3*V0)^(3/5)*(2/3*qe*ns/eps0*sqrt(pi*Te/2/lambda_i))^(-2/5);
        % plasma density at the sheath boundary:
        ns=n0*exp(-1/2);
 		% % y=0 is the planar electrode which is well defined. Set the 
 		% % sheath boundary (s) above this. See Chabert's 2011 textbook.
 		% % What follows in this if statement is the plasma solution.
 		if y>=s;
 			% % for y>sheath edge location, the plasma is quasineutral. One
 			% % still needs to obtain the n0(y), however.

            % % Uncomment the lines below to get rid of complications:
            V=0;
            % electric field in the presheath is half an electron
            % temperature. This field must be present to accelerate ions to
            % the bohm speed at the plasma sheath boundary. I have this
            % field going in the negative direction. The potential drop is
            % Te/2, so you have to divide by the length of the presheath
            % (the size of the plamsa)
            %E_y=-Te/L;
            % No ion flow in x direction, although a cylindrical glow
            % discharge might have a radial ion flow.
            vi_x=0;
            % consider using the Tonks Langmuir solution!
            % first, compute the sheath length:
            %ys=2/3*ub/nneut/Ki;
            % or, just use:
            ys=L/2-s;
            ynorm=2/3*(L/2-y)/ys;
            % ion flow is given by the commented equation below, but must
            % be rewritten so that matlab can handle it.
            %vi_y=2*ub*cos(1/3*(4*pi-atan(sqrt(4/9/ynorm/ynorm-1))));
            % the rewritten version.
            vi_y=2*ub*cos(4*pi/3)*cos(atan(sqrt(4/9/ynorm/ynorm-1))/3)-...
            	2*ub*sin(4*pi/3)*sin(atan(sqrt(4/9/ynorm/ynorm-1))/3);
            %vi_y=0;
            %V=-mi/2/qe*vi_y.^2+mi/2/qe*ub^2;
 
            ni=n0*sqrt(1-(2*(y-s)/L).^2);
            ne=ni;
            V= Te0*log(ne/ns);
            E_y=0;
            E_x=0;
            %ne=n0;ni=n0;
            % use boundary values? Or find the plasma solution, otherwise
            % the densities are discontinuous.
            %ne=ns;ni=ns;
            ve_x=0;ve_y=0;
            vn_x=0;vn_y=0;
 		% below the boundary.	
 		else
 			% I added the extra -mi/2/qe*(ub).^2 term so that the 
 			% potentials match on either side of the boundary! Check this 
            % in the future.
 			V=-3/5*(2/3*qe*ns/eps0*sqrt(pi*Te/2/lambda_i))^(2/3)*...
                (s-y)^(5/3);
 			E_y=-(2/3*qe*ns/eps0*sqrt(pi*Te/2/lambda_i))^(2/3)*...
                (s-y)^(2/3);
           	% I take away the term -mi/2/qe*ub.^2 so that density is 
            % continuous.
            ne=ns*exp(V/Te0);
            %ni=ji*sqrt(-mi/2/qe/V)/qe 	% huge problems when V~0 at sheath 
                                        % edge!
            % I take away the term -mi/2/qe*ub.^2 so that density is 
            % continuous.                                  
            ni=ns*(1-2*qe*V/mi/ub^2).^(-1/2);
            % specify the ion flow at this spatial location, assuming the 
            % sheath is located at y=s. Also, ions should be flowing
            % downward:
            vi_y=-ns*ub/ni;
            E_x=0;
            vi_x=0;
            % % Uncomment to get rid of ion streaming
            %vi_y=0;
            % electron streaming?
            ve_x=0;ve_y=0;
            % should there be neutrals streaming??
            vn_x=0;vn_y=0;
 		end
 	
 		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
 		% % specify acceleration due to gravity; can put this in x or y 
 		% % direction but generally g_y=9.8 m/s/s on Earth experiments
 		g_x=0;
 		g_y=-9.8;
        
        % % try no magnetic field first.
 		B=0;
        %B=4;
        % UV illumination??
        alph=alph0;
        % corot_period is unused for this profile.
        corot_period=0;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
 	case 'linear_profile'
        alph0=0;
        %L=75;    %% the scaling for density gradient in m
        L=0.025;
        % % input the desired electric field here:
        %E_base=-100;
        E_base=0;
 		eta=n0/L;	%% density gradient??
 		% % the cutoff:
 		ratio=10;	% <<<---- factor by which one species is greater than 
                    % the other at the cutoff; to be determined.
 		% set fact=1 if you want ni=ne=n0 at the cuttoff.
        %cutoff=(ratio-1)/(ratio+1)*n0/eta;
        % same thing as the above:
        cutoff=(ratio-1)/(ratio+1)*L;
        % wouldn't this just be:
        %cutoff=(ratio-1)/(ratio+1)*L;
 		
 		%cutoff=L;
 		ni=n0+eta*x;	
 		ne=n0-eta*x;
 		% % To reverse the profile, comment the above 2 lines and uncomment 
        % % the following 2 lines:
 		%ne=n0+eta*x;
 		%ni=n0-eta*x;
 		V=-qe*(eta*x^3)/3/eps0;
 		% % put in the following lines to make sure the densities don't go 
        % % negative.
 		if x<=-cutoff;
 			
 			ni=n0+eta*cutoff;
 			ne=n0-eta*cutoff;
 			%ni=n0;
 			%ne=n0;
 			%disp('outside the cutoff');
 			
 			%% still need to fix the next two lines:
 			V0=-qe*(eta*cutoff^3)/3/eps0;
 			V=V0+qe*(ratio-1)*n0*x^2/(ratio+1)/eps0;
 			
 			
 		end
 		if x>=cutoff;
 	
 			ni=n0+eta*cutoff;
 			ne=n0-eta*cutoff;
 			%ni=n0;
 			%ne=n0;
 			%disp('outside the cutoff');
 			
 			% % still need to fix the next two lines:
 			V0=qe*(eta*cutoff^3)/3/eps0;
 			V=V0-qe*(ratio-1)*n0*x^2/(ratio+1)/eps0;
 			
        end
        % E=0 for now; can add linear and inhomog. profile later.
 		E_x=E_base*(x>-cutoff)*(x<cutoff);	
 		E_y=0;
        alph=alph0;
        vi_x=0;vi_y=0;
        % % not sure to do with the electron flow; set equal to zero for
        % % now.
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        %B=4;
        B=75;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'parabolic_temp'
        alph0=0;
        V=0;
		E_x=0;
		E_y=0;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
		L=0.01;	% half width at full max for this profile
		Ti=Ti0;	% ion profile is flat
		Tmax=20.0;	% "peak" electron temperature
		if abs(x)<=L
            % parabolic inside the inhomogeneity
			Te=Tmax-(Tmax-Te0)*(x/L)^2;	
		else
			Te=Te0;	%% flat outside the inhomogeneity
		end
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		B=4;
        % corot_period is unused for this profile.
        corot_period=0;
	case 'linear_temp'
        alph0=0;
        V=0;
		E_x=0;
		E_y=0;
		ni=n0;
		ne=n0;
		alph=alph0;
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		L=0.01;	% half width at full max for this profile
		Ti=Ti0;	% ion profile is flat
		T_L=1.0/L;	% electron temp goes "down" this amount in a length L; 
                    % units: eV/m
		
		if x<-L
			Te=T_L*L+Te0;
		end
		if abs(x)<=L
			Te=Te0+T_L*L/2-T_L*x/2;
		end
		if x>L
			Te=Te0;
        end
        % % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
		% corot_period is unused for this profile.
        corot_period=0;
    case 'linear_profile_w_ion_drift'
        alph0=0;
        %L=75;    %% the scaling for density gradient in m
        L=0.025;
        % % input the desired electric field here:
        %E_base=-100;
        E_base=0;
 		eta=n0/L;	%% density gradient??
 		% % the cutoff:
 		ratio=10;	% <<<---- factor by which one species is greater than 
                    % the other at the cutoff; to be determined.
        % set fact=1 if you want ni=ne=n0 at the cuttoff.
 		cutoff=(ratio-1)/(ratio+1)*n0/eta;		
 		
 		%cutoff=L;
 		ni=n0+eta*x;	
 		ne=n0-eta*x;
 		% % To reverse the profile, comment the above 2 lines and uncomment 
        % % the following 2 lines:
 		%ne=n0+eta*x;
 		%ni=n0-eta*x;
 		% % put in the following lines to make sure the densities don't 
        % % go negative.
 		if x<=-cutoff;
 			
 			ni=n0+eta*cutoff;
 			ne=n0-eta*cutoff;
 			%ni=n0;
 			%ne=n0;
 			%disp('outside the cutoff');
 			
 		end
 		if x>=cutoff;
 	
 			ni=n0+eta*cutoff;
 			ne=n0-eta*cutoff;
 			%ni=n0;
 			%ne=n0;
 			%disp('outside the cutoff');
 			
 		end
 		E_x=E_base*(x>-cutoff)*(x<cutoff);	% E=0 for now; can add linear 
                                            % and inhomog. profile later.
 		E_y=0;
 		V=0;
        alph=alph0;
        % % If electric field is along -x, then there should be an ion 
        % % drift along -x.
        vi_x=-sqrt(qe*Te/40/1.67e-19);
        vi_y=0;
        Ti=Ti0;Te=Te0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
        % % not sure what to do with electron flow; set to 0 for now.
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'uniform_E_cyl'
		phi=improved_arctan(x,y);
		% % radius of the experimental volume, in meters
		R=0.225;	
        % This corresponds exactly with the "central experimental volume"
        % on the Auburn machine; the smaller, "uniform region" radius is 
        % 0.1
		
		Er=-100;
        ni=n0;
        ne=n0;
        radius=sqrt(x^2+y^2);
        alph=0;
        V=0;
        E_x=Er*cos(phi);
        E_y=Er*sin(phi);
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        % % If you want to use a sheath:
        % assume the bohm speed for ions at sheath for now, check this 
        % later. remember that ve_x is vi_z!!!
        ve_x=sqrt(qe*Te0/mi);
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        Ti=Ti0;Te=Te0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
	case 'cylindrical_profile'
        B=4;
		% % radius of the experimental volume, in meters
		R=0.225;	
        % This corresponds exactly with the "central experimental volume"
        % % on the Auburn machine; the smaller, "uniform region" radius is 
        % % 0.1
		phi=improved_arctan(x,y);
		
		Er=100; %% radial electric field in V/m
		
		radius=sqrt(x^2+y^2);
		%%~!!! 	Attempting to put E-field in only a specific volume.
		%if radius>=0.225 || radius<=0.05
			%Er=0;
		%else 	
			%Er=Te0/2;
		%end
	
		
		E_x=Er*cos(phi);
		E_y=Er*sin(phi);
		
		% % don't worry about V for now.
		V=0;
	
		alph=0;
		% % assuming B is in the +z direction, the ion drifts should be 
		% % simply written as E/B, in the -phi direction.
		% % 4/18/2013: There should be an additional component, probably 
		% % radially inwardcorresponding to the ion flow of the charge 
        % % imbalance.
		vi_x=Er*sin(phi)/B;
		vi_y=-Er*cos(phi)/B;
        ve_x=Er*sin(phi)/B;
        ve_y=-Er*cos(phi)/B;
        vn_x=0;vn_y=0;
		% % Uncomment the line below to turn off ion drag
		vi_x=0;vi_y=0;
		% % Uncomment the line below to turn on an ion drag with no E-field 
		% % term. This will produce an outward, FxB drift
		%E_temp=1;vi_x=E_temp*sin(phi)/B;vi_y=-E_temp*cos(phi)/B;
		ne=n0;ni=n0;
		Ti=Ti0;Te=Te0;
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
	
        % corot_period is unused for this profile.
        corot_period=0;
    case 'linear_cylindrical_profile'
        % % This is like the linear_profile, but instead this is in 
        % % cylindrical geometry
        % % radius of the experimental volume, in meters
		R=0.225;	
        % This corresponds exactly with the "central experimental volume" 
        % on the Auburn machine; the smaller, "uniform region" radius is 
        % 0.1
		phi=improved_arctan(x,y);
        radius=sqrt(x^2+y^2);
		
        
        R1=0.15;
        R2=0.175;
        dR=R2-R1;   % R1 and R2 have been chosen such that dR=0.025 meters.
                    % keep in mind that the center of the gradient, where
                    % ne=ni is at R0=0.1625; use this as an initial
                    % condition in dust_trajectory.m
                    
        R0=R1+dR/2; % center of the density gradient
        ratio=10;   % this ratio signifies by how many times larger is the
                    % electron density than the ion density at the lower
                    % cutoff, and by what factor the ion density is larger 
                    % than the electron density at the upper cutoff.
        % lower radius cutoff:
        r_lower=(R1+(1-ratio)*dR+ratio*R2)/(1+ratio);
        % upper radius cutoff;
        r_upper=2*R0-r_lower;
        
        ne=-2*n0/dR*(radius-R1)+2*n0;
        ni=2*n0/dR*(radius-R2)+2*n0;
        V=0;
        %Er=0;
        Er=(2*qe*n0/eps0/dR)*(2/3*radius.^2-2/3*r_lower^3/radius-...
            (R1+R2)/2*radius+(R1+R2)/2*r_lower^2/radius);
        Er=0;
        % electron and ion densities are flat outside of the inhomogeneity
        if radius<r_lower
            % alternatively, set ne=ni=n0
            ne=-2*n0/dR*(r_lower-R1)+2*n0;
            ni=2*n0/dR*(r_lower-R2)+2*n0;
            ne=n0;ni=n0;
            V=0;
            %Er=(-4*radius+(R1+R2))*log(r_upper/r_lower)+...
            %   4*(r_upper-r_lower);
            Er=0;
        end
        % electron and ion densities are flat outside of the inhomogeneity
        if radius>r_upper
            % alternatively, set ne=ni=n0
            ne=-2*n0/dR*(r_upper-R1)+2*n0;
            ni=2*n0/dR*(r_upper-R2)+2*n0;
            ne=n0;ni=n0;
            V=0;
            Er=(2*qe*n0/eps0/dR)*(1/radius)*(2/3*(r_upper^3-r_lower^3)-...
                (R1+R2)/2*(r_upper^2-r_lower^2));
        end
        
		E_x=Er*cos(phi);
		E_y=Er*sin(phi);
		
		
		% % if this is uncommented, don't worry about E-field for now.
        %E_x=0;E_y=0;
        % % if this is uncommented, don't worry about drifts for now.
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        % % no temperature gradients or UV illumination
        Ti=Ti0;Te=Te0;alph=0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
	case 'two_gaussian_cylindrical_profile'
        phi=improved_arctan(x,y);
		% % radius of the experimental volume, in meters
		R=0.225;	
        % This corresponds exactly with the "central experimental volume"
        % on the Auburn machine; the smaller, "uniform region" radius is 
        % 0.1 m
		
        % This is the center of the gaussian E-field in radial direction 
        % (two peaks)
		r1=R/4;	
        % The center for the gaussian ion density in radial direction 
        % (also two peaks)
		r2=R/2;	
        r1=0.11;
        r2=0.14;
		
        % Width of the gaussian E-field in meters
		a1=0.03;
        % Width of the gaussian ion density "perturbation" to background
		a2=0.06;	
		
		E0=100;	% strength of the electric field in V/m
		radius=sqrt(x^2+y^2);
		% compute the radial drift component of ion velocity by the 
        % Einstein Relation
		mu_i=(1/Ti0)*(3*pi*vthi*lambda_i)/(16*sqrt(2));
        % see 2006 Reynolds phys. plasmas for my inspiration on this
		Er=E0*exp(-((radius-r1)/a1).^2);		
                                                 
		% see 2006 Reynolds phys. plasmas
		ni=n0*(1+2*exp(-((radius-r2)/a2).^2));	
						% for my inspiration on this. I put a factor
						% of 2 here
		% To get a self-consistent electron density, take a spatial 
		% derivative of Er; multiply by eps0/qe Then add ni.
		
		ne=(2*eps0/qe)*Er*(radius-r1)/a1/a1+ni;
		
		% Uncomment the line below to have a constant density profile.
		%ne=n0;ni=n0;	
		E_x=Er*cos(phi);
		E_y=Er*sin(phi);
		% don't worry about V for now.
		V=0;
	
		alph=0;
		% % assuming B is in the +z direction, the ion drifts should be 
		% % simply written as E/B, in the -phi direction.
		% % 4/18/2013: There should be an additional component, probably 
		% % radially inward corresponding to the ion flow of the charge 
		% % imbalance. This comes from ion mobility and this change has 
        % % hopefully been made correctly. (5/13/2013)
		vi_x=Er*sin(phi)/B+mu_i*E_x;	
                    % % I think the above expressions should be fine, 
					% % because the inhomogeneity scale length is much 
					% % larger than the ion gyro-radii.
		vi_y=-Er*cos(phi)/B+mu_i*E_y;
		% % Uncomment the lines below to turn off ion drag and Electric 
        % % field independently
		vi_x=0;vi_y=0;
        % % Not sure what to do with electron terms; temporarily set to
        % % zero.
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        %E_x=0;E_y=0;
		Ti=Ti0;Te=Te0;
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		B=4;
        % corot_period is unused for this profile.
        corot_period=0;
	case 'parabolic_cylindrical_profile'
		phi=improved_arctan(x,y);
		% % radius of the experimental volume, in meters
        % % This corresponds exactly with the "central experimental volume"
		R=0.225;	
		% % on the Auburn machine; the smaller, "uniform region" radius is 
        % % 0.1
        
        % % The radius at which E~0. If r>r_field, E=0, space potential is 
        % % constant.
		r_field=0.075;	
		
        % % the space potential at the edge of the parabolic profile.
		V_edge=2.5;	
        % % space potential when r=0. 
		V_max=25;		
		
		%E0=1e2;	%% strength of the electric field in V/m
		radius=sqrt(x^2+y^2);
		if radius<=r_field
			Er=-2*(V_edge-V_max)*radius/r_field/r_field;
			ne=n0;
			ni=n0;
			% % assuming B is in the +z direction, the ion drifts should be 
			% % simply written as E/B, in the -phi direction.
			% % 4/18/2013: There should be an additional component,
			% % probably radially inward corresponding to the ion flow of 
            % % the charge imbalance.
			vi_x=Er*sin(phi)/B;	% % I think these expressions should be 
                                % % fine, because the inhomogeneity scale 
                                % % length is much larger than the ion 
                                % % gyro-radii.
			vi_y=-Er*cos(phi)/B;
			% % Uncomment the line below to turn off ion drag
			%vi_x=0;vi_y=0;		
		else	
			Er=0;
			% % 4/18/2013: There should be an additional component, 
			% % probably radially inward corresponding to the ion flow of 
            % % the charge imbalance.
			vi_x=0;
			vi_y=0;
			ne=n0;
			ni=n0;
		end
		
		% % Uncomment the line below to have a constant density profile.
		%ne=n0;ni=n0;	
		E_x=Er*cos(phi);
		E_y=Er*sin(phi);
		% % don't worry about V for now.
		V=0;
	
		alph=0;
		Ti=Ti0;Te=Te0;
		% % gravity unimportant since it is not in the plane or too weak.
		g_x=0;g_y=0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % % not sure what to do with electron terms; temporarily set to
        % % zero.
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
	case 'uv_cyl_profile_ysection'
		phi=improved_arctan(x,y);
		% % radius of the experimental volume, in meters
        % % This corresponds exactly with the "central experimental volume"
        % % on the Auburn machine; the smaller, "uniform region" radius is 
        % % 0.1
		R=0.225;	
		
		% % try UV on for all y, but with a beam centered at x=0, and 
        % % -L<=x<=L
		L=0.05;
		alph0=0.25;
        ni=n0;
        ne=n0;
        alph=(x>-L)*alph0*(x<L);
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'uv_step_cyl_radial'
		phi=improved_arctan(x,y);
		% % radius of the experimental volume, in meters
		R=0.225;	% This corresponds exactly with the "central 
                    % experimental volume" on the Auburn machine; the 
                    % smaller, "uniform region" radius is 0.1
		
		% try UV on for all y, but with a beam centered at x=0, and 
        % -L<=x<=L
		L=0.05;
		alph0=0.25;
        ni=n0;
        ne=n0;
        radius=sqrt(x^2+y^2);
        alph=(radius<L)*alph0;
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    
        case 'uv_step_sheath_cyl_radial'
		phi=improved_arctan(x,y);
		% % radius of the experimental volume, in meters
		R=0.225;	% This corresponds exactly with the "central 
                    % experimental volume" on the Auburn machine; the 
                    % smaller, "uniform region" radius is 0.1
		
		% try UV on for all y, but with a beam centered at x=0, and 
        % -L<=x<=L
		L=0.05;
		alph0=0.0125;
        ni=n0;
        ne=n0;
        radius=sqrt(x^2+y^2);
        alph=(radius<L)*alph0;
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_y=0;
        % assume the bohm speed for ions at sheath for now, check this 
        % later. remember that ve_x is vi_z!!!
        ve_x=sqrt(qe*Te0/mi);
        vn_x=0;vn_y=0;
        Ti=Ti0;Te=Te0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    
    case 'cyl_stellarator'
        % obviously a stellarator has a much more complicated geometry, but
        % this is a simple attempt.
        phi=improved_arctan(x,y);
        radius=sqrt(x^2+y^2);
        R=0.6;  % % this corresponds to the beginning of the  
                % % sheath/SOL/density  drop-off in the Large Helical 
                % % Device
                
        % IN THIS CASE, IT IS THE TEMPERATURE PROFILES THAT ARE 
        % INHOMOGENEOUS!!
        Te_max=3.5e3;   %% Te at r=0 in eV.
        Ti_max=7e3;     %% Ti at r=0 in eV.
        T_edge=1e3;     %% Te~Ti at the edge.
        % Using a parabola for now, try gaussian or sech later.
        Te=Te_max-(Te_max-T_edge)*(radius/R).^2;
        Ti=Ti_max-(Ti_max-T_edge)*(radius/R).^2;
        if radius>R
            Te=T_edge;
            Ti=T_edge;
        end
        % the density profiles are nearly flat, provided r<=R.
        ni=n0;
        ne=n0;
        % fill in electric field and ion drift information later.
        V=0;
        E_x=0;
        E_y=0;
        vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
        alph=0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
        lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        B=4;
        % corot_period is unused for this profile.
        corot_period=0;
    
    case 'enceladus_cyl_uv_only' 
        % please note: this profile is set up to be used with
        % "corotating_boris_pusher.m" exclusively! Things will be different
        % if you work in a non-corotating frame.
        phi=improved_arctan(x,y);
        radius=sqrt(x^2+y^2);
        % radius of saturn, in meters, +/- 4000 m:
        r_sat=60268e3; 
		
        % space potential is set to zero, for now
		V=0;
        % assume homogeneous plasma, for now. For reference, n0~4e7 m^-3.
        % also, neutral density is 10^10 m^-3, convert this to Pascals
		ni=n0;
		ne=n0;

        % here's a stab at an electron density gradient: exponential up to
        % 10 saturn radii, but back down to a constant Temp after passing 
        % this radial distance.
        %Te = Te0*exp((radius-r_sat)/(7*r_sat))*(radius<10*r_sat)+...
        %    Te0*exp((9*r_sat)/(7*r_sat))*(radius>10*r_sat);
        Te=Te0;
        % assume no temperature gradients, for now. For reference, Te~1-10
        % eV near saturn, and Ti~10-20 eV?
        %Te=Te0;
        Ti=Ti0;
        
        % comment the line below if you don't want a log temp gradient
        %Te=log(1);
        % no flows in the co-rotating frame with saturn?
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		
		% compute the local debye length
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % check 2009? Farrell et. al. to get an estimate for B during 
        % closest approach of Cassini. Use homogeneous B-field at first,
        % but use a Dipole as simulations become more sophisticated. Keep
        % in mind also, that saturn's magnetic north pole is also its
        % geographical north pole, so if the field is a dipole the magnetic
        % field direction is pointing along the -z direction.
        %B=-3.4e-7;
        % optional: dipole field, z component only for now:
        B=4*pi*1e-7/4/pi*(-4.5793e25/(radius.^3));
        
        % We will work in Saturn's corotating frame, so use the rotation 
        % period of saturn, in seconds:
        corot_period=10.57*3600; % hours*seconds/hour
        % gravity of saturn is not time-dependent
        G = 6.67384e-11;    % gravitational constant
        m_sat = 5.6846e26;  % mass of saturn in kg
        m_enc = 1.08e20;    % mass of enceladus in kg
        m_rhea = 2.3e21;    % mass of rhea in kg
        % distance of grain from saturn
        
        g_sat_r = -G*m_sat/radius^2;
        %g_enc_r = -G*m_enc/
        g_x=g_sat_r*cos(phi);
        g_y=g_sat_r*sin(phi);
        
        % perterbation due to enceladus:
        x_enc = 0;
        y_enc = 0;
        g_enc_r = -G*m_enc/(sqrt((x_enc-x).^2+y_enc-y).^2);
        
        % UV profile, need to calculate appropriate alph0 again, taking
        % into account distance of the grain from the sun. calculated in a
        % seemingly weird way due to how the photo-current is computed in
        % the various charge models.
        solar_distance = 9.5;   % distance from the sun in AU
        effic = 1;  % conductors have efficiencies of ~1, oxides have ~0.1
        % this uv photon flux is for solar radiation and regolith
        f_UV = 2.8e13*effic/solar_distance.^2;
        % coefficient for solar radiation and regolith; consider redo-ing 
        % for water
        alph0 = 0.25*sqrt(4*pi)*f_UV/n0/vthe;  
        % OR: uncomment the line below to turn off UV charging.
        %alph0=0;
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % UV DARKHOUSE FUNCTION
        % UV is absent when ever the grain is behind saturn in the
        % co-rotating frame. Assume light from the sun travels from
        % positive y values toward smaller y values at first; allow for
        % time-dependent correction later
        alph=alph0-alph0*(abs(x)<r_sat&y<0);
        
        % Need to invent "lighthouse function"!
        % This is my "lighthouse function", see my notebook #6 pg 142-143
        %theta_trailing = -mod(2*pi/corot_period*t,2*pi);
        theta_trailing=2*pi-2*pi/corot_period*t+floor(1/corot_period*t);

        % theta_trailing refers to the trailing edge of saturn's shadow. I
        % guess another term for this might be the dawn side.
        if theta_trailing==0 || theta_trailing==pi || theta_trailing==2*pi
            % if the above statement is true, then the leading and trailing
            % edges of saturn's shadow are on the lines are located at |x|
            % = r_sat. 
            alph=alph0-alph0*(phi<theta_trailing && phi>theta_trailing ...
                && abs(x)<r_sat);
            
        else
            % figure out what the equations are for the lines representing
            % the leading and trailing edges of saturn's shadow. see my
            % notebook 36 pg 142-143 for more details.
            m = 1/(cos(theta_trailing-pi/2));   % "slope" of the shadow
            % intercept for the trailing edge
            b_I = r_sat*sin(theta_trailing)- ...
            m*r_sat*cos(theta_trailing);    
            y_I = m*x+b_I;
            % intercept for the leading edge
            b_II = r_sat*sin(theta_trailing-pi) - ...
                m*r_sat*cos(theta_trailing-pi);
            y_II = m*x+b_II;
            
            if theta_trailing > pi
                
                alph=alph0-alph0*(phi<theta_trailing && ... 
                    phi>theta_trailing && y>y_I && y<y_II);
                %disp('theta_t > pi');
            else
                alph=alph0-alph0*(phi<theta_trailing && ... 
                    phi>theta_trailing && y<y_I && y>y_II);
                %disp('theta_t < pi')
            end
            
        end
        %alph=alph0-alph0*(phi<theta_trailing && phi>theta_trailing);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        % in corotational "mode", E_x and E_y are not needed.
        % using E = B x (omega x r)
%       E_x=-x*(B*2*pi/corot_period);
% 		E_y=-y*(B*2*pi/corot_period);
        
        E_x = 0;
        E_y = 0;
        
        % allow for grain capture by saturn:
        if abs(x)<r_sat & abs(y) <r_sat
            disp('grain has been captured by saturn.')
            
            pause;
        end
    
    
    case 'enceladus_cyl' 
    %case 'enceladus'      
        % please note: this profile is set up to be used with
        % "corotating_boris_pusher.m" exclusively! Things will be different
        % if you work in a non-corotating frame.
        phi=improved_arctan(x,y);
        radius=sqrt(x^2+y^2);
        % radius of saturn, in meters, +/- 4000 m:
        r_sat=60268e3; 
		
        % space potential is set to zero, for now
		V=0;
        % assume homogeneous plasma, for now. For reference, n0~4e7 m^-3.
        % also, neutral density is 10^10 m^-3, convert this to Pascals
		ni=n0;
		ne=n0;

        % here's a stab at an electron density gradient: exponential up to
        % 10 saturn radii, but back down to a constant Temp after passing 
        % this radial distance.
        %Te = Te0*exp((radius-r_sat)/(7*r_sat))*(radius<10*r_sat)+...
        %    Te0*exp((9*r_sat)/(7*r_sat))*(radius>10*r_sat);
        Te=Te0;
        % assume no temperature gradients, for now. For reference, Te~1-10
        % eV near saturn, and Ti~10-20 eV?
        %Te=Te0;
        Ti=Ti0;
        
        % comment the line below if you don't want a log temp gradient
        %Te=log(1);
        % no flows in the co-rotating frame with saturn?
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		
		% compute the local debye length
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % check 2009? Farrell et. al. to get an estimate for B during 
        % closest approach of Cassini. Use homogeneous B-field at first,
        % but use a Dipole as simulations become more sophisticated. Keep
        % in mind also, that saturn's magnetic north pole is also its
        % geographical north pole, so if the field is a dipole the magnetic
        % field direction is pointing along the -z direction.
		B=-3.4e-7;	
        %B=-3.4e-5;
        % optional: dipole field, z component only for now:
        B=4*pi*1e-7/4/pi*(-4.5793e25/(radius.^3));
        
        % We will work in Saturn's corotating frame, so use the rotation 
        % period of saturn, in seconds:
        corot_period=10.57*3600; % hours*seconds/hour
        % gravity of saturn is not time-dependent
        G = 6.67384e-11;    % gravitational constant
        m_sat = 5.6846e26;  % mass of saturn in kg
        m_enc = 1.08e20;    % mass of enceladus in kg
        m_rhea = 2.3e21;    % mass of rhea in kg
        % distance of grain from saturn
        
        g_sat_r = -G*m_sat/radius^2;
        %g_enc_r = -G*m_enc/
        g_x=g_sat_r*cos(phi);
        g_y=g_sat_r*sin(phi);
        
        % perterbation due to enceladus:
        x_enc = 0;
        y_enc = 0;
        g_enc_r = -G*m_enc/(sqrt((x_enc-x).^2+y_enc-y).^2);
        
        % UV profile, need to calculate appropriate alph0 again, taking
        % into account distance of the grain from the sun. calculated in a
        % seemingly weird way due to how the photo-current is computed in
        % the various charge models.
        solar_distance = 9.5;   % distance from the sun in AU
        effic = 1;  % conductors have efficiencies of ~1, oxides have ~0.1
        % this uv photon flux is for solar radiation and regolith
        f_UV = 2.8e13*effic/solar_distance.^2;
        % coefficient for solar radiation and regolith; consider redo-ing 
        % for water
        alph0 = 0.25*sqrt(4*pi)*f_UV/n0/vthe;  
        % OR: uncomment the line below to turn off UV charging.
        %alph0=0;
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % UV DARKHOUSE FUNCTION
        % UV is absent when ever the grain is behind saturn in the
        % co-rotating frame. Assume light from the sun travels from
        % positive y values toward smaller y values at first; allow for
        % time-dependent correction later
        alph=alph0-alph0*(abs(x)<r_sat&y<0);
        
        % Need to invent "lighthouse function"!
        % This is my "lighthouse function", see my notebook #6 pg 142-143
        %theta_trailing = -mod(2*pi/corot_period*t,2*pi);
        theta_trailing=2*pi-2*pi/corot_period*t+floor(1/corot_period*t);

        % theta_trailing refers to the trailing edge of saturn's shadow. I
        % guess another term for this might be the dawn side.
        if theta_trailing==0 || theta_trailing==pi || theta_trailing==2*pi
            % if the above statement is true, then the leading and trailing
            % edges of saturn's shadow are on the lines are located at |x|
            % = r_sat. 
            alph=alph0-alph0*(phi<theta_trailing && phi>theta_trailing ...
                && abs(x)<r_sat);
            
        else
            % figure out what the equations are for the lines representing
            % the leading and trailing edges of saturn's shadow. see my
            % notebook 36 pg 142-143 for more details.
            m = 1/(cos(theta_trailing-pi/2));   % "slope" of the shadow
            % intercept for the trailing edge
            b_I = r_sat*sin(theta_trailing)- ...
            m*r_sat*cos(theta_trailing);    
            y_I = m*x+b_I;
            % intercept for the leading edge
            b_II = r_sat*sin(theta_trailing-pi) - ...
                m*r_sat*cos(theta_trailing-pi);
            y_II = m*x+b_II;
            
            if theta_trailing > pi
                
                alph=alph0-alph0*(phi<theta_trailing && ... 
                    phi>theta_trailing && y>y_I && y<y_II);
                %disp('theta_t > pi');
            else
                alph=alph0-alph0*(phi<theta_trailing && ... 
                    phi>theta_trailing && y<y_I && y>y_II);
                %disp('theta_t < pi')
            end
            
        end
        %alph=alph0-alph0*(phi<theta_trailing && phi>theta_trailing);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        % in corotational "mode", E_x and E_y are not needed.
        % using E = B x (omega x r)
%       E_x=-x*(B*2*pi/corot_period);
% 		E_y=-y*(B*2*pi/corot_period);
        
        E_x = 0;
        E_y = 0;
        
        % allow for grain capture by saturn:
        if abs(x)<r_sat & abs(y) <r_sat
            disp('grain has been captured by saturn.')
            
            pause;
        end
     
    case 'enceladus_temp_grad_cyl' 
        % please note: this profile is set up to be used with
        % "corotating_boris_pusher.m" exclusively! Things will be different
        % if you work in a non-corotating frame.
        phi=improved_arctan(x,y);
        radius=sqrt(x^2+y^2);
        % radius of saturn, in meters, +/- 4000 m:
        r_sat=60268e3;
        r_enc=237948000;
        r_rhea=8.76*r_sat;
        
		
        % space potential is set to zero, for now
		V=0;
        % assume homogeneous plasma, for now. For reference, n0~4e7 m^-3.
        % also, neutral density is 10^10 m^-3, convert this to Pascals
        n_scale=(r_rhea-r_enc)/log(1/20);
        n0=n0*(radius<=r_enc)+...
            n0*exp((radius-r_enc)/n_scale)*(radius<r_rhea)*(radius>r_enc)+...
            (1/20)*n0*(radius>=r_rhea);
		ni=n0;
		ne=n0;

        % here's a stab at an electron density gradient: exponential up to
        % 10 saturn radii, but back down to a constant Temp after passing 
        % this radial distance.
        
        % scale length of temperature inhomogeneity
        Re_scale=(r_rhea-r_enc)/log(10);
        Te=Te0;
        Te = Te0*(radius<=r_enc)+...
            Te0*exp((radius-r_enc)/(Re_scale))*(radius<r_rhea)*(radius>r_enc)+...
            10*Te0*(radius>=r_rhea);
        
        % assume no temperature gradients, for now. For reference, Te~1-10
        % eV near saturn, and Ti~10-20 eV?
        %Te=Te0;
        % scale length of temperature inhomogeneity
        Ri_scale=(r_rhea-r_enc)/log(100/30);
        Ti=Ti0;
        Ti=Ti0*(radius<=r_enc)+...
            Ti0*exp((radius-r_enc)/(Ri_scale))*(radius<r_rhea)*(radius>=r_enc)+...
            100/30*Ti0*(radius>=r_rhea);
        
        % comment the line below if you don't want a log temp gradient
        %Te=log(1);
        % no flows in the co-rotating frame with saturn?
		vi_x=0;vi_y=0;
        ve_x=0;ve_y=0;
        vn_x=0;vn_y=0;
		
		% compute the local debye length
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % check 2009? Farrell et. al. to get an estimate for B during 
        % closest approach of Cassini. Use homogeneous B-field at first,
        % but use a Dipole as simulations become more sophisticated. Keep
        % in mind also, that saturn's magnetic north pole is also its
        % geographical north pole, so if the field is a dipole the magnetic
        % field direction is pointing along the -z direction.
		%B=-3.4e-7;	
        %B=-3.4e-5;
        % optional: dipole field, z component only for now:
        B=4*pi*1e-7/4/pi*(-4.5793e25/(radius.^3));
        
        % We will work in Saturn's corotating frame, so use the rotation 
        % period of saturn, in seconds:
        corot_period=10.57*3600; % hours*seconds/hour
        % gravity of saturn is not time-dependent
        G = 6.67384e-11;    % gravitational constant
        m_sat = 5.6846e26;  % mass of saturn in kg
        m_enc = 1.08e20;    % mass of enceladus in kg
        m_rhea = 2.3e21;    % mass of rhea in kg
        % distance of grain from saturn
        
        g_sat_r = -G*m_sat/radius^2;
        %g_enc_r = -G*m_enc/
        g_x=g_sat_r*cos(phi);
        g_y=g_sat_r*sin(phi);
        
        % perterbation due to enceladus:
        x_enc = 0;
        y_enc = 0;
        g_enc_r = -G*m_enc/(sqrt((x_enc-x).^2+y_enc-y).^2);
        
        % UV profile, need to calculate appropriate alph0 again, taking
        % into account distance of the grain from the sun. calculated in a
        % seemingly weird way due to how the photo-current is computed in
        % the various charge models.
        solar_distance = 9.5;   % distance from the sun in AU
        effic = 1;  % conductors have efficiencies of ~1, oxides have ~0.1
        % this uv photon flux is for solar radiation and regolith
        f_UV = 2.8e13*effic/solar_distance.^2;
        % coefficient for solar radiation and regolith; consider redo-ing 
        % for water
        alph0 = 0.25*sqrt(4*pi)*f_UV/n0/vthe;  
        % OR: uncomment the line below to turn off UV charging.
        %alph0=0;
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % UV DARKHOUSE FUNCTION
        % UV is absent when ever the grain is behind saturn in the
        % co-rotating frame. Assume light from the sun travels from
        % positive y values toward smaller y values at first; allow for
        % time-dependent correction later
        alph=alph0-alph0*(abs(x)<r_sat&y<0);
        
        % Need to invent "lighthouse function"!
        % This is my "lighthouse function", see my notebook #6 pg 142-143
        %theta_trailing = -mod(2*pi/corot_period*t,2*pi);
        theta_trailing=2*pi-2*pi/corot_period*t+floor(1/corot_period*t);

        % theta_trailing refers to the trailing edge of saturn's shadow. I
        % guess another term for this might be the dawn side.
        if theta_trailing==0 || theta_trailing==pi || theta_trailing==2*pi
            % if the above statement is true, then the leading and trailing
            % edges of saturn's shadow are on the lines are located at |x|
            % = r_sat. 
            alph=alph0-alph0*(phi<theta_trailing && phi>theta_trailing ...
                && abs(x)<r_sat);
            
        else
            % figure out what the equations are for the lines representing
            % the leading and trailing edges of saturn's shadow. see my
            % notebook 36 pg 142-143 for more details.
            m = 1/(cos(theta_trailing-pi/2));   % "slope" of the shadow
            % intercept for the trailing edge
            b_I = r_sat*sin(theta_trailing)- ...
            m*r_sat*cos(theta_trailing);    
            y_I = m*x+b_I;
            % intercept for the leading edge
            b_II = r_sat*sin(theta_trailing-pi) - ...
                m*r_sat*cos(theta_trailing-pi);
            y_II = m*x+b_II;
            
            if theta_trailing > pi
                
                alph=alph0-alph0*(phi<theta_trailing && ... 
                    phi>theta_trailing && y>y_I && y<y_II);
                %disp('theta_t > pi');
            else
                alph=alph0-alph0*(phi<theta_trailing && ... 
                    phi>theta_trailing && y<y_I && y>y_II);
                %disp('theta_t < pi')
            end
            
        end
        %alph=alph0-alph0*(phi<theta_trailing && phi>theta_trailing);
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % in corotational "mode", E_x and E_y are not needed.
        % using E = B x (omega x r)
%       E_x=-x*(B*2*pi/corot_period);
% 		E_y=-y*(B*2*pi/corot_period);
        
        % don't need electric field in the co-rotating frame!
        E_x = 0;
        E_y = 0;
        
        % allow for grain capture by saturn:
        if abs(x)<r_sat & abs(y) <r_sat
            disp('grain has been captured by saturn.')
            
            pause;
        end
        
	case 'grad_B'
        % magnetic field on axis. For Te=1.6 eV, strength of the gradient 
        % in T/m; maximum gradient is 2 T/m in Auburn machine
		B0=5;	
		beta_B=2;	
		
        ni=n0;
        ne=n0;
        V=0;
        E_x=0;
        E_y=0;
        Ti=Ti0;Te=Te0;
        % IN A SHEATH PROFILE, USE ve_x TO REPRESENT vi_z!!!
        % Also, use ve_y to represent ve_y.
   
        % assume the bohm speed for now, check this later.
        ve_x=sqrt(qe*Te/mi);
       	% Need to account for diamagnetic electron/ion currents
        vi_x=0;vi_y=0;
        ve_y=0;
        vn_x=0;vn_y=0;
        
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
		V=0;
		ni=n0;
		ne=n0;
		alph=0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
		xmax=beta_B/B0;
		% Assume a simple, linear decrease in magnetic field with radial 
		% distance from x=0. Make B=B0 for x<0.
		B=(B0-beta_B*x)*(x<xmax)*(x>=0)+B0*(x<0);	
        % corot_period is unused for this profile.
        corot_period=0;
	case 'grad_B_cyl'
        % magnetic field on axis. For Te=1.6 eV, strength of the gradient 
        % in T/m; maximum gradient is 2 T/m in Auburn machine
		B0=5;	
		beta_B=2;	
		phi=improved_arctan(x,y);
		% radius of the experimental volume, in meters. This corresponds 
        % exactly with the "central experimental volume"
		R=0.225;	 
        % on the Auburn machine; the smaller, "uniform region" radius is 
        % 0.1
		
        % assume that grains levitate near the sheath (not entirely
        % realistic, but better than nothing.)
        n0=n0*exp(-1/2);
        ni=n0;
        ne=n0;
        radius=sqrt(x^2+y^2);
        V=0;
        E_x=0;
        E_y=0;
        Ti=Ti0;Te=Te0;
        vi_x=0;vi_y=0;
    	% IN A SHEATH PROFILE, USE ve_x TO REPRESENT vi_z!!!
        % assume the bohm speed for now if we are at the sheath edge; 
        % check this later. negative because the ions flow in the opposite
        % direction of Z (toward the planar electrode.)
        ve_x=-sqrt(qe*Te/mi);
        
        %ve_x=0;
      	% Need to account for diamagnetic electron/ion currents
        ve_y=0;
        vn_x=0;vn_y=0;
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
		V=0;
		ni=n0;
		ne=n0;
		alph=0;
		lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % NEED TO CHECK THIS LINE BELOW! I AM BASING IT OFF OF
        % 1992_DAUGHERTY_JAP. THE IDEA HERE IS THAT NEAR-MONO-ENERGETIC
        % IONS THAT FALL OUT OF THE SHEATH HAVE AN ENERGY COMPARABLE TO Te,
        % SO THIS ENLARGES THE SIZE OF THE SHEATH.
        %lambda_D=sqrt((eps0*Te*Te)/qe/(ne*Te+ni*Te));
        
		% Assume a simple, linear decrease in magnetic field with radial 
        % distance from r=0.
		B=B0-beta_B*radius;
        % corot_period is unused for this profile.
        corot_period=0;
    case 'auburn_cyl'
    %case 'auburn'
        % magnetic field on axis. For Te=1.6 eV, strength of the gradient 
        % in T/m; maximum gradient is 2 T/m in Auburn machine
		%B0=5;
        B0=4;
		%beta_B=2;	
		phi=improved_arctan(x,y);
		% radius of the experimental volume, in meters. This corresponds 
        % exactly with the "central experimental volume"
		R=0.225;	 
        % on the Auburn machine; the smaller, "uniform region" radius is 
        % 0.1
		
        
        radius=sqrt(x^2+y^2);
        % temperature profiles: assumed constant
        Ti=Ti0;Te=Te0;
        ni=exp(-1/2)*n0*sech(radius/(1*R));
        ne=exp(-1/2)*n0*sech(radius/(1*R));
        V=0;
        % radial electric field in V/m:
        Er=-100;
        Er=-1;
        E_x=Er*cos(phi);
        E_y=Er*sin(phi);
        % If an electric field is present, then ion flow is present
        vi_x=Er*sin(phi)/B0;
        vi_y=-Er*cos(phi)/B0;
        vi_x=0;
        vi_y=0;
        % IN A SHEATH PROFILE, USE ve_x TO REPRESENT vi_z!!!
        ve_x=sqrt(qe*Te/mi);
        ve_y=-Er/B0;
        ve_x=0;
        ve_y=0;
        vn_x=0;
        vn_y=0;
        
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
		alph=0;
        
		%lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % NEED TO CHECK THIS LINE BELOW! I AM BASING IT OFF OF
        % 1992_DAUGHERTY_JAP. THE IDEA HERE IS THAT NEAR-MONO-ENERGETIC
        % IONS THAT FALL OUT OF THE SHEATH HAVE AN ENERGY COMPARABLE TO Te,
        % SO THIS ENLARGES THE SIZE OF THE SHEATH.
        lambda_D=sqrt((eps0*Te*Te)/qe/(ne*Te+ni*Te));
		
        % Assume a simple, linear decrease in magnetic field with radial 
        % distance from r=0.
		%B=B0-beta_B*radius;
        B=B0;
        % corot_period is unused for this profile.
        corot_period=0;
        
        case 'auburn_ring_cyl'
        % magnetic field on axis. For Te=1.6 eV, strength of the gradient 
        % in T/m; maximum gradient is 2 T/m in Auburn machine
		%B0=5;
        B0=4;
		%beta_B=2;	
		phi=improved_arctan(x,y);
		% radius of the experimental volume, in meters. This corresponds 
        % exactly with the "central experimental volume"
		R=0.225;	 
        % on the Auburn machine; the smaller, "uniform region" radius is 
        % 0.1
		
        
        radius=sqrt(x^2+y^2);
        % temperature profiles: assumed constant
        Ti=Ti0;Te=Te0;
        
   
        % define radial electric field stuff first
        E0=500;
        % assume an off-center gaussian, specified by a "sharpness" R0 and 
        % by a location parameter r0
        r0=0.05;
        %R0=0.02;
        FWHM=0.04;
        R0=FWHM/2/sqrt(log(2));
        Er=E0*exp(-(radius-r0).^2/R0.^2);
        E_x=Er*cos(phi);
        E_y=Er*sin(phi);
        % compute space potential, assuming V0=V(r=0) is known
        V=-E0*R0*(erf(r0/R0)+erf((radius-r0)/R0));
        %V=-E0*R0*erfc((radius-r0/R0));
        % density profiles, assume boltzmann electrons?:
        ns=exp(-1/2)*n0;
        ne=ns*exp(V/Te);
        % I don't think the next line is correct
        %ni=-eps0/qe*(2*(radius-r0)/R.^2.*Er)+ne;
        ni=ne+(eps0/qe)*(Er./radius-2/R0.^2*(radius-r0).*Er);
        % If an electric field is present, then ion flow is present. The
        % correct convention is that if the electric field is in the
        % positive direction, then ion flow is in the phi direction.
        vi_x=Er*sin(phi)/B0;
        vi_y=-Er*cos(phi)/B0;
        % IN A SHEATH PROFILE, USE ve_x TO REPRESENT vi_z!!!
        ve_x=sqrt(qe*Te/mi);
        %ve_x=0;
        % use ve_y for the phi direction!
        ve_y=-Er/B0;
        vn_x=0;
        vn_y=0;
        
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
		alph=0;
        
		%lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % NEED TO CHECK THIS LINE BELOW! I AM BASING IT OFF OF
        % 1992_DAUGHERTY_JAP. THE IDEA HERE IS THAT NEAR-MONO-ENERGETIC
        % IONS THAT FALL OUT OF THE SHEATH HAVE AN ENERGY COMPARABLE TO Te,
        % SO THIS ENLARGES THE SIZE OF THE SHEATH.
        lambda_D=sqrt((eps0*Te*Te)/qe/(ne*Te+ni*Te));
		
        % Assume a simple, linear decrease in magnetic field with radial 
        % distance from r=0.
		%B=B0-beta_B*radius;
        B=B0;
        % corot_period is unused for this profile.
        corot_period=0;
        
        
        
        
    case 'auburn_cyl_parabolic'
        % magnetic field on axis. For Te=1.6 eV, strength of the gradient 
        % in T/m; maximum gradient is 2 T/m in Auburn machine
		%B0=5;
        B0=4;
		%beta_B=2;	
		phi=improved_arctan(x,y);
		% radius of the experimental volume, in meters. This corresponds 
        % exactly with the "central experimental volume"
		R=0.225;	 
        % on the Auburn machine; the smaller, "uniform region" radius is 
        % 0.1
		
        
        radius=sqrt(x^2+y^2);
        % temperature profiles: assumed constant
        Ti=Ti0;Te=Te0;
        
        % compute space potential, assuming V0=V(r=0) is known:
        V0=25;  % keep in mind this is a "negative" parabolic potential 
        r0=0.1; % where the potential stops, and the plasma potential is 
                % just flat
        % density profiles, assume boltzmann electrons?:
        ns=exp(-1/2)*n0;
        if radius<=r0;
           	V=-abs(V0)+abs(V0)*radius.^2/r0.^2;
            Er=-2*abs(V0)*radius/r0.^2;
            ne=ns*exp(V/Te);
            ni=ne-(eps0/qe)*(4*abs(V0)/r0.^2);
        else
            V=0;
            Er=0;
            ne=ns*exp(V/Te);
            ni=ne-(eps0/qe)*(4*abs(V0)/r0.^2);
        end
        E_x=Er*cos(phi);
        E_y=Er*sin(phi);
        
                
        
        % If an electric field is present, then ion flow is present. The
        % correct convention is that if the electric field is in the
        % positive direction, then ion flow is in the phi direction.
        vi_x=Er*sin(phi)/B0;
        vi_y=-Er*cos(phi)/B0;
        %vi_x=0;
        %vi_y=0;
        % IN A SHEATH PROFILE, USE ve_x TO REPRESENT vi_z!!!
        ve_x=sqrt(qe*Te/mi);
        %ve_x=0;
        % use ve_y for the phi direction!
        ve_y=-Er/B0;
        %ve_y=0;
        vn_x=0;
        vn_y=0;
        
        % % gravity unimportant since it is not in the plane or too weak.
        g_x=0;g_y=0;
		alph=0;
        
		%lambda_D=sqrt((eps0*Ti*Te)/qe/(ne*Ti+ni*Te));
        % NEED TO CHECK THIS LINE BELOW! I AM BASING IT OFF OF
        % 1992_DAUGHERTY_JAP. THE IDEA HERE IS THAT NEAR-MONO-ENERGETIC
        % IONS THAT FALL OUT OF THE SHEATH HAVE AN ENERGY COMPARABLE TO Te,
        % SO THIS ENLARGES THE SIZE OF THE SHEATH.
        lambda_D=sqrt((eps0*Te*Te)/qe/(ne*Te+ni*Te));
		
        % Assume a simple, linear decrease in magnetic field with radial 
        % distance from r=0.
		%B=B0-beta_B*radius;
        B=B0;
        % corot_period is unused for this profile.
        corot_period=0;
% use endswitch with octave, end with matlab.
%endswitch
end

% the more correct description for the grain capacitance:
% C=4*pi*eps0*a*exp(-a/lambda_D)
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




