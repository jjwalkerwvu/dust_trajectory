% oml_monoenergetic_ions.m

function [Itot,qf,Kn_R0,P0,P1,Pg1]=...
    oml_monoenergetic_ions(qflag,a,alph,Ti,n_e,n_i,Te,B,Z,C,qd,lambda_D,...
    lambda_i,w,species);
%   explanation of inputs:
%	  qflag 	= whether or not to evaluate equilibrium charge for a given 
%		          model; qflag=1 means do calculate q_eq, qflag=0 means do not 
%             calculate q_eq.
%   a       = dust grain size, in meters
%   alph    = coefficient of UV illumination, something like f_uv/ne/vthe
%   Ti      = ion temperature, in eV
%   ne      = electron density in m^-3
%   ni      = ion density in m^-3
%   Te      = electron temperature, in eV
%   B       = Magnetic field strength, in Tesla
%   Z       = charge state of ions; use 1 for singly ionized plasma
%   C       = Capacitance of dust grain in Farads
%   qd      = grain charge, in Coloumbs
%   lambda_D    = linearized debye length, in meters
%   lambda_i    = mean free path of ions, in meters
%   w       = w is a two element array, given by [we wi], where we is the 
%             grain speed relative to the electrons, and wi is the grain
%             speed relative to the ions; both are in units of m/s
%   species     = this is the mass number of the ion species.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Explanation of outputs:
% Itot    = total current to the dust grain, in amps
% qf      = equilibrium charge on the dust grain, in coloumbs
% Kn_R0   = Knudsen capture radius of grain, dimensionless
% P0      = probability ion has no collisions in grain sheath
% P1      = probability ion has 1 collision in grain sheath
% Pg1     = probability ion has many collisions in grain sheath
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

eta=n_e/n_i;
Tau=Te/Ti;
mr=me/mi;
vthe=sqrt(2*qe*Te/me);      % local electron thermal speed, m/s
vthi=sqrt(mr/Tau)*vthe;		% local ion (proton) thermal speed, m/s

% electron thermal mach number:
Me=w(1)/vthe;
% ion thermal mach number:
Mi=w(2)/vthi;


% % magnetization parameters
if B==0
	% % error checking, in case the magnetic field is identically zero.
	e_mag=0;
	i_mag=0;
else
    % to match up with Patacchini and Hutchinson 2007
	Rle=sqrt(pi/4)*me*vthe/qe/B;	
    % to match up with Patacchini and Hutchinson 2007
	Rli=sqrt(pi/4)*mi*vthi/Z/qe/B;	
	
	e_mag=a/Rle;
	i_mag=a/Rli;
	% % The size of the dust grain and debye sheath combined
	a_and_s=a+2.5*lambda_D;
	% % The following may be more appropriate for magnetization ratios:
	%e_mag=a_and_s/Rle;
	%i_mag=a_and_s/Rli;
end



% % the following parameters are unneeded for the OML model, so they are 
% % just set to zero for when they are needed in the call to 
% % dust_bisection.m
%lambda_i=inf;
mu_i=0;


% % Knudsen capture radius 
Kn_R0=0;
P0=1;
P1=0;
Pg1=0;


% Calculate plasma currents to particle
if qd<=0	% negative dust potential (phi=qd/C)

	if e_mag<1	%% UNMAGNETIZED ELECTRONS	
        if Me==0
        % derived via integration of Maxwellian from vmin=sqrt(2*qe*phi/me) 
        % to infinity, given 4*pi*a^2 collection area
        %Ie=-ne*qe*(4*pi*a^2)*vthe*exp(qd/C/Te);
            Ie=-n_e*qe*(4*pi*a^2)/sqrt(4*pi)*vthe*exp(qd/C/Te);
        	% expression above matches Patacchini and Hutchinson, 2007
        else
        	%Ie=-n_e*qe*(4*pi*a^2)/sqrt(4*pi)*vthe*exp(qd/C/Te);
        	Ie=-.5*sqrt(pi)*a.^2*n_e*qe*vthe/Me*(...
                (Me.^2+.5+qd/C/Te)*sqrt(pi)*(erf(Me+sqrt(-qd/C/Te))+...
                erf(Me-sqrt(-qd/C/Te)))+...
                (sqrt(-qd/C/Te)+Me)*exp(-(Me-sqrt(-qd/C/Te)).^2)-...
                (sqrt(-qd/C/Te)-Me)*exp(-(Me+sqrt(-qd/C/Te)).^2));
        end
        		
	else	%% MAGNETIZED ELECTRONS; this 		
        % % same as unmagnetized case, except collection area reduced due 
        % % to magnetization; may be off by some constant factor due to 
        % % cos-dependence of incidence angle
       	if Me==0
        % derived via integration of Maxwellian from vmin=sqrt(2*qe*phi/me) 
        % to infinity, given 4*pi*a^2 collection area
        %Ie=-ne*qe*(4*pi*a^2)*vthe*exp(qd/C/Te);
            Ie=-.5*n_e*qe*(4*pi*a^2)/sqrt(4*pi)*vthe*exp(qd/C/Te);
        	% expression above matches Patacchini and Hutchinson, 2007
        else
        	%Ie=-n_e*qe*(4*pi*a^2)/sqrt(4*pi)*vthe*exp(qd/C/Te);
        	Ie=-.25*sqrt(pi)*a.^2*n_e*qe*vthe/Me*(...
                (Me.^2+.5+qd/C/Te)*sqrt(pi)*(erf(Me+sqrt(-qd/C/Te))+...
                erf(Me-sqrt(-qd/C/Te)))+...
                (sqrt(-qd/C/Te)+Me)*exp(-(Me-sqrt(-qd/C/Te)).^2)-...
                (sqrt(-qd/C/Te)-Me)*exp(-(Me+sqrt(-qd/C/Te)).^2));
        end
       		
	end
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if i_mag<1	%% UNMAGNETIZED IONS
		
        	% % OML current, cf. Allen, Phys. Scr. 45 (1992), eq. 51	
        		%Ii=ni*qe*(4*pi*a^2)*sqrt(-2*qe*qd/C/mi);
        	% % I'm replacing the above line with the current I see more  
        	% % often in the literature (Shukla 2001 pop, for example)
        	% Include a statement for w==0 so that there is no division by 
            % zero.
		if Mi==0
            % this is OML theory for mono-energetic ions, so if there is no
            % ion flow, there can be no ion current.
            Ii=0;
			%Ii=n_i*qe*(4*pi*a^2)/sqrt(4*pi)*vthi*(1-qd/C/Ti);
			% expression above matches Patacchini and Hutchinson, 2007
        else
            Ii=pi*n_i*qe*a^2*w(2)*(1-qe*qd/C/(mi*w(2).^2));
% 			Ii=sqrt(pi)*a.^2*n_i*Z*qe*vthi*(...
%                 (Mi.^2+0.5-qd/C/Ti)*sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));			
			% % see 1992_northrop_ps or 1981_whipple_repprogphys for the 
			% % above. Also, 1996 Horanyi and 1996 Northrop.
		end
	
	% %  Call the Newton-Raphson method to find the equilibrium charge FOR 
    % % THE INPUT CONDITIONS.
        if qflag==1
			% % a charge between 0 and 1e6 elementary charges is a good 
            % % search interval.
			[output]=dust_bisection(1e6,'oml_monoenergetic_ions',(qd/qe),eta,alph,Ti,Te,...
                e_mag,i_mag,C,lambda_D,lambda_i,a,w,species);
			qf=round(output);
			qf=qf*qe;
		else
			% % return 0 if you do not wish to calculate equilibrium charge 
            % % (qflag=0.)
			qf=0;
            % % just return zero, value does not matter we just want to 
            % % skip the charge calculation step.
        end
        
        % now that charge, potential profile, and ion mean free path are 
        % known, the Capture Radius can be calculated at the mean ion 
        % thermal kinetic energy.
        R0 = (abs(qd/C)*a*(1+a/lambda_D))/(1.5*Ti+abs(qd/C)*a/lambda_D);
        if R0==0 || lambda_i==inf
            % If R0=0, use OML currents??
            P0=1;
            P1=0;
            Pg1=0;
            Kn_R0=0;    				
        else
            % compute the Knudsen Capture radius; see 2008 Gatti PRE for 
            % details
            Kn_R0 = lambda_i/(2*1.22*R0);
            P0=exp(-1/(Kn_R0));
        	P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
            Pg1=1-(P0+P1);
        end
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % % MAGNETIZED IONS
	else					    
       	% % simple thermal flux of ions, assuming ballistic trajectories  
        % % along field lines w/ reduced collection area due to 
        % % magnetization
        Ii=n_i*qe*(2*pi*a^2)*vthi;	
        % reference for the above?? Maybe 1982 Rubenstein phys. fluids
        
        % % for magnetized ions and electrons, a simple charge model is 
        % % sufficient.
		if qflag==1
			
			qf=0.5*C*Te*log(mr/Tau/eta/eta);
			qf=round(qf/qe);
			qf=qf*qe;
		else
			% return 0 if you do not wich to calculate equilibrium charge 
            % (qflag=0.)
			qf=0;	
            % just return zero, value does not matter we just want to skip
			% the charge calculation step.
        end
        % I think Kn_R0 is just zero if there the ions are magnetized?
        Kn_R0=0;
	end
	ph_flux=alph*4*n_i*eta*vthe;	%% just set ph_flux=4*n0*vthe for now!
	% FOR ENCELADUS, OR OTHER SOLAR SYSTEM RELEVANT SITUATIONS 
	% (Horanyi 1996): ph_flux=2.5e10 K/d/d; this has units cm^-2 s^-1
	% where d is the distance from the sun in AU
	% and K is the efficiency factor; ~1 for conductors and ~0.1 for 
    % dielectrics.
	Ip=qe*pi*a^2*ph_flux;		%% only valid for q<0! (Horanyi 1996)
    % % Not sure why I have to multiply by sqrt(2*pi)... (august 2013)
    Ip=1/sqrt(4*pi)*Ip;
    % better expression: (4/24/2014)
    Ip=qe*n_i*eta*vthe*alph*pi*a^2;
    
	Itot=Ii+Ie+Ip;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% can add equations for positive dust potential here; I think the below is 
% valid only for unmagnetized electrons, ions
else
	% Whipple 1981, reviews of geophysics; 1995 Cui and Goree IEEE
	Ie = -n_e*qe*(4*pi*a^2)*vthe/sqrt(4*pi)*(1+qd/C/Te);	
	Ii = n_i*qe*(4*pi*a^2)*vthi/sqrt(4*pi)*exp(-qd/C/Ti);
	% Need to call the function for finding equilibrium charge:
	% Call the Newton-Raphson method to find the equilibrium charge FOR THE 
    % INPUT CONDITIONS.
	if qflag==1
		% a charge between 0 and 1e6 elementary charges is a good search 
        % interval.
		[output]= dust_bisection(1e6,'oml_monoenergetic_ions',(qd/qe),eta,alph,Ti,Te,...
            e_mag,i_mag,C,lambda_D,lambda_i,mu_i,a,w,species);
		qf=round(output);
		qf=qf*qe;
	else
		% return 0 if you do not wich to calculate equilibrium charge 
        % (qflag=0.)
		qf=0;	
        % just return zero, value does not matter we just want to skip
        % the charge calculation step.
	end
%%%%%%	%%	FOLLOWING NEEDS TO BE CHANGED FOR POSITIVE DUST POTENTIAL!!	
	ph_flux=alph*4*n_i*eta*vthe;	%% just set ph_flux=4*n0*vthe for now!
	% FOR ENCELADUS, OR OTHER SOLAR SYSTEM RELEVANT SITUATIONS 
	% (Horanyi 1996): ph_flux=2.5e10 K/d/d; this has units cm^-2 s^-1
	% where d is the distance from the sun in AU
	% and K is the efficiency factor; ~1 for conductors and ~0.1 for 
    % dielectrics.
	Ip=qe*pi*a^2*ph_flux;		% only valid for q<0! (Horanyi 1996)
    % better expression: (4/24/2014)
    Ip=qe*n_i*eta*vthe*alph*pi*a^2;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	%Ip=qe*pi*a^2*ph_flux;		%% only valid for q<0! (Horanyi 1996)
	%Qabs=1;		% absorption efficiency is ~1 for 2*pi*a/lambda_uv>1.
    % the yield of photoelectrons; perhaps use Y=1? so we get 1 
    % photoelectron for every uv photon?
	%Y=1;		
	%juv=alph*1e15;	%% UV photon flux, in units of m^-2 s^-1?
    % Thermal energy of photo electrons; for now I'm assuming Tph=1 eV for 
    % convenience.
	%Tph=1;	
    % % only valid for q>0! (Shukla 2001)
	%Ip=qe*pi*a^2*Qabs*Y*juv*exp(-qe*qd/C/Tph);	
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	Itot=Ii+Ie+Ip;
    % Knudsen capture radius does not exist.
    Kn_R0=0;
end