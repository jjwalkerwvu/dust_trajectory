function [Itot,qf,Kn_R0,P0,P1,Pg1]=hutchinson(qflag,a,alph,species,...
    Ti,n_e,n_i,Te,B,Z,C,qd,lambda_D,lambda_i,w);

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
qe=1.6e-19;
me=9.1e-31;
mp=1.67e-27;
eps0=8.854e-12;
mi=species*mp;

eta=n_e/n_i;
Tau=Te/Ti;
mr=me/mi;
vthe=sqrt(2*qe*Te/me);	% local electron thermal speed, m/s
vthi=sqrt(mr/Tau)*vthe;	% local ion (proton) thermal speed, m/s

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
	%Rle=me*vthe/qe/B;
	Rle=sqrt(pi/4)*me*vthe/qe/B;	%% to match up with Patacchini and Hutchinson 2007
	% % I think to match up with Patacchini and Hutchinson 2007, just use exactly their formula:
	% e_mag=a/sqrt(pi*Te*me/2/qe/B/B);	% Sept 2013
	
	%Rli=mi*vthi/Z/qe/B;
	Rli=sqrt(pi/4)*mi*vthi/Z/qe/B;	%% to match up with Patacchini and Hutchinson 2007
	% % I think to match up with Patacchini and Hutchinson 2007, just use exactly their formula:
	% i_mag=a/sqrt(pi*Ti*mi/2/qe/B/B);		% Sept 2013
	
	e_mag=a/Rle;
	i_mag=a/Rli;
	% % The size of the dust grain and debye sheath combined
	a_and_s=a+2.5*lambda_D;
	% % The following may be more appropriate for magnetization ratios:
	%e_mag=a_and_s/Rle;
	%i_mag=a_and_s/Rli;
end
%% the following parameters are not need yet, so they are just set to zero for when they
%% are needed in the call to dust_bisection.m
%lambda_i=inf;
mu_i=0;
% % Knudsen capture radius 
Kn_R0=0;
P0=1;
P1=0;
Pg1=0;

% % FOR LAMBDA_D --> INFINITE???
% z is the function of magnetization used by Patacchini and Hutchinson, 2007
%z=e_mag/(1+e_mag);

% iota* in Patacchini and Hutchinson:
%iota=1-0.0946*z-0.305*z.^2+0.95*z.^3-2.2*z.^4+1.15*z.^5;

%% lower bound on electron current:
%Ie_low=4*pi*a.^2*vthe/2/sqrt(pi)*iota*exp(qd/C/Te);
%
%% eta, which is (Vp/Te)/e_mag:
%eta_mag=-qd/C/Te/e_mag;
%
%% w, which is eta/(1+eta):
%w_mag=eta_mag/(1+eta_mag);
%
%% A, the fitting polynomial, a function of w:
%A_fit=0.678*w_mag+1.543*w_mag.^2-1.212*w_mag.^3;
%
%% Ie*, which is the empirical formula for electron current as a function of magnetization
%Ie=4*pi*a.^2*vthe/2/sqrt(pi)*exp(qd/C/Te)*(A_fit+(1-A_fit)*iota);

% % FOR LAMBDA_D = FINITE, AND DEBYE-HUCKEL POTENTIAL:
% eta, which is now dependent on grain sheath size:
%eta_mag=-qd/C/Te/e_mag*(1+a/lambda_D);
%eta_mag=-qd/C/Te/e_mag*(1+e_mag/4*(1-exp(-4*a/lambda_D/e_mag)));

% w, which is eta/(1+eta):
%w_mag=eta_mag/(1+eta_mag);

% A, the fitting polynomial, a function of w:
%A_fit=0.678*w_mag+1.543*w_mag.^2-1.212*w_mag.^3;

% Ie*, which is the empirical formula for electron current as a function of magnetization
%Ie=4*pi*a.^2*vthe/2/sqrt(pi)*exp(qd/C/Te)*(A_fit+(1-A_fit)*iota);

% Patacchini-Hutchinson model can include the spatial dependence of electron flux to the sphere.
% I have chosen not to put this in at the current time, since it is not necessary for determining the
% total grain charge.

% % Check magnetization of ions, ion current is the same as OML.
% % Comment out the electron currents from OML model
if qd<=0	% negative dust potential (phi=qd/C)
% % Don't have to split electrons into magnetized and unmagnetized regimes; 
% % this is taken care of in one expression for the electron current in the 
% % Patacchini-Hutchinson charge model.
	z=e_mag/(1+e_mag);

	% iota* in Patacchini and Hutchinson:
	iota=1-0.0946*z-0.305*z.^2+0.95*z.^3-2.2*z.^4+1.15*z.^5;
	% % FOR LAMBDA_D = FINITE, AND DEBYE-HUCKEL POTENTIAL:
	% eta, which is now dependent on grain sheath size:
	eta_mag=-qd/C/Te/e_mag*(1+e_mag/4*(1-exp(-4*a/lambda_D/e_mag)));

	% w, which is eta/(1+eta):
	w_mag=eta_mag/(1+eta_mag);

	% A, the fitting polynomial, a function of w:
	A_fit=0.678*w_mag+1.543*w_mag.^2-1.212*w_mag.^3;

	% Ie*, which is the empirical formula for electron current as a 
    % function of magnetization
	Ie=-qe*4*pi*a.^2*vthe/2/sqrt(pi)*n_e*exp(qd/C/Te)*...
        (A_fit+(1-A_fit)*iota);

% Patacchini-Hutchinson model can include the spatial dependence of 
% electron flux to the sphere. I have chosen not to put this in at the 
% current time, since it is not necessary for determining the total grain 
% charge.
	
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if i_mag<1	%% UNMAGNETIZED IONS
        
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
		
        % ion current, for all thermal mach numbers.
		if Mi==0
			Ii=n_i*qe*(4*pi*a^2)/sqrt(4*pi)*vthi*(1-qd/C/Ti);
			% expression above matches Patacchini and Hutchinson, 2007
		else
			Ii=sqrt(pi)*a.^2*n_i*Z*qe*vthi*((Mi.^2+0.5-qd/C/Ti)*...
                sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));			
			% % see 1992_northrop_ps or 1981_whipple_repprogphys for the 
			% % above. Also, 1996 Horanyi and 1996 Northrop.
		end
	
	% %  Call the Newton-Raphson method to find the equilibrium charge FOR 
    % % THE INPUT CONDITIONS.
		if qflag==1
			% % a charge between 0 and 1e4 elementary charges is a good 
            % % search interval.
			[output]= dust_bisection(1e6,'hutchinson',(qd/qe),eta,alph,...
                Ti,Te,e_mag,i_mag,C,lambda_D,lambda_i,a,w,species);
			qf=round(output);
			qf=qf*qe;
		else
			% % return 0 if you do not wish to calculate equilibrium charge 
            % % (qflag=0.)
			qf=0;	% % just return zero, value does not matter we just  
                    % % want to skip the charge calculation step.
        end
        
        
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
	else	%% MAGNETIZED IONS				    
       	% % simple thermal flux of ions, assuming ballistic trajectories 
    	% % along  field lines w/ reduced collection area due to magnetization
    	Ii=n_i*qe*(2*pi*a^2)*vthi;	% reference for this?? Maybe 1982 Rubenstein phys. fluids
	%% for magnetized ions and electrons, a simple charge model is sufficient.
		if qflag==1
			
			qf=0.5*C*Te*log(mr/Tau/eta/eta);
			qf=round(qf/qe);
			qf=qf*qe;
		else
			%% return 0 if you do not wich to calculate equilibrium charge (qflag=0.)
			qf=0;	%% just return zero, value does not matter we just want to skip
				%% the charge calculation step.
		end
	end
	ph_flux=alph*4*n_i*eta*vthe;	%% just set ph_flux=4*n0*vthe for now!
	%% FOR ENCELADUS, OR OTHER SOLAR SYSTEM RELEVANT SITUATIONS (Horanyi 1996):
	%% ph_flux=2.5e10 K/d/d; this has units cm^-2 s^-1
	%% where d is the distance from the sun in AU
	%% and K is the efficiency factor; ~1 for conductors and ~0.1 for dielectrics.
	Ip=qe*pi*a^2*ph_flux;		%% only valid for q<0! (Horanyi 1996)
   	% % Not sure why I have to multiply by sqrt(2*pi)... (august 2013)
  	Ip=1/sqrt(4*pi)*Ip;
    % better expression: (4/24/2014)
    Ip=qe*n_i*eta*vthe*alph*pi*a^2;

	Itot=Ii+Ie+Ip;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else
	% can add equations for positive dust potential here; I think the below is valid only for unmagnetized
	% electrons, ions
	Ie = -n_e*qe*(4*pi*a^2)*vthe/sqrt(4*pi)*(1+qd/C/Te);	%% Whipple 1981, reviews of geophysics; 1995 Cui and Goree IEEE
	Ii = n_i*qe*(4*pi*a^2)*vthi/sqrt(4*pi)*exp(-qd/C/Ti);
	%% Need to call the function for finding equilibrium charge:
	%%  Call the Newton-Raphson method to find the equilibrium charge FOR THE INPUT CONDITIONS.
	if qflag==1
		%% a charge between 0 and 1e4 elementary charges is a good search interval.
		[output]= dust_bisection(1e6,'oml',(qd/qe),eta,alph,Ti,Te,e_mag,i_mag,C,lambda_D,lambda_i,mu_i,a,w,species);
		%[output]=dust_newt_meth(1e4,'oml',(qd/qe),eta,alph,Ti,Te,e_mag,i_mag,C);
		qf=round(output);
		qf=qf*qe;
	else
		%% return 0 if you do not wich to calculate equilibrium charge (qflag=0.)
		qf=0;	%% just return zero, value does not matter we just want to skip
			%% the charge calculation step.
	end
%%%%%%	%%	FOLLOWING NEEDS TO BE CHANGED FOR POSITIVE DUST POTENTIAL!!	
	ph_flux=alph*4*n_i*eta*vthe;	%% just set ph_flux=4*n0*vthe for now!
	% % FOR ENCELADUS, OR OTHER SOLAR SYSTEM RELEVANT SITUATIONS 
	% % (Horanyi 1996): ph_flux=2.5e10 K/d/d; this has units cm^-2 s^-1
	% % where d is the distance from the sun in AU
	% % and K is the efficiency factor; ~1 for conductors and ~0.1 for 
    % % dielectrics.
	Ip=qe*pi*a^2*ph_flux;		% only valid for q<0! (Horanyi 1996)
    % % Not sure why I have to multiply by sqrt(2*pi)... (august 2013)
    %Ip=1/sqrt(4*pi)*Ip;
    % better expression: (4/24/2014)
    Ip=qe*n_i*eta*vthe*alph*pi*a^2;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	%Ip=qe*pi*a^2*ph_flux;		%% only valid for q<0! (Horanyi 1996)
	%Qabs=1;		%% absorption efficiency is ~1 for 2*pi*a/lambda_uv > 1.
	%Y=1;		%% the yield of photoelectrons; perhaps use Y=1? so we get 1 photoelectron for every uv 		photon?
	%juv=alph*1e15;	%% UV photon flux, in units of m^-2 s^-1?
	%Tph=1;		%% Thermal energy of photo electrons; for now I'm assuming Tph=1 eV for convenience.
	%Ip=qe*pi*a^2*Qabs*Y*juv*exp(-qe*qd/C/Tph);	%% only valid for q>0! (Shukla 2001)
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	Itot=Ii+Ie+Ip;

end