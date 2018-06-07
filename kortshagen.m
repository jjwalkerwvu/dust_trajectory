function [Itot,qf,Kn_R0,P0,P1,Pg1]=kortshagen(qflag,a,alph,Ti,ni,Te,ne,...
    B,Z,C,qd,lambda_D,lambda_i,w,species);
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
% I've decided to get rid of global vars; they are commented if you feel 
% like using them again.
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
    % Use below to match up with Patacchini and Hutchinson 2007
	Rle=sqrt(pi/4)*me*vthe/qe/B;	
	%Rli=mi*vthi/Z/qe/B;
    % Use below to match up with Patacchini and Hutchinson 2007
	Rli=sqrt(pi/4)*mi*vthi/Z/qe/B;	
	e_mag=a/Rle;
	i_mag=a/Rli;
end

% % The size of the dust grain and debye sheath combined
a_and_s=a+2.5*lambda_D;
% % The following may be more appropriate for magnetization ratios:
%e_mag=a_and_s/Rle;
%i_mag=a_and_s/Rli;

% now that charge, potential profile, and ion mean free path are known,
% the Capture Radius can be calculated at the mean ion thermal
% kinetic energy.
R0 = (abs(qd/C)*a*(1+a/lambda_D))/(1.5*Ti+abs(qd/C)*a/lambda_D);


% % calculate the ion mobility. confusingly, DO NOT MULTIPLY BY qe!!!
mu_i=(1/Ti)*(3*pi*vthi*lambda_i)/(16*sqrt(2));
if mu_i==inf;
    % obviously, if mu_i is infinite, or the plasma is essentially 
    % collisionless, we are in the OML regime so the hydrodynamic current 
    % is basically zero. To enforce this, just set mu_i=0 if lambda_i=inf
    mu_i=0;
end
    

% Calculate plasma currents to particle
if qd<=0	% negative dust potential (phi=qd/C)
    if e_mag<1		% UNMAGNETIZED ELECTRONS		
            if Me==0
            % derived via integration of Maxwellian from  
            % vmin=sqrt(2*qe*phi/me) to infinity, given 4*pi*a^2 collection 
            % area
                Ie=-ne*qe*(4*pi*a^2)/sqrt(4*pi)*vthe*exp(qd/C/Te);
                % expression above matches Patacchini and Hutchinson, 2007
            else
                Ie=-.5*sqrt(pi)*a.^2*ne*qe*vthe/Me*(...
                    (Me.^2+.5+qd/C/Te)*sqrt(pi)*(erf(Me+sqrt(-qd/C/Te))+...
                    erf(Me-sqrt(-qd/C/Te)))+...
                    (sqrt(-qd/C/Te)+Me)*exp(-(Me-sqrt(-qd/C/Te)).^2)-...
                    (sqrt(-qd/C/Te)-Me)*exp(-(Me+sqrt(-qd/C/Te)).^2));
            end
    % MAGNETIZED ELECTRONS        
    else
        	% same as unmagnetized case, except collection area reduced due 
        	% to  magnetization; may be off by some constant factor due to 
        	% cos-dependence of incidence angle
       		Ie=-ne*qe*(2*pi*a^2)*vthe/2/sqrt(pi)*exp((qd/C)/Te);
    end
    
	if i_mag<1					
        %%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        % KORTSHAGEN CURRENT IS A COMBINATION OF OML, HYD., AND CEC ION
        % CURRENTS.
        
        % compute the probabilities for various numbers of collisions an
        % ion undergoes in the dust grain sheath. P0 means no collisions,
        % P1 is one collsion, and Pg1 is greater than 1 collision.
       
        % R0 can very easily be zero if the charge on the grain is zero!
        % The statement below is used to prevent any problems
       
        if Mi==0
            if R0==0
                % If R0=0, use OML currents??
                P0=1;
                P1=0;
                Pg1=0;
                Kn_R0=0;
                ioml=ni*qe*(4*pi*a^2)*vthi/2/sqrt(pi)*(1-qd/C/Ti);
                icec=0;
                ihyd=0;
            		
    				
            else
                % compute the Knudsen Capture radius; see 2008 Gatti PRE  
                % for details
    		
                Kn_R0 = lambda_i/(2*1.22*R0); 
                P0=exp(-1/(Kn_R0));
                P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
                Pg1=1-(P0+P1);
                ioml=ni*qe*(4*pi*a^2)*vthi/2/sqrt(pi)*(1-qd/C/Ti);
                icec=ni*qe*(4*pi*(1.22*R0)^2)*vthi/2/sqrt(pi);
                % calculate the ion current due to hyd. charge model, from
                % Gotti et. al. Phys Rev E 2008, might want to find other 
                % sources too
                ihyd=16*pi*a*qe*ni*mu_i*abs((qd)/C)/2/sqrt(pi);
                					
            end
        		
        
        
            Ii=P0*ioml+P1*icec+Pg1*ihyd;
            if qflag==1
                % a charge between 0 and 1e4 elementary charges is a good 
                % search interval.
                [output]=dust_bisection(1e6,'kortshagen',(qd/qe),eta,...
                    alph,Ti,Te,e_mag,i_mag,C,lambda_D,lambda_i,a,w,...
                    species);
                qf=round(output);
                qf=qf*qe;
            else
                % return 0 if you do not wish to calculate equilibrium charge 
                % (qflag=0.)
                qf=0;
                % just return zero, value does not matter we just want to skip
                % the charge calculation step.
            end
            
        % Mi~=0, or for ion flow cases:    
        else
            if R0==0
                % If R0=0, use OML currents??
                P0=1;
                P1=0;
                Pg1=0;
                Kn_R0=0;
                % the best I can do for flowing ions right now is to just
                % modify the oml current:
                ioml=sqrt(pi)*a.^2*ni*Z*qe*vthi*(...
                    (Mi.^2+0.5-qd/C/Ti)*sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));
                %ioml=ni*qe*(4*pi*a^2)*vthi/2/sqrt(pi)*(1-qd/C/Ti);
                icec=0;
                ihyd=0;
            		
    				
            else
                % compute the Knudsen Capture radius; see 2008 Gatti PRE  
                % for details
    		
                Kn_R0 = lambda_i/(2*1.22*R0); 
                P0=exp(-1/(Kn_R0));
                P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
                Pg1=1-(P0+P1);
                % the best I can do for flowing ions right now is to just
                % modify the oml current:
                ioml=sqrt(pi)*a.^2*ni*Z*qe*vthi*(...
                    (Mi.^2+0.5-qd/C/Ti)*sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));
                %ioml=ni*qe*(4*pi*a^2)*vthi/2/sqrt(pi)*(1-qd/C/Ti);
                icec=ni*qe*(4*pi*(1.22*R0)^2)*vthi/2/sqrt(pi);
                % calculate the ion current due to hyd. charge model, from
                % Gotti et. al. Phys Rev E 2008, might want to find other 
                % sources too
                ihyd=16*pi*a*qe*ni*mu_i*abs((qd)/C)/2/sqrt(pi);
                5+5;					
            end
        		
        
        
            Ii=P0*ioml+P1*icec+Pg1*ihyd;
            if qflag==1
                % a charge between 0 and 1e4 elementary charges is a good 
                % search interval.
                [output]= dust_bisection(1e6,'kortshagen',(qd/qe),eta,...
                    alph,Ti,Te,e_mag,i_mag,C,lambda_D,lambda_i,a,w,...
                    species);
                qf=round(output);
                qf=qf*qe;
            else
                % return 0 if you do not wish to calculate equilibrium charge 
                % (qflag=0.)
                qf=0;
                % just return zero, value does not matter we just want to skip
                % the charge calculation step.
            end
            
        end
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 	% MAGNETIZED CASE?? (B -> inf)
	else				    
       		% simple thermal flux of ions, assuming ballistic trajectories 
        	% along  field lines w/ reduced collection area due to 
            % magnetization
      	Ii=ni*qe*(2*pi*a^2)*vthi;	% reference for this??
        if qflag==1
        			
			% a charge between 0 and 1e4 elementary charges is a good 
            % search interval.
			[output]= dust_bisection(1e6,'kortshagen',(qd/qe),eta,alph,...
                Ti,Te,e_mag,i_mag,C,lambda_D,lambda_i,a,w,species);
            disp(output)
			qf=round(output);
			qf=qf*qe;
		else
			% return 0 if you do not wich to calculate equilibrium charge 
            % (qflag=0.)
			qf=0;	
			% just return zero, value does not matter we just want to skip
			% the charge calculation step.
        end
        		
	end
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	ph_flux=alph*4*ni*eta*vthe;	%% just set ph_flux=4*n0*vthe for now!
	% FOR ENCELADUS, OR OTHER SOLAR SYSTEM RELEVANT SITUATIONS 
	% (Horanyi 1996): ph_flux=2.5e10 K/d/d; this has units cm^-2 s^-1
	% where d is the distance from the sun in AU
	% and K is the efficiency factor; ~1 for conductors and ~0.1 for 
    % dielectrics.
	Ip=qe*pi*a^2*ph_flux;		% only valid for q<0! (Horanyi 1996)
    % % Not sure why I have to multiply by sqrt(2*pi)... (august 2013)
    Ip=1/sqrt(4*pi)*Ip;
    % better expression: (4/24/2014)
    Ip=qe*ni*eta*vthe*alph*pi*a^2;
    
	Itot=Ii+Ie+Ip;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
% equations for positive dust potential:
else
	% these might need to be changed for Kortshagen charge model, but I 
    % have the expressions from Whipple 1981 here for now
	Ie = -ne*qe*(4*pi*a^2)*vthe/2/sqrt(pi)*(1+(qd/C)/Te);	
	Ii = ni*qe*(4*pi*a^2)*vthi/2/sqrt(pi)*exp(-(qd/C)/Ti);
    % These are from %% Whipple 1981, reviews of geophysics; also 
    % 1995 Cui and Goree IEEE
	%%~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%	%%	FOLLOWING NEEDS TO BE CHANGED FOR POSITIVE DUST POTENTIAL!!
	ph_flux=alph*4*ni*eta*vthe/2/sqrt(pi);	%% just set ph_flux=4*n0*vthe for now!
	% FOR ENCELADUS, OR OTHER SOLAR SYSTEM RELEVANT SITUATIONS 
    % (Horanyi 1996):
	% ph_flux=2.5e10 K/d/d; this has units cm^-2 s^-1
	% where d is the distance from the sun in AU
	% and K is the efficiency factor; ~1 for conductors and ~0.1 for 
    % dielectrics.
	Ip=qe*pi*a^2*ph_flux;		% only valid for q<0! (Horanyi 1996)
    % better expression: (4/24/2014)
    Ip=qe*ni*eta*vthe*alph*pi*a^2;
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	%Ip=qe*pi*a^2*ph_flux;		%% only valid for q<0! (Horanyi 1996)
	%Qabs=1;		% absorption efficiency is ~1 for 2*pi*a/lambda_uv > 1.
	%Y=1;		
    % the yield of photoelectrons; perhaps use Y=1? so we get 1 
    % photoelectron for every uv photon?
	%juv=alph*1e15;	% UV photon flux, in units of m^-2 s^-1?
	%Tph=1;		% Thermal energy of photo electrons; for now I'm assuming 
    %Tph=1 eV for convenience.
	%Ip=qe*pi*a^2*Qabs*Y*juv*exp(-qe*qd/C/Tph);	% only valid for q>0! 
    %(Shukla 2001)
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end
