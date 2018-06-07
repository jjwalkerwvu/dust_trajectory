% % 	charging_models.m
% %     Jeffrey Walker
% %	
		      
function [Itot,q,Kn_R0,P0,P1,Pg1]=...
    charging_models(qflag,ch_model,a,alph,Te,Ti,ne,ni,B,Z,C,qd,...
    lambda_D,lambda_i,w,species)
	
%   explanation of inputs:
%	  qflag 	= whether or not to evaluate equilibrium charge for a given 
%		          model; qflag=1 means do calculate q_eq, qflag=0 means do not 
%             calculate q_eq.
%	  ch_model        =   string which specifies charging model, current options
%                       (all listed in charging_models.m):
%                       'northrop','oml','simple','kortshagen','hutchinson',
%                       'oml_monoenergetic_ions', and 'constant_q'.
%   a       = dust grain size, in meters
%   alph    = coefficient of UV illumination, something like f_uv/ne/vthe
%   Te      = electron temperature, in eV
%   Ti      = ion temperature, in eV
%   ne      = electron density in m^-3
%   ni      = ion density in m^-3
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
% q       = charge on the dust grain, in coloumbs
% Kn_R0   = Knudsen capture radius of grain, dimensionless
% P0      = probability ion has no collisions in grain sheath
% P1      = probability ion has 1 collision in grain sheath
% Pg1     = probability ion has many collisions in grain sheath
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% some error handling, in case a is specified as a string 'electron' or
% 'ion':
if strcmp(a,'electron')==1 || strcmp(a,'ion')==1
    % just to make sure that if 'electron' or 'ion' is chosen for the dust
    % grain size, the charge model MUST be 'constant_q'.
    ch_model='constant_q';
end

% should be some error handling for qflag?
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


eta=ne/ni;
Tau=Te/Ti;
mr=me/mi;
vthe=sqrt(2*qe*Te/me);	% local electron thermal speed, m/s
vthi=sqrt(mr/Tau)*vthe;		% local ion (proton) thermal speed, m/s

% % short circuit the e_mag, i_mag check:
%e_mag=10;
%i_mag=10;

switch ch_model
    
    % I still have not written the "northrop" charging model function, but
    % it is the limiting case of high dust grain velocity relative ions
    % when compared to the ion thermal speed. Electrons are treated as
    % being much faster than the dust grain velocity relative to electrons.
    case 'northrop'
        [Itot,q]=northrop(a,alph,species,mi,...
            Ti,vthi,Rli,ni,Te,vthe,Rle,ne,B,Z,C,qd,er,vx,vy,V_space);
     
    case 'oml'
        [Itot,q,Kn_R0,P0,P1,Pg1]=oml(qflag,a,alph,Ti,ne,ni,Te,B,Z,C,qd,...
            lambda_D,lambda_i,w,species);

    case 'simple'
        [Itot,q] = simple(qflag,alph,ne,ni,eta,vthe,vthi,a,qd,C);
        
    case 'kortshagen'
        [Itot,q,Kn_R0,P0,P1,Pg1]=kortshagen(qflag,a,alph,Ti,ni,Te,ne,...
            B,Z,C,qd,lambda_D,lambda_i,w,species);
    case 'hutchinson'
        [Itot,q,Kn_R0,P0,P1,Pg1]=hutchinson(qflag,a,alph,species,...
            Ti,ne,ni,Te,B,Z,C,qd,lambda_D,lambda_i,w);
    case 'oml_monoenergetic_ions'
        [Itot,q,Kn_R0,P0,P1,Pg1]=oml_monoenergetic_ions(qflag,a,alph,Ti,ne,ni,Te,B,Z,C,qd,...
            lambda_D,lambda_i,w,species);
        
    % constant_q charge model fixes the grain charge to a given number.
    % This can be used to show particle trajectories for ions or electrons
    % in a profile specified by profiles.m.
    case 'constant_q'
        % not using Knudsen capture parameters for this profile.
     	Kn_R0=0;P0=0; P1=0;Pg1=0;
        % no current, because the charge is specified to be constant.
        Itot=0;
        if ischar(a)==1
            if strcmp(a,'electron')==1
                charge_state=-1;
            end
            if strcmp(a,'ion')==1
                charge_state=Z;
            end
        % IF a IS NOT A CHARACTER, WE ARE NOT TALKING ABOUT AN ELECTRON OR
        % AN ION, SO THIS IS WHERE YOU INPUT THE NUMBER OF ELECTRONS ON A 
        % GRAIN. USE A NEGATIVE NUMBER FOR charge_state TO INDICATE A
        % NEGATIVE CHARGE (charge_state IS THAT MANY ELECTRONS.)
        else
            % charge state, i.e., number of electrons. A negative number 
            % indicates positive charge, while a positive number indicates
            % negative charge.
            charge_state=Z;
            % use this instead?
            %charge_state=Z;
        end
        % Now specify the charge on the grain/elementary particle.
        q=qe*charge_state;
            

end

