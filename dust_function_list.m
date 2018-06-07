function [f_x,f_deriv] = dust_function_list(case_label,x,eta,alph,Ti,Te,...
    e_mag,i_mag,C,lambda_D,lambda_i,a,w,species)
% % 3/11/13 note: NEED to add alph to the oml cases!!!

% 	Explanation of inputs:
%   case_label  = the equation you want to solve using bisection method or other
%                 purpose; this is a string specifying a charge model or other
%                 equation
%   x           = the number of charges on the grain, positive or negative
%   eta         = the ratio of ne/ni
%   alph        = coefficient of UV illumination, something like f_uv/ne/vthe
%   Ti          = ion temperature, in eV
%   Te          = electron temperature, in eV
%   e_mag       = dimensionless parameter, dust grain size divided by 
%                 electron Larmor radius, in meters
%   i_mag       = dimensionless parameter, dust grain size divided by 
%                 ion Larmor radius
%   C           = dust grain capacitance, in Farads
%   lambda_D    = linearized debye length, in meters
%   lambda_i    = mean free path of ions, in meters
%   a_grain     = the radius of the dust grain, in meters
%   w           = w is a two element array, given by [we wi], where we is the 
%                 grain speed relative to the electrons, and wi is the grain
%                 speed relative to the ions; both are in units of m/s
%   species     = this is the mass number of the ion species.
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Explanation of outputs:
% f_x     = value of the function; usually dimensionless current
% f_deriv = derivative of the function
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

Tau=Te/Ti;
mr=me/mi;
vthe=sqrt(2*qe*Te/me);	% local electron thermal speed, m/s
vthi=sqrt(mr/Tau)*vthe;	% local ion (proton) thermal speed, m/s

% % split up the w-array into:
% % we - grain velocity relative to an electron flow
% % wi - grain velocity relative to an ion flow
we=w(1);
wi=w(2);
% thermal mach numbers:
Me=we/vthe;
Mi=wi/vthi;



% % A useful list of functions that can be called by other external 
% % functions. This outputs the function and the derivative of the 
% % function; when calling this function from a main program just ignore 
% % f_deriv. case_label is a string input that must match one of these 
% % functions
switch case_label
    case '1a'
        f_x = exp(x).*cos(4*x);
        % figure out the derivative and do this later; 
        % just picking f_deriv=0 for now.
        f_deriv=0;
    case '1b'
        f_x = x.^(5/2);
        % figure out the derivative and do this later; 
        % just picking f_deriv=0 for now.
        f_deriv=0;
    case '1c'
        f_x = exp(cos(x));
        % figure out the derivative and do this later; 
        % just picking f_deriv=0 for now.
        f_deriv=0;
    % % gaussian function; for HW #3, problem 3 A)
    case 'gaussian'
        f_x = exp(-x.^2);
        % figure out the derivative and do this later; 
        % just picking f_deriv=0 for now.
        f_deriv=0;
    % % 1-v^2 function in the interval -1<v<1 for HW #3, problem 3 B)
    case '1-v^2'
        f_x = (x>=-1).*(x<=1).*(1-x.^2);
        % figure out the derivative and do this later; 
        % just picking f_deriv=0 for now.
        f_deriv=0;
    % % v*exp(-v) function for HW #3, problem 3 C), I called it EEDF 
    % % because it has the same functional form as an eedf.
    case 'eedf'
        f_x = x.*exp(-x); 
        % figure out the derivative and do this later; 
        % just picking f_deriv=0 for now.
        f_deriv=0;
    case 'transcendental'
        Ti=0.025;
        Te=1.6;
        E0=100;
        species=40; % Argon
        me=9.1e-31;
        mi=1.67e-27*species;
        % x_coord tells you where you are in space in order to calculate
        % the space potential. (in m)
        x_coord=-.1;
        % calculate Vspace in the prescribed way for a linear E-field
        V_space = 100*(x_coord); 
        % % x here represents surface potential.
        f_x=sqrt(Ti/mi)*(1-(x/Ti - V_space/Ti))-...
            sqrt(Te/me)*exp((x-V_space)/Te);
    case 'oml_monoenergetic_ions'
        if x>=0
            f_x=sqrt(pi)*(wi/vthe)/2*(1-qe*qe*x/C/(mi*wi.^2))-...
                eta*exp(qe*x/C/Te);
            f_deriv=0; 
        else
            f_x=sqrt(pi)*(wi/vthe)/2*exp(-qe*qe*x/C/(mi*wi.^2))-...
                eta*(1+qe*x/C/Te);
            f_deriv=0;        
        end
    case 'oml'	%%<----- Still need to put in photo-current!!!
        % % Note: x=elementary charges on dust grain.
    	% % both ions and electrons are unmagnetized if temp==1.
    	temp=e_mag<1 && i_mag<1;
        % if x<=0, then the grain is negatively charged.
        if temp==1;
            oml_case='unmag';
        end
        % % electrons = magnetized, ions = unmagnetized.
        temp=e_mag>1 && i_mag<1;
        if temp==1;
            oml_case='mag_e';
        end
        % % electrons = magnetized, ions = magnetized.
        temp=e_mag>1 && i_mag>1;
        if temp==1;
            oml_case='mag';
        end
    		% % the case where ions are magnetized and electrons are not so 
    		% % I'm should almost never happen, so I am not treating it.
        if x<=0;
            switch oml_case
    			case 'unmag'
                    if wi==0 || we==0
                        if wi==0 && we==0
                            f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)...
                                -eta*exp(qe*x/C/Te);
                            f_deriv=0;
                            
                        end
                        if wi==0 && we~=0
                            % first term comes from the ion current, second 
                            % comes from electron current
                            f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)...
                                -.5*eta/Me*(...
                                (Me.^2+.5+qe*x/C/Te)*sqrt(pi)*...
                                (erf(Me+sqrt(-qe*x/C/Te))+...
                                erf(Me-sqrt(-qe*x/C/Te)))+...
                                (sqrt(-qe*x/C/Te)+Me)*...
                                exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                                (sqrt(-qe*x/C/Te)-Me)*...
                                exp(-(Me+sqrt(-qe*x/C/Te)).^2));
                            f_deriv=0;
                        end
                        if wi~=0 && we==0
                            % first term is the electron term, the other
                            % terms comprise the ion term.
                            f_x=-eta*exp(qe*x/C/Te)...
                                +.5*sqrt(mr/Tau)*((Mi.^2+1/2-qe*x/C/Ti)*...
                                sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));
                            f_deriv=0;
                           
                        end
                    else 
                    % Lines below are for a flow shifted maxwellian 
                    % population of ions. This is how it is written in 1996
    				% Horanyi araa. The first set of terms is the ion term,
    				% which comprise the first 3 lines. The  electron term
    				% starts on the fourth line, and comprises the rest of
    				% the expression. 
                       f_x=.5*sqrt(mr/Tau)*(...
                           (Mi.^2+1/2-qe*x/C/Ti)*sqrt(pi)/Mi*erf(Mi)+...
                           exp(-Mi.^2))...
                           -.25*eta/Me*((Me.^2+.5+qe*x/C/Te)*sqrt(pi)*...
                           (erf(Me+sqrt(-qe*x/C/Te))+...
                           erf(Me-sqrt(-qe*x/C/Te)))+...
                           (sqrt(-qe*x/C/Te)+Me)*...
                           exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                           (sqrt(-qe*x/C/Te)-Me)*...
                           exp(-(Me+sqrt(-qe*x/C/Te)).^2));
    				% f_deriv is needed if you want to use this .m file for 
    				% the newton method, but dust_bisection.m works much 
    				% faster and with greater accuracy. Hence, I have just 
                    % set f_deriv=0 for this case.
    					f_deriv=0;
                        
                    end
	    		case 'mag_e'
	    			if wi==0 || we==0
                        if wi==0 && we==0
                            f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)...
                                -.5*eta*exp(qe*x/C/Te);
                            f_deriv=0;
                        end
                        if wi==0 && we~=0
                            % first term comes from the ion current, second 
                            % comes from electron current
                            f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)...
                                -.25*eta/Me*(...
                                (Me.^2+.5+qe*x/C/Te)*sqrt(pi)*...
                                (erf(Me+sqrt(-qe*x/C/Te))+...
                                erf(Me-sqrt(-qe*x/C/Te)))+...
                                (sqrt(-qe*x/C/Te)+Me)*...
                                exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                                (sqrt(-qe*x/C/Te)-Me)*...
                                exp(-(Me+sqrt(-qe*x/C/Te)).^2));
                            f_deriv=0;
                        end
                        if wi~=0 && we==0
                            % first term is the electron term, the other
                            % terms comprise the ion term.
                            f_x=-eta*exp(qe*x/C/Te)...
                                +sqrt(mr/Tau)*((Mi.^2+1/2-qe*x/C/Ti)*...
                                sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));
                            f_deriv=0;
                        end
                    else 
                    % Lines below are for a flow shifted maxwellian 
                    % population of ions. This is how it is written in 1996 
    				% Horanyi araa. The first set of terms is the ion term,
    				% which comprise the first 3 lines. The  electron term
    				% starts on the fourth line, and comprises the rest of
    				% the expression. The factor of 1/8 is needed due to 
                    % the factor of 1/2 that comes from electron 
                    % magnetization.
                       f_x=.5*sqrt(mr/Tau)*(...
                           (Mi.^2+1/2-qe*x/C/Ti)*sqrt(pi)/Mi*erf(Mi)+...
                           exp(-Mi.^2))...
                           -.125*eta/Me*((Me.^2+.5+qe*x/C/Te)*sqrt(pi)*...
                           (erf(Me+sqrt(-qe*x/C/Te))+...
                           erf(Me-sqrt(-qe*x/C/Te)))+...
                           (sqrt(-qe*x/C/Te)+Me)*...
                           exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                           (sqrt(-qe*x/C/Te)-Me)*...
                           exp(-(Me+sqrt(-qe*x/C/Te)).^2));
    				% f_deriv is needed if you want to use this .m file for 
    				% the newton method, but dust_bisection.m works much 
    				% faster and with greater accuracy. Hence, I have just 
                    % set f_deriv=0 for this case.
    					f_deriv=0;
	    			end
					
			
	    		case 'mag'
	    		% % I'll do this later.
	    		f_x=0;
	    		f_deriv=0;
            end
            
        % % POSITIVE grain potential equations go here.    
        else
            switch oml_case
    			case 'unmag'
    				if wi==0 || we==0
                        if wi==0 && we==0
                            f_x=sqrt(mr/Tau)*exp(-qe*x/C/Ti)...
                                -eta*(1+qe*x/C/Te);
                            f_deriv=0;
                        end
                        if wi==0 && we~=0
                            % first term comes from the ion current, second 
                            % comes from electron current
                            f_x=sqrt(mr/Tau)*exp(-qe*x/C/Ti)...
                                -.5*eta*((Me.^2+1/2+qe*x/C/Te)*...
                                sqrt(pi)/Me*erf(Me)+exp(-Me.^2));
                            f_deriv=0;
                        end
                        if wi~=0 && we==0
                            f_x=-eta*(1+qe*x/C/Te)...
                                +.5*sqrt(mr/Tau)/Mi*(...
                                (Mi.^2+.5-qe*x/C/Ti)*sqrt(pi)*...
                                (erf(Mi+sqrt(qe*x/C/Ti))+...
                                erf(Mi-sqrt(qe*x/C/Ti)))+...
                                Mi*(sqrt(qe*x/C/Ti/Mi)+1)*...
                                exp(-(Mi-sqrt(qe*x/C/Ti)).^2)-...
                                Mi*(sqrt(qe*x/C/Ti/Mi)-1)*...
                                exp(-(Mi+sqrt(qe*x/C/Ti)).^2));
                            f_deriv=0;
                        end
                    else 
                    % Lines below are for a flow shifted maxwellian 
                    % population of ions. This is how it is written in 1996
    				% Horanyi araa. The first set of terms is the electron 
    				% term, (1-qe*x/C/Ti) which comprise the first 3 lines.  
    				% The ion term starts on the fourth line, and comprises 
    				% the rest of the expression. 
                       f_x=-.5*eta*(...
                           (Me.^2+1/2+qe*x/C/Te)*sqrt(pi)/Me*erf(Me)+...
                           exp(-Me.^2))...
                           +.25*sqrt(mr/Tau)/Mi*(...
                           (Mi.^2+.5+qe*x/C/Ti)*sqrt(pi)*...
                           (erf(Mi+sqrt(qe*x/C/Ti))+...
                           erf(Mi-sqrt(qe*x/C/Ti)))+...
                           (sqrt(qe*x/C/Ti)+Mi)*...
                           exp(-(Mi-sqrt(qe*x/C/Ti)).^2)-...
                           (sqrt(qe*x/C/Ti)-Mi)*...
                           exp(-(Mi+sqrt(qe*x/C/Ti)).^2));
  
    				% f_deriv is needed if you want to use this .m file for 
    				% the newton method, but dust_bisection.m works much 
    				% faster and with greater accuracy. Hence, I have just 
                    % set f_deriv=0 for this case.
    					f_deriv=0;
    				end
    				
	    		case 'mag_e'
	    			if wi==0 || we==0
                        if wi==0 && we==0
                            f_x=sqrt(mr/Tau)*exp(-qe*x/C/Ti)...
                                -.5*eta*(1+qe*x/C/Te);
                            f_deriv=0;
                        end
                        if wi==0 && we~=0
                            % first term comes from the ion current, second 
                            % comes from electron current
                            f_x=sqrt(mr/Tau)*exp(-qe*x/C/Ti)...
                                -.25*eta*((Me.^2+1/2+qe*x/C/Te)*...
                                sqrt(pi)/Me*erf(Me)+exp(-Me.^2));
                            f_deriv=0;
                        end
                        if wi~=0 && we==0
                            f_x=-eta*(1+qe*x/C/Te)...
                                +sqrt(mr/Tau)/Mi*(...
                                (Mi.^2+.5-qe*x/C/Ti)*sqrt(pi)*...
                                (erf(Mi+sqrt(qe*x/C/Ti))+...
                                erf(Mi-sqrt(qe*x/C/Ti)))+...
                                (sqrt(qe*x/C/Ti)+Mi)*...
                                exp(-(Mi-sqrt(qe*x/C/Ti)).^2)-...
                                (sqrt(qe*x/C/Ti)-Mi)*...
                                exp(-(Mi+sqrt(qe*x/C/Ti)).^2));
                            f_deriv=0;
                        end
                    else  
                    % Lines below are for a flow shifted maxwellian 
                    % population of ions. This is how it is written in 1996
    				% Horanyi araa. The first set of terms is the electron 
    				% term, (1-qe*x/C/Ti) which comprise the first 3 lines.  
    				% The ion term starts on the fourth line, and comprises 
    				% the rest of the expression. 
                       f_x=-eta*(...
                           (Me.^2+1/2+qe*x/C/Te)*sqrt(pi)/Me*erf(Me)+...
                           exp(-Me.^2))...
                           +sqrt(mr/Tau)/Mi*(...
                           (Mi.^2+.5-qe*x/C/Ti)*sqrt(pi)*...
                           (erf(Mi+sqrt(qe*x/C/Ti))+...
                           erf(Mi-sqrt(qe*x/C/Ti)))+...
                           (sqrt(qe*x/C/Ti)+Mi)*...
                           exp(-(Mi-sqrt(qe*x/C/Ti)).^2)-...
                           (sqrt(qe*x/C/Ti)-Mi)*...
                           exp(-(Mi+sqrt(qe*x/C/Ti)).^2));
    				% f_deriv is needed if you want to use this .m file for 
    				% the newton method, but dust_bisection.m works much 
    				% faster and with greater accuracy. Hence, I have just 
                    % set f_deriv=0 for this case.
    					f_deriv=0;
                    end
	    		case 'mag'
	    		% % I'll do this later.
	    		f_x=0;
	    		f_deriv=0;
            end
        end
    case 'kortshagen'	%%<----- Still need to put in photo-current!!!
    % % kortshagen case will find the equilibrium charge for kortshagen 
	% % charge model. x=elementary charges on the grain. both ions and 
    % % electrons are unmagnetized if temp==1.
        if x<=0;
            if e_mag<1 && i_mag<1;
    			k_case='unmag';
            end
    	 	% % electrons = magnetized, ions = unmagnetized.
    		if e_mag>1 && i_mag<1;
    			k_case='mag_e';
    		end
    		% % electrons = magnetized, ions = magnetized.
    		if e_mag>1 && i_mag>1;
    			k_case='mag';
    		end
    		% % the case where ions are magnetized and electrons are not 
    		% % should almost never happen, so I'm not treating it.
            switch k_case
                case 'unmag'
    			
    				R0=a*abs(qe*x/C)*(1+a/lambda_D)...
                        /(1.5*Ti+abs(qe*x/C)*a/lambda_D);
    				% R0 can very easily be zero if the charge on the 
                    % grain is zero! There is no sheath yet; just use the 
                    % OML currents. The statement below is used to prevent  
                    % any problems;  also check to make sure lambda_i~=inf,
                    % or collisionless!
                    if R0==0 || lambda_i==inf
        			% If R0=0, use OML currents??
                        P0=1;
                        P1=0;
                        Pg1=0;
    					%disp('you are in the case R0=0')
                        if wi==0 || we==0
                            if wi==0 && we==0
                                f_x=-eta*exp(qe*x/C/Te)...
                                    +sqrt(mr/Tau)*(1-qe*x/C/Ti);
                                f_deriv=0;
                            end
                            if wi==0 && we~=0
                            % first term comes from the ion current, second 
                            % comes from electron current
                                f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)...
                                    -.5*eta/Me*(...
                                    (Me.^2+.5-qe*x/C/Te)*sqrt(pi)*...
                                    (erf(Me+sqrt(-qe*x/C/Te))+...
                                    erf(Me-sqrt(-qe*x/C/Te)))+...
                                    (sqrt(-qe*x/C/Te)+Me)*...
                                    exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                                    (sqrt(-qe*x/C/Te)-Me)*...
                                    exp(-(Me+sqrt(-qe*x/C/Te)).^2));
                                f_deriv=0;
                            end
                            if wi~=0 && we==0
                                f_x=-eta*exp(qe*x/C/Te)...
                                    +.5*sqrt(mr/Tau)*(...
                                    (Mi.^2+1/2-qe*x/C/Ti)...
                                    *sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));
                                f_deriv=0;
                            end
                            
                        % Lines below are for a flow shifted maxwellian 
                        % population of ions and electrons. This is how it 
                        % is written in 1996 Horanyi araa. 
                        else
                            f_x=.5*sqrt(mr/Tau)*(...
                                (Mi.^2+1/2-qe*x/C/Ti)*sqrt(pi)/Mi*...
                                erf(Mi)+ exp(-Mi.^2))...
                                -.25*eta/Me*(...
                                (Me.^2+.5+qe*x/C/Te)*sqrt(pi)*...
                                (erf(Me+sqrt(-qe*x/C/Te))+...
                                erf(Me-sqrt(-qe*x/C/Te)))+...
                                (sqrt(-qe*x/C/Te)+Me)*...
                                exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                                (sqrt(-qe*x/C/Te)-Me)*...
                                exp(-(Me+sqrt(-qe*x/C/Te)).^2));
  
                        % f_deriv is needed if you want to use this .m file for 
                        % the newton method, but dust_bisection.m works much 
                        % faster and with greater accuracy. Hence, I have just 
                        % set f_deriv=0 for this case.
                            f_deriv=0;
                        end
                            
                    % R0~=0, so there IS a capture radius.
                    else
                        % check to see if ions or electrons are flowing.
                        if wi==0 || we==0
                            if wi==0 && we==0
                                Ie=-eta*sqrt(Tau/mr)*exp(qe*x/C/Te); 
                            
                                Kn_R0 = lambda_i/(2*1.22*R0);
                                P0=exp(-1/(Kn_R0));
                                P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
                                Pg1=1-(P0+P1);
                		
                                % first term is the electron current, then  
                                % i_oml, then i_cec, and finally i_hyd 
                                Ioml=P0*(1-qe*x/C/Ti);
                                Icec=P1*(1.22*R0/a)^2;
                                Ihyd=Pg1*3*pi*(lambda_i)*...
                                    abs(qe*x)/C/Ti/a/4/sqrt(2);	
    			
                                f_x=Ie+Ioml+Icec+Ihyd;  
                        % % Do not need f_deriv for dust_bisection.m! Go 
    					% % back and finish if this is needed for doing the 
                        % % newton root-finding method.
                                f_deriv=0;
                            end
                            if wi==0 && we~=0
                                % electron current:
                             	Ie=-.25*eta/Me*sqrt(Tau/mr)*(...
                                    (Me.^2+.5+qe*x/C/Te)*sqrt(pi)*...
                                    (erf(Me+sqrt(-qe*x/C/Te))+...
                                    erf(Me-sqrt(-qe*x/C/Te)))+...
                                    (sqrt(-qe*x/C/Te)+Me)*...
                                    exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                                    (sqrt(-qe*x/C/Te)-Me)*...
                                    exp(-(Me+sqrt(-qe*x/C/Te)).^2));
                        	
                       
                                Kn_R0 = lambda_i/(2*1.22*R0);
                                P0=exp(-1/(Kn_R0));
                                P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
                                Pg1=1-(P0+P1);
                		
                                % i_oml, then i_cec, and finally i_hyd
                                Ioml=P0*(1-qe*x/C/Ti);
                                Icec=P1*(1.22*R0/a)^2;
                                Ihyd=Pg1*3*pi*(lambda_i)*...
                                    abs(qe*x)/C/Ti/a/4/sqrt(2);	
    			
                                f_x=Ie+Ioml+Icec+Ihyd;  
                        % % Do not need f_deriv for dust_bisection.m! Go 
    					% % back and finish if this is needed for doing the 
                        % % newton root-finding method.
                                f_deriv=0;
                            end
                            if wi~=0 && we==0
                              	Ie=-eta*sqrt(Tau/mr)*exp(qe*x/C/Te); 
                            
                                Kn_R0 = lambda_i/(2*1.22*R0);
                                P0=exp(-1/(Kn_R0));
                                P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
                                Pg1=1-(P0+P1);
                		
                                % first term is the electron current, then  
                                % i_oml, then i_cec, and finally i_hyd
                                
                                Ioml=.5*P0*((Mi.^2+1/2-qe*x/C/Ti)*...
                                    sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));
                                %Ioml=P0*(1-qe*x/C/Ti);
                                Icec=P1*(1.22*R0/a)^2;
                                Ihyd=Pg1*3*pi*(lambda_i)*...
                                    abs(qe*x)/C/Ti/a/4/sqrt(2);	
    			
                                f_x=Ie+Ioml+Icec+Ihyd;  
                        % % Do not need f_deriv for dust_bisection.m! Go 
    					% % back and finish if this is needed for doing the 
                        % % newton root-finding method.
                                f_deriv=0;
                            end
                            
                        % for this case, the capture radius exists and both
                        % ions and electrons are flowing.
                        else
                          	Ie=-.25*eta/Me*sqrt(Tau/mr)*(...
                                (Me.^2+.5+qe*x/C/Te)*sqrt(pi)*...
                                (erf(Me+sqrt(-qe*x/C/Te))+...
                                erf(Me-sqrt(-qe*x/C/Te)))+...
                                (sqrt(-qe*x/C/Te)+Me)*...
                            	exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                            	(sqrt(-qe*x/C/Te)-Me)*...
                             	exp(-(Me+sqrt(-qe*x/C/Te)).^2));
                            
                        	Kn_R0 = lambda_i/(2*1.22*R0);
                         	P0=exp(-1/(Kn_R0));
                           	P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
                          	Pg1=1-(P0+P1);
                		
                            % i_oml, then i_cec, and finally i_hyd
                                
                         	Ioml=.5*P0*((Mi.^2+1/2-qe*x/C/Ti)*...
                                sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2));
                         	%Ioml=P0*(1-qe*x/C/Ti);
                          	Icec=P1*(1.22*R0/a)^2;
                          	Ihyd=Pg1*3*pi*(lambda_i)*...
                                    abs(qe*x)/C/Ti/a/4/sqrt(2);	
    			
                           	f_x=Ie+Ioml+Icec+Ihyd;  
                        end
                    end
    				
    				
    				
    				
                % fix this case to allow for wi=0, and we=0 possibility!    
	    		case 'mag_e'
	    			R0=a*abs(qe*x/C)*(1+a/lambda_D)/...
                        (1.5*Ti+abs(qe*x/C)*a/lambda_D);
    				if R0==0
    					P0=1;
    					P1=0;
    					Pg1=0;
    					%disp('you are in the case R0=0')
                        if wi==0	
    					% first term comes from the ion current, second 
                        % comes from electron current
                            f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)-...
                                eta*exp(qe*x/C/Te);
                            f_deriv=-sqrt(mr*Tau)-eta*exp(qe*x/C/Te);
                        else
                        % Lines below are for a flow shifted maxwellian 
                        % population of ions. This is how it is written in  
                        % 1996 Horanyi 1996 araa. The first line below is 
                        % the ion term. The next 2 lines are the electron 
                        % term.
                            f_x=.5*sqrt(mr/Tau)*((Mi.^2+1/2-qe*x/C/Ti)*...
                                sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2))+...
                                 -.25*eta/Me*((Me.^2+.5+qe*x/C/Te)*...
                                 sqrt(pi)*(erf(Me+sqrt(-qe*x/C/Te))+...
                                 erf(Me-sqrt(-qe*x/C/Te)))+...
                                 (sqrt(-qe*x/C/Te)+Me)*...
                                 exp(-(Me-sqrt(-qe*x/C/Te)).^2)-...
                                 (sqrt(-qe*x/C/Te)-Me)*...
                                 exp(-(Me+sqrt(-qe*x/C/Te)).^2));
                
    					% f_deriv is needed if you want to use this .m file for the newton method,
    					% but dust_bisection.m works much faster and with greater accuracy.
    					% Hence, I have just set f_deriv=0 for this case.
                            f_deriv=0;
                        end
    				else
    					%Kn_a=lambda_i/a;
    					Kn_R0 = lambda_i/(2*1.22*R0);
    					P0=exp(-1/(Kn_R0));
        				P1=(1/(Kn_R0))*exp(-1/(Kn_R0));
                		Pg1=1-(P0+P1);
                		%disp('You are in the case R0!=0')
                		% % first term is the electron current, then i_oml, 
                        % % then i_cec, and finally i_hyd
    				
    					% Ie is reduced by factor of 1/2 for magnetization
    					Ie=-.5*eta*sqrt(Tau/mr)*exp(qe*x/C/Te);	
    					Ioml=P0*(1-qe*x/C/Ti);
    					Icec=P1*(1.22*R0/a)^2;
    					Ihyd=Pg1*3*pi*(lambda_i)*abs(qe*x)/C/Ti/a/4/sqrt(2);	
    			
    					f_x=Ie+Ioml+Icec+Ihyd;
	    			
	    				f_deriv=0;
                    end
    				
					
			
	    		case 'mag'
	    		% % I'll do this later.
                    f_x=0;
                    f_deriv=0;
                    
            % I think this hanging end below corresponds to the switch
            % statement for the different magnetization cases.
            end
            
        % % THIS IS CURRENTLY THE UNMAGNETIZED, Q_dust>0 OML CASE!!!
		% % THIS WILL NEED TO BE REDONE FOR KORTSHAGEN MODEL.    
        else
            f_x=-eta*(1+qe*x/C/Te)+sqrt(mr/Tau)*exp(-qe*x/C/Ti);
            f_deriv=-eta*sqrt(1/mr/Tau) - exp(-qe*x/C/Ti);
		
        end
    				
    case 'hutchinson'   
    % % fill this out for OML charge model to find equilibrium charge. x=elementary charges on dust grain.
    % % both ions and electrons are unmagnetized if temp==1.
        temp=e_mag<1 && i_mag<1;
        if x<=0
        % % negative grain potential.
            z=e_mag/(1+e_mag);

            % iota* in Patacchini and Hutchinson:
            iota=1-0.0946*z-0.305*z.^2+0.95*z.^3-2.2*z.^4+1.15*z.^5;
            % % FOR LAMBDA_D = FINITE, AND DEBYE-HUCKEL POTENTIAL:
            % eta, which is now dependent on grain sheath size:
            eta_mag=-qe*x/C/Te/e_mag*(1+e_mag/4*(1-exp(-4*a/lambda_D/e_mag)));

            % w, which is eta/(1+eta):
            w_mag=eta_mag/(1+eta_mag);

            % A, the fitting polynomial, a function of w:
            A_fit=0.678*w_mag+1.543*w_mag.^2-1.212*w_mag.^3;
            if i_mag<1
            % % unmagnetized ions:
                if Mi==0
                    % expression above Patacchini and Hutchinson, 2007
                    f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)-eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);
                    f_deriv=0;
                else
                    f_x=.5*sqrt(mr/Tau)*((Mi.^2+1/2-qe*x/C/Ti)*sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2))+...
                    -eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);			
                    % % see 1992_northrop_ps or 1981_whipple_repprogphys for the above.
                    % % Also, 1996 Horanyi and 1996 Northrop.
                end
            else
                % % magnetized ions:
                if Mi==0
                % % I should think carefully about Mi=0 vs Mi~=0. Does this
                % % change anything?
                    f_x=sqrt(pi)/2*sqrt(mr/Tau)+...
                    -eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);
                else	
                    f_x=sqrt(pi)/2*sqrt(mr/Tau)+...
                    -eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);
                    % % see 1992_northrop_ps or 1981_whipple_repprogphys for the above.
                    % % Also, 1996 Horanyi and 1996 Northrop.
                end
                
            end
        else
        % % positive grain potential.
        % % THIS IS CURRENTLY THE UNMAGNETIZED, Q_dust>0 OML CASE!!!
		% % THIS WILL NEED TO BE REDONE FOR HUTCHINSON MODEL.
            f_x=-eta*sqrt(Tau/mr)*(1+qe*x/C/Te)+exp(-x/C/Ti);
            f_deriv=-eta*sqrt(1/mr/Tau) - exp(-x/C/Ti);
        end

    case 'phgk'	% still need to actually finish this! 10/9/2013
    % % The Patacchini-Hutchinson and Gatti-Kortshagen model uses the Hutchinson electron current,
    % % and the Gatti-Kortshagen ion current. Consider a rewrite when I figure out how to introduce a
    % % flow shifted Maxwellian into all of this.
    
    % % fill this out for OML charge model to find equilibrium charge. x=elementary charges on dust grain.
    % % both ions and electrons are unmagnetized if temp==1.
        temp=e_mag<1 && i_mag<1;
        if x<=0
        % % negative grain potential.
            z=e_mag/(1+e_mag);

            % iota* in Patacchini and Hutchinson:
            iota=1-0.0946*z-0.305*z.^2+0.95*z.^3-2.2*z.^4+1.15*z.^5;
            % % FOR LAMBDA_D = FINITE, AND DEBYE-HUCKEL POTENTIAL:
            % eta, which is now dependent on grain sheath size:
            eta_mag=-qe*x/C/Te/e_mag*(1+e_mag/4*(1-exp(-4*a/lambda_D/e_mag)));

            % w, which is eta/(1+eta):
            w_mag=eta_mag/(1+eta_mag);

            % A, the fitting polynomial, a function of w:
            A_fit=0.678*w_mag+1.543*w_mag.^2-1.212*w_mag.^3;
            if i_mag<1
            % % unmagnetized ions:
                if Mi==0
                    % expression above Patacchini and Hutchinson, 2007
                    f_x=sqrt(mr/Tau)*(1-qe*x/C/Ti)-eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);
                    f_deriv=0;
                else
                    f_x=.5*sqrt(mr/Tau)*((Mi.^2+1/2-qe*x/C/Ti)*sqrt(pi)/Mi*erf(Mi)+exp(-Mi.^2))+...
                    -eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);			
                    % % see 1992_northrop_ps or 1981_whipple_repprogphys for the above.
                    % % Also, 1996 Horanyi and 1996 Northrop.
                end
            else
                % % magnetized ions:
                if Mi==0
                % % I should think carefully about Mi=0 vs Mi~=0. Does this
                % % change anything?
                    f_x=sqrt(pi)/2*sqrt(mr/Tau)+...
                    -eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);
                else	
                    f_x=sqrt(pi)/2*sqrt(mr/Tau)+...
                    -eta*(A_fit+(1-A_fit)*iota)*exp(qe*x/C/Te);
                    % % see 1992_northrop_ps or 1981_whipple_repprogphys for the above.
                    % % Also, 1996 Horanyi and 1996 Northrop.
                end
                
            end
        % % positive grain potential.
        % % THIS IS CURRENTLY THE UNMAGNETIZED, Q_dust>0 OML CASE!!!
		% % THIS WILL NEED TO BE REDONE FOR HUTCHINSON MODEL.    
        else
        
            f_x=-eta*sqrt(Tau/mr)*(1+qe*x/C/Te)+exp(-x/C/Ti);
            f_deriv=-eta*sqrt(1/mr/Tau) - exp(-x/C/Ti);
        end

    
    case 'nunomura'
        Ti=0.025;   %% does not matter anyway; ions are at bohm speed or greater
        Te=1.6;
        species=40; % Argon
        me=9.1e-31;
        mi=1.67e-27*species;
        qe=1.6e-19;
        %% you have to pick an x_coord, really a vertical position above
        %% the sheath in order to get V_space and E_ion, the kinetic energy
        %% of the ions at that spatial location (they are accelerated with
        %% decreasing x due to the electric field of the planar sheath.)
        x_coord =1;
        E_ion =1;
        V_space =1;
        
        %% x here represents surface potential.
        f_x = sqrt(qe*Te/mi)*(1-qe*x/(E_ion+abs(qe*V_space)))-sqrt(8*qe*Te/pi/me)*exp(V_space+x);  
	%% figure out the derivative and do this later; just picking f_deriv=0 for now.
        f_deriv=0;
%% error handling for incorrect string input:
% mistake=(case_label=='1a')||(case_label=='1b')||(case_label=='1c')||...
%     (case_label=='gaussian')||(case_label=='1-v^2')||(case_label=='eedf');
% if mistake==0
% 	exception = 'You must pick a valid function!';
% 	error(exception)
% end     
end



