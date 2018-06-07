% % Jeffrey J. Walker
% % This is a simple program for finding the root of a specific quadratic 
% % using the bisection method.

function [c]= dust_bisection(domain,inp_func,guess,eta,alph,Ti,Te,...
    e_mag,i_mag,C,lambda_D,lambda_i,a_grain,w,species);
% 	Explanation of inputs:
% 	domain      = the size of the search domain; you must pick a size
%	  inp_func    = a string specifying which function to call from 
%                 dust_function_list
%                 c is the output of the function, which is/are the root(s).  
%                 +/- domain are the endpoints of the search
%   guess       = an initial guess for the dust grain charge, I usually just
%                 pick 0
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
% Explanation of output:
% c is the charge on the dust grain, given in units of elementary charges,
% NOT coloumbs!
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
vthi=sqrt(mr/Tau)*vthe;		% local ion (proton) thermal speed, m/s
%% error handling:
if sign(domain)<0
	exception = 'YOU MUST USE A POSITIVE DEFINITE VALUE FOR THE DOMAIN!!!';
	error(exception)
end

a=domain;
b=-domain;

% The idea in the following set of statements is to only look in a small 
% interval around the current value of the dust charge. It should not 
% change too much between updates, so look in this smaller interval.
%if guess ~=0
	%a=guess+50;
	%b=guess-50;
%else
	%a=domain;
	%b=-domain;
%end
% Nmax is the max number of iterations, to prevent an infinite loop
Nmax = 1e6;
% tol is the tolerance allowed. this could potentially be an input for the 
% function
tol = 1e-2;
% n is the iteration counter, initialized to 1
n=1;

    err = b-a;

	%while n<=Nmax && abs(err(i)) > tol
    while n<=Nmax && abs(err) > tol
	% calculate the error at the beginning of the loop. Use a and b 
	% instead of b and c, because if the error is less than the tolerance,
	% the current value of c from the last time through the loop will be
	% our root.
		%err(i) = 0.5*b(i)-0.5*a(i);
        err = 0.5*b-0.5*a;
	% chop the interval in half
		%c(i) = 0.5*a(i)+0.5*b(i);
        c = 0.5*a+0.5*b;
		%f_a = C*(a(i))^4 - B*(a(i))^2 + A;
		%f_c = C*(c(i))^4 - B*(c(i))^2 + A;
        [f_a] = dust_function_list(inp_func,a,eta,alph,Ti,Te,...
            e_mag,i_mag,C,lambda_D,lambda_i,a_grain,w,species);
        
       
        [f_c] = dust_function_list(inp_func,c,eta,alph,Ti,Te,...
            e_mag,i_mag,C,lambda_D,lambda_i,a_grain,w,species);
        
		if sign(f_a)==sign(f_c)
			
            			a=c;
		else
			
            			b=c;
		end
	
		n=n+1;
    end
    %disp(strcat('number of steps needed to converge:',' ',num2str(n)));
    %disp(c);
end