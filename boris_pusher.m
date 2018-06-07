% % boris_pusher.m
% %
% % boris_pusher.m is a function that uses the boris algorithm to advance 
% % the positions and velocities of a particle. This function was written 
% % for the purpose of dust grains.

function [x,y,vx,vy,w]=...
    boris_pusher(dtNwt,md,q,x0,y0,vx0,vy0,E_x,E_y,B,nu_dn,g_x,g_y,...
    vex,vey,vix,viy,vnx,vny)
% below is the version of inputs that I intend to use. (sept 2013)
% I have probably more inputs than necessary, because I intend on treating 
% the linear approximation to ion drag 2010 Bacharis PRE??
%[x,y,vx,vy,w]=boris_pusher(dtNwt,md,q,x0,y0,vx0,vy0,E_x,E_y,B,g_x,g_y,ni,nu_dn,vnx,vny,vex,vey,vix,viy,Ti,lambda_D,ch_model);

%   explanation of inputs:
%	  dtNwt     = the newton timestep, in s
%   md        = dust grain mass, in kg
%   q         = dust grain charge, in coloumbs
%   x0        = dust grain x-position from last full time step, in m
%   y0        = dust grain y-position from last full time step, in m
%   vx0       = dust grain velocity in x-direction from last half time step, m
%   vy0       = dust grain velocity in y-direction from last half time step, m
%   Ex        = local electric field in x-direction, in V/m
%   Ey        = local electric field in y-direction, in V/m
%   B         = local magnetic field in z-direction, in T
%   gx        = local acceleration due to gravity in x-direction, in m/s^2
%   gy        = local acceleration due to gravity in y-direction, in m/s^2
%   vnx       = neutral gas velocity in m/s
%   vex       = electron flow velocity in x-direction in m/s
%   vey       = electron flow velocity in y-direction in m/s
%   vix       = ion flow velocity in x-direction in m/s
%   viy       = ion flow velocity in x-direction in m/s 
%   nu_dn     = dust neutral collision frequency
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   explanation of outputs:
% x   = new dust grain x-position at the full time step, in meters
% y   = new dust grain y-position at the full time step, in meters
% vx0 = new dust grain velocity in x-direction at the half time step, m/s
% vy0 = new dust grain velocity in y-direction at the half time step, m/s
% w   = a two element array, w=[we wi], the grain speed relative to 
%       [electron ion] flow, calculated at the half time step, in m/s
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % uncomment the line below if you want to turn neutral drag off.
%nu_dn=0;

% % compute grain speed relative to an electron flow:
we=sqrt((vx0-vex).^2+(vy0-vey).^2);
% % compute grain speed relative to an ion flow:
% % velocity of the grain is relative to the ions:
% % e.g., if vx = 0, but the ions are streaming towards the grain in the 
% % -x direction with velocity vix, then it is equivalent to the grain 
% % moving at a velocity +vix in the +x-direction.
wi=sqrt((vx0-vix).^2+(vy0-viy).^2);
% % make a w-vector; the first element is the grain speed relative to 
% % electron flow, the second is the grain speed relative to ion flow.
w=[we wi];

% % begin definition of some quantities for the linear ion drag. These are
% % from 2010 Baccharis PRE.
%beta_T=a*abs(q/C)/lambda_D/Ti; % from 1992 Barnes PRL and 2005 Fortov
% % modified Coloumb logarithm
%C_log_mod=-exp(beta_T/2)*expint(-beta_T/2); from 2005 Fortov    
% % Damping frequency from ion-collection force:
%nu_ic=pi*a.^2*mi*ni*sqrt(8*qe*Ti/pi/mi)*(1-q/C/Ti); 
% % Damping frequency from ion-orbit force: DOES NOT CURRENTLY HAVE THE
% % RIGHT UNITS!!! (OCT. 14 2013)
%nu_io=sqrt(32*pi)/3*sqrt(mi/2/qe/Ti)*eps0*Ti.^2*C_log_mod*beta_T.^2;

%%~~~~#1   	
% % first step in boris method: calculate vx,vy at the half step, apply 
% % half the electric impulse. Half step calculations of velocity are  
% % offset by dt/2 from quantities calculated at spatial locations;
% % velocities are dt/2 BEHIND position calculations. 
% % ~3/18/2013: I've make an adjustment to allow for Ex and Ey componenets; 
% % I think this is correct.
% % ~4/17/2013: I've included the ion drag force term, assuming that it  
% % has no dependence on the grain velocity. Essentially, this means that 
% % ion drag is being calculated at a specific spatial location.
% % ~7/18/2013: fi_x or fi_y are the ion drag force components, and
% % fn_x and fn_y are the neutral drag force components
vx_minus = vx0+dtNwt*q*E_x/2/md + dtNwt*g_x/2+...
    dtNwt*nu_dn*vnx/2;%+dtNwt*nu_ic*vix/2+dtNwt*nu_io*vix/2;
vy_minus = vy0+dtNwt*q*E_y/2/md + dtNwt*g_y/2+...
    dtNwt*nu_dn*vny/2;%+dtNwt*nu_ic*viy/2+dtNwt*nu_io*viy/2; 
% % Next four lines are old, vestigal code that includes drag force in a
% % different, but incorrect way.
%vx_minus = vx0+dtNwt*q*E_x/2/md +dtNwt*fi_x/2/md+dtNwt*fn_x/2/md+...
%dtNwt*g_x/2;
%vy_minus = vy0+dtNwt*q*E_y/2/md +dtNwt*fi_y/2/md+dtNwt*fn_y/2/md+...
%dtNwt*g_y/2;   
% A is simply a factor related to the gyro frequency
A = dtNwt*q*B/(2*md);   
v1 = ((1-dtNwt*(nu_dn)/2)*vx_minus+A*vy_minus);
v2 = ((1-dtNwt*(nu_dn)/2)*vy_minus-A*vx_minus);
% % at the end of this step, calculate the velocities and apply the
% % other half of the electric impulse
vx=((1+dtNwt*(nu_dn)/2)*v1+A*v2)/((1+dtNwt*(nu_dn)/2).^2+A.^2)+...
    dtNwt*q*E_x/2/md + dtNwt*g_x/2 +... 
    dtNwt*nu_dn*vnx/2 ;%+dtNwt*nu_ic*vix/2+dtNwt*nu_ic*vix/2;
vy=((1+dtNwt*(nu_dn)/2)*v2-A*v1)/((1+dtNwt*(nu_dn)/2).^2+A.^2)+...
    dtNwt*q*E_y/2/md + dtNwt*g_y/2 +... 
    dtNwt*nu_dn*vny/2 ;%+dtNwt*nu_ic*viy/2+dtNwt*nu_ic*viy/2; 

% % compute grain speed:
%w=sqrt(vx.^2+vy.^2);  
% % compute grain speed relative to an electron flow:
we=sqrt((vx-vex).^2+(vy-vey).^2);
% % compute grain speed relative to an ion flow:
% % velocity of the grain relative to the ions:
% % e.g., if vx = 0, but the ions are streaming towards the grain in the 
% % -x direction with velocity vix, then it is equivalent to the grain 
% % moving at a velocity +vix in the +x-direction.
wi=sqrt((vx-vix).^2+(vy-viy).^2);
% % make a w-vector; the first element is the grain speed relative to 
% % electron flow, the second is the grain speed relative to ion flow.
w=[we wi];
    
%%~~~~#2
% % second step in boris method: calculate positions at the full timestep
% % based on the velocities calculated at the half timestep. (The positions  
% % will be half a timestep ahead of the velocities)
x=dtNwt*vx+x0;
y=dtNwt*vy+y0;

end