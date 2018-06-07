% % improved_arctan.m
% % May 2013, Jeffrey Walker
% %
% % This function computes the arctangent of an (x,y) coordinate pair
% % originally in Cartesian coordinates. The atan function has been 
% % expanded here to allow for angles greater than 90 degrees
% % 
% % I think the function atan2 actually does the same thing, but I still use 
% % this anyway since I know exactly how it works

	
function [phi]=improved_arctan(x,y);
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Brief Explanation of inputs and how this program works
% %	x - x-coordinate
% %	y - y-coordinate
% % 	
% % The output phi will be between 0 and 2*pi RADIANS. If you want to have 
% % it in degrees,	you will have to do the necessary coding where/whenever 
% % you call this function.		
			
if x>=0 && y==0
	phi=0;
end
if x==0 && y>=0
	phi=pi/2;
end
if x<=0 && y==0
	phi=pi;
end
if x==0 && y<=0
	phi=1.5*pi;
end
	
	
if x>0 && y>0	
	quadrant='upper_right';
	phi=atan(y/x);
			
end
if x<0 && y>0
	quadrant='upper_left';
	phi=atan(abs(x)/y)+pi/2;
			
end
if x<0 && y<0
	quadrant='lower_left';
	phi=atan(abs(y)/abs(x))+pi;
			
end
if x>0 && y<0
	quadrant='lower_right';
	phi=atan(abs(x)/abs(y))+1.5*pi;
end