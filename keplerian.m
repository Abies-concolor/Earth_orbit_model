function [t, M] = keplerian(T,e,nu)
% [t, M] = keplerian(T,e,nu)
% Solves the forward Kepler problem, i.e. given true anomaly, returns time of flight
%
%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%
%T - orbital period in arbitrary time units
%e - orbtal eccentricity
%nu - true anomaly in degrees - measured CCW from periapsis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%% OUTPUT %%%%%%%%%%%%%%
%t - period of time elapsed since last periapsis passage (perihelion for Earth)
%    Returned in the units of T. Also called time of flight.
%M - mean anomaly (in radians)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dr. T. S. Kostadinov, Nov. 7 , 2006 - January 2013

%%%%%%%%%%%% References/Works consulted %%%%%%%%%%%%%%%%%%%%%%%%
% Meeus, J. (1998), Astronomical Algorithms, Willmann-Bell, Richmond, VA. (2009 ed.)
% http://scienceworld.wolfram.com/physics/EccentricAnomaly.html
% http://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion
% http://en.wikipedia.org/wiki/Mean_anomaly
% http://en.wikipedia.org/wiki/Eccentric_anomaly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Treat special cases first, insure that T (instead of 0) is returned if nu==360
%This ensures that orbit.m code that calculates season lengths treats the
%cases when the longitude of perihelion is an integer multiple of pi/2
%radians correctly. 

%For nu>360 and negative nu's, return nu to remainder angle only, measured CCW from perihelion. 
if nu>360 | nu < 0
    nu = mod(nu,360);
end

if nu == 0
    t = 0; M = 0; 
elseif nu ==180
    t = T/2; M = pi; 
elseif nu ==360
    t = T; M = 2*pi; 
else  
    E = 2*atan(tand(nu/2)*sqrt((1-e)/(1+e))); %See Meeus Eq. 30.1
    
    M = E-e*sin(E); %Kepler's Equation, e.g. Meeus 30.5
    t = (M*T)/(2*pi);
    
    if M < 0 %t < 0 as well, orbit past aphelion
        M = M + 2*pi;
        t= T+t;
    end
end