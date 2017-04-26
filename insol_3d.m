function [sol,daylength] = insol_3d(a,AU,T,e,obliquity,precession,Fo,dayofyear,lats)
% [sol,daylength] = insol_3d(a,AU,T,e,obliquity,precession,Fo,dayofyear,lats)
%
%Built June 19-26, 2012, HJA, Blue river, Oregon, using code snippets from orbit.m
% Minor revisions February 2013, moved plotting to main GUI function
% Loops through days of year (counting starts Jan. 1, and vernal equinox is fixed to be March 20, at the 
% beginning of the day) and latitudes and calculated daily insolation for
% each latitude band and each day, given solutions of the Milankovitch
% orbital parameters for a given year, or demo/arbitrary values. 

%See orbit.m for citations and more detailed explanation of input and output

% Dr. Tihomir S. Kostadinov, Roy Gilb, November 2006 - January 2013

if nargin ==0
 % !!! IMPORTANT!!! The constants and variables below are NOT normally used in the model,
    %instead they are passed as arguments by the caller functions of this function. 
    %These are the default values provided here for easy use of this
    %function in stand-alone mode. 
    AU = 149.597870700; %in millions of km, 2013 Selected Astronomical Constants, from the online Astronomical Almanac at http://asa.usno.navy.mil/, section K
    a  = 1.00000261*AU; %(Standish, E. Myles; Williams, James C.. "Orbital Ephemerides of the
         %Sun, Moon, and Planets" (PDF). International Astronomical Union Commission 4: (Ephemerides).
    e = 0.01670236225492288; %Laskar 2004 solution for J2000.0 (year 0  for him)
    obliquity = 0.4090928042223415; %radians,  %Laskar 2004 solution for J2000.0 (year 0  for him)
    precession = pi-1.796256991128036;%radians, Laskar 2004 solution for J2000.0 (year 0  for him); complementary angle calculated here for internal use only
    mean_anomaly = 90; %in degrees, CCW from perihelion, tied to time passage and via Kpler's II lawto true anomaly
    Fo = 1366;
    dayofyear = [1:5:365, 365];
    lats = 90:-5:-90;
    T = 365.256363; %Sidereal year length in days. (Should equinox-to-equinox year be used?)
        %T is prescribed a-priori, as Kepler's III Law is not in the model
else
    obliquity = obliquity*(pi/180);
    precession = precession*(pi/180);
end

dayofequinox = 31+28+19;%March 20, beginning; %78
[time_since_perihelion, mean_anomaly_equinox] = keplerian(T,e,precession*(180/pi));

M = NaN*ones(size(dayofyear));
for i = 1:length(dayofyear)
    days_since_spring = dayofyear(i) - dayofequinox; %time units (days); -77
    if days_since_spring < 0
        days_since_spring = 365 + days_since_spring; %288
    end
    %determine mean anomaly of spring equinox - this wil be March 21 always and will depedn on precession's value;
    time_elapsed = mod(days_since_spring + time_since_perihelion-1,T); %288+77
    %!!!We are now making March 19 be a day that is 1.25636042*86400 seconds long.
    %   starting the whole day count at the time of vernal equinox
    M(i) = 360*(time_elapsed)/T;
end

%Preallocate arrays for speed
sol = NaN(length(lats),length(M));
daylength = sol;

for dd = 1:length(M)
    for ss = 1:length(lats)        
        position_e = (pi/180)*keplerian_inverse(T,e,M(dd));        
        r = a*(1-e^2)./(1+e*cos(position_e)); 
        r_vector = [r*cos(position_e), r*sin(position_e),0];
                
        tilt_m = generate_rot_m(obliquity,precession);
        k=[0 0 1]';
        kk = tilt_m*k;     
        sun_declination = (180/pi)*acos(dot(r_vector, kk)/(norm(kk)*norm(r_vector)))-90;      
        [sol(ss,dd), daylength(ss,dd)] = insolation(Fo,r,lats(ss),sun_declination,AU);
    end
end