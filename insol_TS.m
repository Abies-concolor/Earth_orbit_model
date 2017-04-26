function [sol, chosen_day_sol, annual_mean_sol] = insol_TS(a,AU,T,years,dayofyear, e,obliquity,precession,Fo,lat,month,day)
% [sol, chosen_day_sol, annual_mean_sol] = insol_TS(a,AU,T,years,dayofyear, e,obliquity,precession,Fo,lat,month,day)
%
% Loops through days of year (counting starts Jan. 1, and vernal equinox is fixed to be March 20, at the
% beginning of the day) and latitudes and calculated daily insolation for
% each latitude band and each day, given solutions of the Milankovitch
% orbital parameters for a given year, or demo/arbitrary values.

%See orbit.m for citations and more detailed explanation of input and output

% Dr. Tihomir S. Kostadinov, Roy Gilb, November 2006 - September 2013

if nargin ==0
    % !!! IMPORTANT!!! The constants and variables below are NOT normally used in the model,
    %instead they are passed as arguments by the caller functions of this function.
    %These are the default values provided here for testing of this
    %function in stand-alone mode.
    AU = 149.597870700; %in millions of km, 2013 Selected Astronomical Constants, from the online Astronomical Almanac at http://asa.usno.navy.mil/, section K
    a  = 1.00000261*AU; %(Standish, E. Myles; Williams, James C.. "Orbital Ephemerides of the
    %Sun, Moon, and Planets" (PDF). International Astronomical Union Commission 4: (Ephemerides).
    %e = 0.01670236225492288; %Laskar 2004 solution for J2000.0 (year 0  for him)
    obliquity = 0.4090928042223415; %radians,  %Laskar 2004 solution for J2000.0 (year 0  for him)
    precession = pi-1.796256991128036;%radians, Laskar 2004 solution for J2000.0 (year 0  for him); complementary angle calculated here for internal use only
    mean_anomaly = 90; %in degrees, CCW from perihelion, tied to time passage and via Kpler's II lawto true anomaly
    Fo = 1366;
    lat = 43;
    T = 365.256363; %Sidereal year length in days. (Should equinox-to-equinox year be used?)
    %T is prescribed a-priori, as Kepler's III Law is not in the model
    years = -500:20:500;
    e = 0.001:0.01:0.51;
    precession = precession*ones(size(e));
    obliquity = obliquity*ones(size(e));
    dayofyear = [1:5:365, 365];
else
    obliquity = obliquity*(pi/180);
    precession = precession*(pi/180);
end

years = years(:);
e = e(:);
precession = precession(:);
obliquity = obliquity(:);

if ~(length(years)==length(e) && length(years)==length(precession) && length(years)==length(obliquity))
    error('Time series of years and Milankovitch parameters passed to insol_TS.m are invalid');
end

dayofequinox = 31+28+19;%March 20, beginning; %78

switch month
        case 1
            ndays = 0;
        case 2
            ndays = 31;
        case 3
            ndays = 59;
        case 4
            ndays = 90;
        case 5
            ndays = 120;
        case 6
            ndays = 151;
        case 7
            ndays = 181;
        case 8
            ndays = 212;
        case 9
            ndays = 243;
        case 10
            ndays = 273;
        case 11
            ndays = 304;
        case 12
            ndays = 334;
        otherwise
            error('')
end
chosen_day = ndays + day; %day of year for the chosen date - used to display TS of ionsolation for a given date

sol = NaN(length(dayofyear),length(years));

for dd = 1:length(dayofyear)
    for yy = 1:length(years)     
        [time_since_perihelion, ~] = keplerian(T,e(yy),precession(yy)*(180/pi));
        
        days_since_spring = dayofyear(dd) - dayofequinox; %time units (days); -77
        if days_since_spring < 0
            days_since_spring = 365 + days_since_spring; %288
        end
        %determine mean anomaly of spring equinox - this wil be March 21 always and will depedn on precession's value;
        time_elapsed = mod(days_since_spring + time_since_perihelion-1,T); %288+77
        %!!!We are now making March 19 be a day that is 1.256363*86400 seconds long.
        %   starting the whole day count at the time of vernal equinox
        M = 360*(time_elapsed)/T;
           
        %begin snippet that is common to this script and insol_3d.m & orbit.m
        position_e = (pi/180)*keplerian_inverse(T,e(yy),M);
        
        r = a*(1-e(yy)^2)./(1+e(yy)*cos(position_e));
        r_vector = [r*cos(position_e), r*sin(position_e),0];
        tilt_m = generate_rot_m(obliquity(yy),precession(yy));
        k=[0 0 1]';
        kk = tilt_m*k;
        
        sun_declination = (180/pi)*acos(dot(r_vector, kk)/(norm(kk)*norm(r_vector)))-90;
        [sol(dd,yy), ~] = insolation(Fo,r,lat,sun_declination,AU);        
    end
end

chosen_day_sol = sol(chosen_day,:); %in W/m^2, averaged over a 24-hr period
annual_mean_sol = mean(sol); %Mean in W/m^2 over the year, approximated only with the represented days. 