function nu = keplerian_inverse(~,e,M)
% nu = keplerian_inverse(~,e,M)
% From Meeus Ch. 30, page 206 - Third method.
% Due to Roger Sinnott, Sky & Telescope, 1985
% Solves the inverse Kepler problem, i.e. given time of flight since last periapsis,
%     returns true anomaly
%
%%%%%% INPUT %%%%%%%
%T - orbital period in arbitrary time units (unused in this version,
%       replaced in argument list by ~)
%e - orbital eccentricity
%M - mean anomaly in degrees
%%%%%%%%%%%%%%%%%%%%
%
%%%%%%  OUTPUT %%%%%
% nu - true anomaly in degrees
%%%%%%%%%%%%%%%%%%%%
%
% Dr. T. S. Kostadinov, Nov. 7 , 2006 - January 2013

M = M*(pi/180);
N = 55; %Number of steps, this will ensure ~16 digits of accuracy, should be very sufficient
%(Meeus recommends 53 for 16-digit precision machine)

%%%%%%%% SINNOTT code(as published by MEEUS), translated from BASIC by T.S. Kostadinov:
F = sign(M);
M = abs(M)/(2*pi);
M = (M-floor(M))*2*pi*F;
if M < 0
    M = M+2*pi;
end
F = 1;
if M > pi
    F = -1;
    M = 2*pi-M;
end
Eo = pi/2;
D = pi/4;
for j = 1:N
    M1 = Eo-e*sin(Eo);
    Eo = Eo + D*sign(M-M1);
    D = D/2;
end
Eo = Eo*F;
%%%%%%%%%%%%%END OF SINNOTT CODE %%%%%%%%%%%%%%%%%%%%

nu = 2*atand(tan(Eo/2)*sqrt((1+e)/(1-e))); %Meeus Eq. 30.1

if nu < 0
    nu = nu + 360;
end
