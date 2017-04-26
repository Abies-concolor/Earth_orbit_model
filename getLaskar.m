function [ecc_sol, obliq_sol, precess_sol] = getLaskar(input_yr, Lask_neg, Lask_pos)
%[ecc_sol, obliq_sol, precess_sol] = getLaskar(input_yr, Lask_neg, Lask_pos)
%Returns the Milankovitch parameters from Laskar's 2004 solution, for agiven year; Uses linear interplation
%where given year falls between table entries.
%
%%%% INPUT %%%%%%%%%%%
%input_yr = year passed by the user from the laskar_year_text box -
%   serves as a search variable to find the orbital paramters in the file
% The following two parameters are passed by the caller, where the Laskar
% 2004 input files are loaded, as follows:
%Lask_pos = load('INSOLP.LA2004.BTL.ASC');%Solutions for the future
%Lask_neg = load('INSOLN.LA2004.BTL.100.ASC'); %Solutions for the past
% See %http://www.imcce.fr/Equipes/ASD/insola/earth/earth.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%% OUTPUT %%%%%%%%%
%Retrieves data from Laskar's data files and returns values for:
% ecc_sol - eccentricity
% obliq_so -  obliquity in degrees
% precess_sol - longitude of perihelion in degrees (This is omega_bar in
%   Berger et al. (2010), Fig. 1)
%%%%%%%%%%%%%%%%%%%%%%%%

%%Author: Roy Gilb
%Date: June 21, 2012

%Remove 0 from negative file
Lask_neg(1,:) = [];

%Flip the negative file
Lask_neg = flipud(Lask_neg);

%Concatenate the files
Laskar_solution = vertcat(Lask_neg, Lask_pos);

%Extract each column and assign variables

%Year values
table_year = Laskar_solution(:,1);

%Eccentricity values - no adjustment needed
ecc = Laskar_solution(:,2);
%obliquity values - no asjustment needed
obliq = (180/pi)*Laskar_solution(:,3);
%Need to adjust precession angle - find where it jumps from 360 to 0
precess = (180/pi)*Laskar_solution(:,4); %This is omega_bar

%Interp
ecc_sol = interp1(table_year, ecc, input_yr,'linear', NaN);
obliq_sol = interp1(table_year, obliq, input_yr, 'linear', NaN);

q = find(table_year==input_yr, 1,'first');
if ~isempty(q)
    precess_sol = precess(q);
else
    qs = find(table_year<input_yr,1,'last');
    qe = find(table_year>input_yr,1,'first');
    if abs(precess(qe)-precess(qs)) >=180 %a jump happens between 360 degrees and a new cycle startign near 0 degrees; 180 chosen to be arbitrarily big, much larger than
        %annual step of ~18 degrees/1,000 years
        precess_sol = interp1([table_year(qs); table_year(qe)],[precess(qs)-360; precess(qe)], input_yr,'linear', NaN);
        %precess(qs) converted to small negative value to deal with incorrect interpolation at jump; See Interpolation chapter in Meeus's Astronomical Algorithms
        if precess_sol<0
            precess_sol = precess_sol+360;
        end
    else
        precess_sol = interp1([table_year(qs); table_year(qe)],[precess(qs); precess(qe)], input_yr,'linear', NaN);
    end
end