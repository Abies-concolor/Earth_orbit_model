This is the Earth orbit model code as described in the following publication:

Kostadinov, T. S. and Gilb, R.: Earth Orbit v2.1: a 3-D visualization and analysis model of Earth's orbit, Milankovitch cycles and insolation, Geosci. Model Dev., 7, 1051-1068, doi:10.5194/gmd-7-1051-2014, 2014. 

The code posted here is the same as the code in the Supplement to this publication (http://www.geosci-model-dev.net/7/1051/2014/gmd-7-1051-2014-supplement.zip), with the exception of a minor fix to Earth_orbit_v2_1.m to correct a plot clipping issue.

For instructions of use, refer to the publication cited above and the original ReadMe.txt file (contents given here immediately below). 

Original ReadMe.txt file: -->  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++  Earth Orbit Model v2.1 Discalimer, License and Instructions for Use  +++++++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++  Authors: Dr. Tihomir Sabinov Kostadinov and Roy Gilb   ++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         
         (c) Dr. Tihomir Kostadiov and Roy Gilb
             University of Richmond
             University of California Santa Barbara
         
         //Author Contact Information: 
         Department of Geography and the Environment
         28 Westhampton Way
         University of Richmond
         Richmond, VA 23173
         USA
         E-mail: tkostadi@richmond.edu
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Reference: When available, consult the final revised version of: 
T. S. Kostadinov and R. Gilb (2013), Earth Orbit v2.1: a 3-D visualization and analysis model of Earth’s orbit, Milankovitch cycles and insolation,
Geosci. Model Dev. Discuss., 6, 5947-5980, 2013, www.geosci-model-dev-discuss.net/6/5947/2013/
doi:10.5194/gmdd-6-5947-2013

++++++++++++++++++++++++++++++++++++++++
+++++  Disclaimer & License ++++++++++++
++++++++++++++++++++++++++++++++++++++++
      (c) Tihomir Kostadinov & Roy Gilb, 2006-2014
               University of Richmond, VA, USA 
               University of California Santa Barbara, CA, USA
 

                                  ***Disclaimer: 
This software is supplied "as is". No warranty is given, express or impllied, of fitness for any purpose. Under no 
circumstances shall the authors or their instituions be liable to anyone for direct, indirect, incidental, consequential, 
special, exemplary, or any other kind of damages (however caused and on any theory of liability, 
and including damages incurred by third parties), arising from or relating to this software, or user's use, inability to use, or misuse of the 
software, or errors of the software.  
This software is not guaranteed to be error-free and is not meant to be used in any mission-critical applications. 
Use at your own risk and verify with other sources when appropriate. 
                                   
                                   *** License:
Unless superceding rules  of the author's institutions apply, this software is free to use, distribute and 
modify for non-commercial purposes under the Creative Commons BY-NC-SA 3.0 license. Authors should 
be notified at <tkostadi@richmond.edu> if any modified copies are to be distributed
to third parties. Any such modifications must be clearly documented and stated properly.     


++++++++++++++++++++++++++++++++++++++++
+++++  Installation  Instructions  +++++
++++++++++++++++++++++++++++++++++++++++
There is no installation per-se required. You need to have MATLAB(r) installed on your system.  Download file Earth_Orbit_v2_1.zip and unzip its contents into a
folder of your choice.  Either add the path to this folder to the MATLAB(r) path or set the current working directory to that folder.
In order to run the model, type the following at the MATLAB(r) command prompt: 
>> Earth_orbit_v2_1

and press Enter.  This should bring up the GUI with all the user controls.  These are fairly intuitive and self-explanatory.  Brief help is provided below.
For details, refer to the publication referenced above. 

If surface plots do not render properly, you may need to change the renderer used by MATLAB(r). At the prompt, type: 
>> set(gcf,'renderer','r_name')
where r_name can be opengl, painters or zbuffer. This will change the renderer for the current figure.

Note: If the GUI window is not sized properly for your system/monitor configuration (should be rare), you can change the window size (e.g. maximize it). If that 
does not work, you may need to use the "guide" tool in MATLAB(r) to resize/rearrange the GUI figure and resave it.  At the prompt, type
>> guide Earth_orbit_v2_1    

++++++++++++++++++++++++++++++++++++++++
+++++  Brief Model Use Help  +++++++++++
++++++++++++++++++++++++++++++++++++++++
Pressing the Help Button in the GUI should open this file in the MATLAB editor or another editor of your choice.
For details, refer to the publication referenced above. 

    (1) Initial Constants and Options: 
		Choose a real astronomical solution or the demo mode (user-selected Milankovitch parameters)
    (2) Milankovitch Orbital Parameters: 
        Select year since J2000 or values for the orbital parameters, depending on chosen mode. 
    
    (3) Calendar/Date, Latitude and Insolation Options:
		Select month and date, latitude (Southern hemisphere latitudes are negative),
		the TSI(S_o) value, as well as the calendar start date - either equinox
		is fixed to be March 20, or perihelion is fixed to be January 3.
		Press the "Plot/Update Orbit" button to produce/update a 3D orbital configuration plot 
		using the currently chosen input parameters, as well as update the text output in the GUI.
    
    (4) Time Series (TS) and Insolation Plotting Options:
		Depending on selected mode, various time series plots can be produced as indicated on each button.
		The plot(s) can be produced if the button is active. 
		For the astronomical solutions modes, users need to select the range of years to plot. Check the 
		corresponding box if you wish to save the plot data as ASCII files for further analysis.
    
    (5) Outputs:
		Numerical outputs update when the "Plot/Update Orbit button is pressed.


Summary of the Contents of Earth_Orbit_v2_1.zip:
--------------------------------------------------------------------------------
 FileName                     Explanation
--------------------------------------------------------------------------------
ReadMe.txt                    This file
Earth_orbit_v2_1.m            Main matlab m-file that needs to be called to run the model; spawns the main GUI control window
Earth_orbit_v2_1.fig          GUI figure associated with the above. Useful if user wished to modify the GUI.  
orbit.m                       Main function thet calculates and plots orbital position 
keplerian.m                   Solves forwardf Keplerian problem
keplerian_inverse.m           Solves inverse Keplerian problem 
insolation.m                  Calculates insolation   
insol_3d.m                    Calculates spatio-temporal distribution of insolation for a single year  
insol_TS.m                    Calculates insolation time series as a function of year since J2000 and latitude
generate_rot_m.m              Helper function to calculare rotation matrix for orbital geometry computation
Berger_orbpar.m               Calculates the Milankovitch parameters according to the Berger 1978 solution
getLaskar.m                   Read in and parse Laskar et al. (2004) Milankovitch parameter solution
orbit_output				  Directory in which numerical model output is saved. Examples are provided. 

INSOLP.LA2004.BTL.ASC         File provided by Jacques Laskar, see below for details
INSOLN.LA2004.BTL.100.ASC     File provided by Jacques Laskar, see below for details
Berger_data.dat               Contains output from Berger_orbpar.m
EPICA_CO2.dat                 EPICA CO2 time series, see file header for details. 
EPICA_deuterium.dat           EPICA deuterium/temperature time series, see file header for details. 
LR04_dO18_benthic_stack.dat   The Lisiecki and Raymo (2005) benthic delta-O-18 stack, see file header for details.  
Zachos01_d18O.dat             The Zachos et al. (2001) benthic delta-O-18 complitaion, see file header for details.
contemp_insol.dat             J2000 insolation pattern for the Laskar et al. (2004) solution and using S_o = 1,366 W m^-2
    

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++  ReadMe for the input files provided by Dr. Jacques Laskar:  ++++
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The following is an excerpt from the Readme.txt located at :
http://www.imcce.fr/Equipes/ASD/insola/earth/La2004/README.TXT
accessed August 1, 2013

This file, provided by Dr. Jacques Laskar at IMCCE, France, 
describes the contents of the following two ASCII files, 
used without modification by the Earth Orbit model: 
Filenames: 
INSOLN.LA2004.BTL.100.ASC
INSOLP.LA2004.BTL.ASC
These files are courtesy of Dr. Jacques Laskar and are available at: 
https://www.imcce.fr/Equipes/ASD/insola/earth/La2004/index.html
The files were downloaded from the above link on July 16, 2012

The main website for the paleoclimate astronomical solutions of Dr. Jacques Laskar is
at http://www.imcce.fr/Equipes/ASD/insola/earth/earth.html

===========================                         ============================
                         ASTRONOMIE ET SYSTEMES DYNAMIQUES
                         
                         INSTITUT DE MECANIQUE CELESTE

                                    la2004

                                2010, January 18 
===========================                         ============================
La2004    : 
==========
********************************************************************************
*  Authors: J. Laskar, M. Gastineau, F. Joutel                                 *
*  (c) Astronomie et Systemes Dynamiques, Institut de Mecanique Celeste,       *
*      Paris (2004)                                                            *
*                                                                              *
*                               Jacques Laskar                                 *
*                               Astronomie et Systemes Dynamiques,             *
*                               Institut de Mecanique Celeste                  *
*                               77 av. Denfert-Rochereau                       *
*                               75014 Paris                                    *
*                       email:  laskar@imcce.fr                                *
*                                                                              *
********************************************************************************


CONTAINS:
=========

--------------------------------------------------------------------------------
 FileName                 Explanations
--------------------------------------------------------------------------------
INSOLP.LA2004.BTL.ASC    Nominal solution La2004
                                    after present years (0 to +21Myr)
INSOLN.LA2004.BTL.100.ASC Nominal solution La2004,
                                    before present years (-101Myr to 0)


Byte-per-byte Description of file: INSOLN.LA2004.BTL.ASC
Byte-per-byte Description of file: INSOLP.LA2004.BTL.ASC
--------------------------------------------------------------------------------
   Bytes Format  Units   Label    Explanations
--------------------------------------------------------------------------------
   1-14   F13.3  1000yr  t        Time from J2000  in 1000 years
  18-39   D25.16 ---     e        eccentricity
  43-64   D25.16 rad     eps      obliquity (radians)
  68-89   D25.16 rad     pibar    longitude of perihelion from moving equinox
                                  (radians)
--------------------------------------------------------------------------------

Byte-per-byte Description of file: INSOLN.LA2004.BTL.100.ASC
Byte-per-byte Description of file: INSOLN.LA2004.BTL.250.ASC
--------------------------------------------------------------------------------
   Bytes Format  Units   Label    Explanations
--------------------------------------------------------------------------------
   1-9    F8.0   1000yr  t        Time from J2000  in 1000 years
  10-18   D9.6   ---     e        eccentricity
  19-27   D9.6   rad     eps      obliquity (radians)
  28-37   D10.6  rad     pibar    longitude of perihelion from moving equinox
                                  (radians)
--------------------------------------------------------------------------------

Description:

  La2004 is the nominal solution for precessional quantities and orbital quantities of 
  the Earth. The solution La2004 is provided with fortran subroutine in order to 
  compute the insolation quantities of Earth.
  


BIBLIOGRAPHY:
============



  La90: Laskar, J.: 1990, The chaotic motion of the solar system.
                  A numerical estimate of the chaotic zones
                  Icarus, 88, 266

  La93: Laskar, J., Joutel, F., Boudin, F.: 1993, Orbital, precessional
                  and insolation quantities for the Earth
                  from -20 Myr to + 10Myr
                  Astron. Astrophys. 270, 522

  La2004 :Laskar, J., Gastineau, M., Joutel, F., Robutel, P., Levrard, B., 
                  Correia, A.,
                  A long term numerical solution for the insolation quantities 
                  of Earth. {\it in preparation}



================================================================================
