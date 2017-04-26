function varargout = Earth_orbit_v2_1(varargin)
% EARTH_ORBIT_V2_1 M-file for Earth_orbit_v2_1.fig
%   This is the function that raises the main model GUI and provides all model controls. 
%   Normally, this is the only function that needs to be called on the command propmt in order
%   to run the model.  For details, see the model publication and the
%   ReadMe.txt file.  The code in the file below is a mix of MATLAB(r)
%   generated code and model author generated code.  Comments are provided
%   where needed. The external functions that are called are necessary for
%   the model to run and are also commented. T.K., Sept. 30, 2013. 
%
%
%      EARTH_ORBIT_V2_1, by itself, creates a new EARTH_ORBIT_V2_1 or raises the existing
%      singleton*.
%
%      H = EARTH_ORBIT_V2_1 returns the handle to a new EARTH_ORBIT_V2_1 or the handle to
%      the existing singleton*.
%
%      EARTH_ORBIT_V2_1('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in EARTH_ORBIT_V2_1.M with the given input arguments.
%
%      EARTH_ORBIT_V2_1('Property','Value',...) creates a new EARTH_ORBIT_V2_1 or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before orbit_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Earth_orbit_v2_1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Earth_orbit_v2_1

% Last Modified by GUIDE v2.5 18-Mar-2014 20:47:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Earth_orbit_v2_1_OpeningFcn, ...
    'gui_OutputFcn',  @Earth_orbit_v2_1_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Earth_orbit_v2_1 is made visible.
function Earth_orbit_v2_1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Earth_orbit_v2_1 (see VARARGIN)

% Choose default command line output for Earth_orbit_v2_1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles, false);

% UIWAIT makes Earth_orbit_v2_1 wait for user response (see UIRESUME)
% uiwait(handles.main_panel);


% --- Outputs from this function are returned to the command line.
function varargout = Earth_orbit_v2_1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function e_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function e_Callback(hObject, eventdata, handles)
% hObject    handle to e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e as text
%        str2double(get(hObject,'String')) returns contents of e as a double
e = str2double(get(hObject, 'String'));
if isnan(e) | ~(e>=get(handles.e_slider,'Min') & e<= get(handles.e_slider,'Max'))
    e = handles.data.e; %assume the old value
    set(hObject, 'String', e);
    errordlg(['Input must be a number b/n ', num2str(get(handles.e_slider,'Min')), ' and ' ,...
        num2str(get(handles.e_slider,'Max'))],'Error');
end
% Save the new e value
handles.data.e = e;
%update slider
set(handles.e_slider, 'Value',e);

guidata(hObject,handles);

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%cla;
try
    clf(1);
catch ME
    %warning('Tried to clear non-existing figure');
end

ecc = handles.data.e;
obliq = handles.data.obliquity;

omega_bar = handles.data.precession;
prec_for_orbit = 180-omega_bar; %The one used internally by keplerian and orbit.m
if prec_for_orbit < 0 %If this is not implemented, prec_for_orbit becomes negative and that
    %messes up some things within orbit.m, particularly the season length calculations
    prec_for_orbit = 360 + prec_for_orbit;
end

omega = 180+omega_bar; %longitude of perigee in Berger 2010 eliptical integrals paper
if omega>360
    omega = omega-360;          %adjust angle so it doesn't surpass 360 degrees
end


if handles.data.cal_mode == 2
    perihelion_dayofyear = 2;
    switch handles.data.month
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
    dayofyear = ndays + handles.data.day -1; %considering the beginnign of that day, so on Feb.1 00 Z
    %there have passed 31 das of teh year, on Dec. 1 - 334, and on Dec. 31 - 364!
    %Setting perihelion on the date of Jan. 3, 00 Z, this means it occurs after
    %two days have passed since Jan. 1, 00 Z.
    
    %!!!We are now making Dec. 31 be a day that is 1.256363*86400 seconds long
    time_elapsed = dayofyear - perihelion_dayofyear; %SINCE PERIHELION!
    if time_elapsed < 0
        time_elapsed = handles.data.period + time_elapsed;
    end
    M = 360*(time_elapsed/handles.data.period); %true anomaly in degrees
    mean_solar = -999; %Not implemented for this calendar case
elseif handles.data.cal_mode == 1
    %compute M for a calendar that is equinox-fixed on March 21 at 0 Z
    
    dayofequinox = 31+28+19;%March 20, beginning;
    
    switch handles.data.month
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
    dayofyear = ndays + handles.data.day-1;
    
    days_since_spring = dayofyear - dayofequinox; %time units (days);
    if days_since_spring < 0
        days_since_spring = 365 + days_since_spring;
    end
    
    mean_solar = 360*(days_since_spring/handles.data.period);
    
    %determine mean anomaly of spring equinox - this will be March 20 always and will depend on precession's value;
    [time_since_perihelion, mean_anomaly_equinox] = keplerian(...
        handles.data.period,ecc,prec_for_orbit);
    time_elapsed = mod(days_since_spring + time_since_perihelion, handles.data.period);
    %!!!We are now making March 19 be a day that is 1.256363*86400 seconds long.
    %starting the whole day count at the time of vernal equinox
    M = 360*(time_elapsed/handles.data.period);
else
    error('Unknown calendar mode. Specify equinox-fixed or perihelion-fixed')
end

%ready to call orbit function
f1h = figure(1);
clf;
set(f1h,'Name','3D Orbital Configuration','NumberTitle','off','Units','normalized','OuterPosition',[0.3 0.15 .7 .85]);
[r_vec_length, period, sun_dec, sol,season_length, daylength,true_anomaly] = orbit(...
    handles.data.sm_axis, handles.data.AU, handles.data.period, ecc, obliq, prec_for_orbit,M,...
    handles.data.Fo, handles.data.latitude);
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The following line was added May 20, 2015 to fix clipping issues that
%%%% occured with Matlab R2014 and later. This line is NOT in the GMD
%%%% publication!!!
set(gca,'Clipping','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate true solar longitude:
true_solar = mod(true_anomaly + omega,360); %From Berger 2010 elliptical paper, Eq. 6

%Write output to the GUI outputs box
set(handles.sun_dec, 'String',sprintf('%5.3f',sun_dec));
set(handles.insolation, 'String',sprintf('%8.4f',sol));
set(handles.daylength, 'String', sprintf('%5.2f',daylength));
set(handles.rvec_length, 'String', sprintf('%12.9f',r_vec_length));
set(handles.distance_au, 'String',sprintf('%12.9f',r_vec_length/handles.data.AU));
set(handles.spring_len, 'String',sprintf('%5.2f',season_length(1)));
set(handles.summer_len, 'String',sprintf('%5.2f',season_length(2)));
set(handles.fall_len, 'String',sprintf('%5.2f',season_length(3)));
set(handles.winter_len, 'String',sprintf('%5.2f',season_length(4)));

%Solar Longitudes
set(handles.mean_solar_text, 'String', sprintf('%7.4f',mean_solar));
set(handles.trueSolar_text, 'String', sprintf('%7.4f', true_solar));

set(handles.e, 'String', sprintf('%8.7f',ecc));
set(handles.obliquity, 'String', sprintf('%6.4f',obliq));
set(handles.precession, 'String', sprintf('%6.4f',omega_bar));
set(handles.lon_perihelion, 'String',sprintf('%3.2f',omega));

guidata(hObject, handles);


% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset)

%Default values
handles.data.AU = 149.597870700; %in millions of km, Astronomical Almanac from aa.usno.navy.mil
handles.data.sm_axis = 1.00000261*handles.data.AU;
%(Standish, E. Myles; Williams, James C.. "Orbital Ephemerides of the
%Sun, Moon, and Planets" (PDF). International Astronomical Union Commission 4: (Ephemerides). Retrieved 2010-04-03. See table 8.10.2.
handles.data.e = 0.01670236225492288; %Laskar 2004 solution for J2000.0 (year 0  for him)
handles.data.obliquity = 0.4090928042223415*(180/pi); %Laskar 2004 solution for J2000.0 (year 0  for him);
handles.data.precession = 1.796256991128036*(180/pi); %Laskar 2004 solution for J2000.0 (year 0  for him)
handles.data.latitude = 43; %degrees
%set initial year values
handles.data.myear = 0;
handles.data.laskar_year_text = 0;
handles.data.period = 365.256363; %Meeus Astronomical Algortihms Appendix I
handles.data.Fo = 1366; %W/m^2, http://www.pmodwrc.ch/pmod.php?topic=tsi/composite/SolarConstant
%Initialize date of model to current date taken from computer clock
handles.data.month = str2double(datestr(now,5));
handles.data.day = str2double(datestr(now,7));
handles.data.cal_mode = 1; %fixed date of equinox
handles.data.solutions_mode = 0; %demo

handles.data.start_year_text = -200;
handles.data.end_year_text = 200;
handles.data.save_insol_data = 0; %Default value is to not save the plot data

%Set initial field values
set(handles.sm_axis, 'String', sprintf('%12.9f',handles.data.sm_axis));
set(handles.myear, 'String', handles.data.myear);
set(handles.laskar_year_text, 'String', handles.data.laskar_year_text);
set(handles.e, 'String', handles.data.e);
set(handles.obliquity, 'String', handles.data.obliquity);
set(handles.precession, 'String', handles.data.precession);
set(handles.latitude, 'String', handles.data.latitude);
set(handles.period, 'String', sprintf('%12.7f',handles.data.period));
set(handles.Fo, 'String', handles.data.Fo);
set(handles.month, 'Value', str2double(datestr(now,5)));
set(handles.day, 'Value',str2double(datestr(now,7)));
set(handles.calendar_mode, 'Value', 1);
set(handles.start_year_text, 'String', handles.data.start_year_text);
set(handles.end_year_text, 'String', handles.data.end_year_text);

set(handles.saving_data_checkbox,'Value',get(handles.saving_data_checkbox,'Min')); %Default checkbox state

%Set sliders properties
set(handles.e_slider, 'Min',0.005);
set(handles.e_slider, 'Max', 0.9);
set(handles.e_slider, 'Value', handles.data.e);

set(handles.obliquity_slider, 'Min',0);
set(handles.obliquity_slider, 'Max', 60);
set(handles.obliquity_slider, 'Value', handles.data.obliquity);

set(handles.precession_slider, 'Min',0);
set(handles.precession_slider, 'Max', 360);
set(handles.precession_slider, 'Value', handles.data.precession);

%Set values for year slider
set(handles.year_slider, 'Min', -1000);   %minimum value allowed by Berger solution in kiloyears since present
set(handles.year_slider, 'Max', 1000);    %maximum value allowed by Berger solution
set(handles.year_slider, 'Value', handles.data.myear);

set(handles.laskar_year_slider, 'Min', -101000);   %minimum value allowed by Laskar solution in kiloyears since present
set(handles.laskar_year_slider, 'Max', 21000);    %maximum value allowed by Laskar solution
set(handles.laskar_year_slider, 'Value', handles.data.laskar_year_text);

%Disable year text box and slider for Berger and Laskar- not part of the default demo
set(handles.year_slider, 'Enable', 'off');
set(handles.myear, 'Enable', 'off');

%Disable laskar slider and text - demo mode is default
set(handles.laskar_year_slider, 'Enable', 'off');
set(handles.laskar_year_text, 'Enable', 'off');

%Disable Time series options - only allowed with Berger or Laskar real astronomical solution
set(handles.timeSeries_button, 'Enable', 'off','Value',0);
set(handles.start_year_text, 'Enable', 'off');
set(handles.end_year_text, 'Enable', 'off');
set(handles.paleo_data_plot_button,'Enable','off','Value',0);

%Load Laskar's input files into the handles structure - will be used by getLaskar within the
%laskar_year_text and laskar_year_slider callbacks
handles.data.laskar_pos = load('INSOLP.LA2004.BTL.ASC');
handles.data.laskar_neg = load('INSOLN.LA2004.BTL.100.ASC');

% Update handles structure
guidata(handles.main_panel, handles);


% --- Executes on key press over calculate with no controls selected.
function calculate_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on slider movement.
function e_slider_Callback(hObject, eventdata, handles)
% hObject    handle to e_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
a = get(hObject, 'Value');
%update value of corresponding field
set(handles.e, 'String', a);
%update value of e;
handles.data.e = a;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function e_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function obliquity_Callback(hObject, eventdata, handles)
% hObject    handle to obliquity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of obliquity as text
%        str2double(get(hObject,'String')) returns contents of obliquity as a double
e = str2double(get(hObject, 'String'));
if isnan(e) | ~(e>=get(handles.obliquity_slider,'Min') & e<= get(handles.obliquity_slider,'Max'))
    e = handles.data.obliquity; %assume the old value
    set(hObject, 'String', e);
    errordlg(['Input must be a number b/n', num2str(get(handles.obliquity_slider,'Min')), ' and ' ,...
        num2str(get(handles.obliquity_slider,'Max'))],'Error');
end
% Save the new e value
handles.data.obliquity = e;
%update slider
set(handles.obliquity_slider, 'Value',e);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function obliquity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to obliquity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function precession_Callback(hObject, eventdata, handles)
% hObject    handle to precession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of precession as text
%        str2double(get(hObject,'String')) returns contents of precession as a double
e = str2double(get(hObject, 'String'));
if isnan(e) | ~(e>=get(handles.precession_slider,'Min') & e<= get(handles.precession_slider,'Max'))
    e = handles.data.precession; %assume the old value
    set(hObject, 'String', e);
    errordlg(['Input must be a number b/n', num2str(get(handles.precession_slider,'Min')), ' and ' ,...
        num2str(get(handles.precession_slider,'Max'))],'Error');
end
% Save the new e value
handles.data.precession = e;
%update slider
set(handles.precession_slider, 'Value',e);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function precession_CreateFcn(hObject, eventdata, handles)
% hObject    handle to precession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function latitude_Callback(hObject, eventdata, handles)
% hObject    handle to latitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latitude as text
%        str2double(get(hObject,'String')) returns contents of latitude as a double
e = str2double(get(hObject, 'String'));
if isnan(e) | ~(e>=-90 & e<= 90)
    e = handles.data.latitude; %assume the old value
    set(hObject, 'String', e);
    errordlg('Input must be a number b/n -90 and +90','Error');
end
% Save the new e value
handles.data.latitude = e;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function latitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function obliquity_slider_Callback(hObject, eventdata, handles)
% hObject    handle to obliquity_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
a = get(hObject, 'Value');
%update value of corresponding field
set(handles.obliquity, 'String', a);
%update value of e;
handles.data.obliquity = a;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function obliquity_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to obliquity_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function precession_slider_Callback(hObject, eventdata, handles)
% hObject    handle to precession_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
a = get(hObject, 'Value');
%update value of corresponding field
set(handles.precession, 'String', a);
%update value of e;
handles.data.precession = a;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function precession_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to precession_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function Fo_Callback(hObject, eventdata, handles)
% hObject    handle to Fo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fo as text
%        str2double(get(hObject,'String')) returns contents of Fo as a double
e = str2double(get(hObject, 'String'));
if isnan(e) | ~(e>=500 & e<= 2500)
    e = handles.data.Fo; %assume the old value
    set(hObject, 'String', e);
    errordlg('Input must be a number b/n 500 and 2500','Error');
end
% Save the new e value
handles.data.Fo = e;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Fo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in month.
function month_Callback(hObject, eventdata, handles)
% hObject    handle to month (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns month contents as cell array
%        contents{get(hObject,'Value')} returns selected item from month
contents = get(hObject,'String');
mo = contents{get(hObject,'Value')};
switch mo
    case 'January'
        handles.data.month = 1;
    case 'February'
        handles.data.month = 2;
    case 'March'
        handles.data.month = 3;
    case 'April'
        handles.data.month = 4;
    case 'May'
        handles.data.month = 5;
    case 'June'
        handles.data.month = 6;
    case 'July'
        handles.data.month = 7;
    case 'August'
        handles.data.month = 8;
    case 'September'
        handles.data.month = 9;
    case 'October'
        handles.data.month = 10;
    case 'November'
        handles.data.month = 11;
    case 'December'
        handles.data.month = 12;
    otherwise
        error('');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function month_CreateFcn(hObject, eventdata, handles)
% hObject    handle to month (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in day.
function day_Callback(hObject, eventdata, handles)
% hObject    handle to day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns day contents as cell array
%        contents{get(hObject,'Value')} returns selected item from day
contents = get(hObject,'String');
day = str2num(contents{get(hObject,'Value')});

switch handles.data.month
    case 1
        handles.data.day = day;
    case 2
        if day == 29
            set(hObject, 'Value', 1);
            errordlg('Invalid day selection. No leap years are allowed in this demo. Dec. 31 is assumed slightly longer instead.  Assuming day 1','Error');
            handles.data.day = 1;
        elseif day > 29
            set(hObject, 'Value', 1);
            errordlg('Invalid day selection. Assuming day 1','Error');
            handles.data.day = 1;
        else
            handles.data.day = day;
        end
    case 3
        handles.data.day = day;
    case 4
        if day > 30
            set(hObject, 'Value', 1);
            errordlg('Invalid day selection. Assuming day 1','Error');
            handles.data.day = 1;
        else
            handles.data.day = day;
        end
    case 5
        handles.data.day = day;
    case 6
        if day > 30
            set(hObject, 'Value', 1);
            errordlg('Invalid day selection. Assuming day 1','Error');
            handles.data.day = 1;
        else
            handles.data.day = day;
        end
    case 7
        handles.data.day = day;
    case 8
        handles.data.day = day;
    case 9
        if day > 30
            set(hObject, 'Value', 1);
            errordlg('Invalid day selection. Assuming day 1','Error');
            handles.data.day = 1;
        else
            handles.data.day = day;
        end
    case 10
        handles.data.day = day;
    case 11
        if day > 30
            set(hObject, 'Value', 1);
            errordlg('Invalid day selection. Assuming day 1','Error');
            handles.data.day = 1;
        else
            handles.data.day = day;
        end
    case 12
        handles.data.day = day;
    otherwise
        error('');
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function day_CreateFcn(hObject, eventdata, handles)
% hObject    handle to day (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in calendar_mode.
function calendar_mode_Callback(hObject, eventdata, handles)
% hObject    handle to calendar_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns calendar_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calendar_mode
contents = get(hObject,'String');
selection = contents{get(hObject,'Value')};
switch selection
    case 'Perihelion: Jan. 3'
        handles.data.cal_mode = 2;
    case 'Vernal equinox: Mar. 20'
        handles.data.cal_mode = 1;
    otherwise
        error('');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function calendar_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calendar_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in close_btn.
function close_btn_Callback(hObject, eventdata, handles)
% hObject    handle to close_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

user_response = questdlg('Really Exit Program?');
switch user_response
    case {'No'}
    case {'Cancel'}
    case {'Yes'}
        delete(handles.main_panel)
        close all
end


% --- Executes during object creation, after setting all properties.
function main_panel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function main_panel_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to main_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when main_panel is resized.
function main_panel_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to main_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close main_panel.
function main_panel_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on slider movement.
function year_slider_Callback(hObject, eventdata, handles)
% hObject    handle to year_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
a = get(hObject, 'Value');
%update value of corresponding field
set(handles.myear, 'String', a);
%update value of e;
handles.data.myear = a;

%update Milankovitch parameters
%calculate

switch handles.data.solutions_mode
    case 0                       %Demo solution chosen
        %msgbox('Demo method')
        %No need to do anything here...
        error('Calling year slider when demo is chosen is incorrect...')
        
    case 1                       %Berger solution chosen
        %msgbox('Berger method')
        %Year between -999999 - +999999 as input by the user
        year_to_plot = handles.data.myear*1000;
        
        ecc = -999;     %*ones(size(t));
        obliq = ecc; omega_bar = ecc;
        
        %Call Berger's orbpar solution here
        [ecc, obliq, omega] = Berger_orbpar(year_to_plot+2000); %Fix to J2000.0 convention (add 2000)
        
        %Adjust the obliquity and omegar_bar values to the convention
        obliq = obliq*180/pi;
        omega_bar = omega*(180/pi)-180;
        if omega_bar <0
            omega_bar = 360+omega_bar;
        end
        %prec_for_orbit = 360-omega_bar*180/pi;
        
        handles.data.e = ecc;
        handles.data.precession = omega_bar;
        handles.data.obliquity = obliq;
        
    case 2                          %Laskar solution chosen
        error('Wrong solution slider chosen...')
    otherwise
        msgbox('Invalid solutions method')
end

set(handles.e, 'String', sprintf('%8.7f',handles.data.e));
set(handles.obliquity, 'String', sprintf('%6.4f',handles.data.obliquity));
set(handles.precession, 'String', sprintf('%6.4f',handles.data.precession));

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function year_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to year_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%Callback for the Berger solution
function myear_Callback(hObject, eventdata, handles)
% hObject    handle to myear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of myear as text
%        str2double(get(hObject,'String')) returns contents of myear as a double
e = str2double(get(hObject, 'String'));
if isnan(e) || ~(e>=get(handles.year_slider,'Min') && e<= get(handles.year_slider,'Max'))
    e = handles.data.myear; %assume the old value
    set(hObject, 'String', e);
    errordlg(['Input must be a number b/n ', num2str(get(handles.year_slider,'Min')), ' and ' ,...
        num2str(get(handles.year_slider,'Max'))],'Error');
end
% Save the new e value
handles.data.myear = e;
%update slider
set(handles.year_slider, 'Value',e);


switch handles.data.solutions_mode
    case 0                       %Demo solution chosen
        %msgbox('Demo method')
        %No need to do anything here...
        error('Calling year slider when demo is chosen is incorrect...')
        
    case 1                       %Berger solution chosen
        %msgbox('Berger method')
        %Year between -999999 - +999999 as input by the user
        year_to_plot = handles.data.myear*1000; %The sliders are now in kiloyears
        
        ecc = -999;     %*ones(size(t));
        obliq = ecc; omega_bar = ecc;
        
        %Call Berger's orbpar solution here
        [ecc, obliq, omega] = Berger_orbpar(year_to_plot+2000); %Fix to J2000.0 convention (add 2000)
        
        %Adjust the obliquity and omegar_bar values to the convention
        obliq = obliq*180/pi;
        omega_bar = omega*(180/pi)-180;
        if omega_bar <0
            omega_bar = 360+omega_bar;
        end
        %prec_for_orbit = 360-omega_bar*180/pi;
        
        handles.data.e = ecc;
        handles.data.precession = omega_bar;
        handles.data.obliquity = obliq;
        
    case 2                          %Laskar solution chosen
        error('Incorrect solution called...')
    otherwise
        msgbox('Invalid solutions method')
end

set(handles.e, 'String', sprintf('%8.7f',handles.data.e));
set(handles.obliquity, 'String', sprintf('%6.4f',handles.data.obliquity));
set(handles.precession, 'String', sprintf('%6.4f',handles.data.precession));


guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function myear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to myear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in solution_menu.
function solution_menu_Callback(hObject, eventdata, handles)
% hObject    handle to solution_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns solution_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from solution_menu

contents = get(hObject,'String');
selection = contents{get(hObject,'Value')};
switch selection
    case 'Orbital Parameters Demo'
        %%%handles.data.cal_mode = 2;
        handles.data.solutions_mode = 0; %demo
        %set(handles.data.myear, 'String', 0);
        %set(handles.data.year_slider,  0;
        set(handles.myear, 'Enable', 'off');
        set(handles.year_slider, 'Enable', 'off');
        set(handles.laskar_year_text, 'Enable', 'off');
        set(handles.laskar_year_slider, 'Enable', 'off');
        
        set(handles.e, 'Enable', 'on');
        set(handles.obliquity, 'Enable', 'on');
        set(handles.precession, 'Enable', 'on');
        
        set(handles.obliquity_slider, 'Enable', 'on');
        set(handles.precession_slider, 'Enable', 'on');
        set(handles.e_slider, 'Enable', 'on');
        
        %Disable Time series options - only allowed with real astro solutions
        set(handles.timeSeries_button, 'Enable', 'off');
        set(handles.start_year_text, 'Enable', 'off');
        set(handles.end_year_text, 'Enable', 'off');
        set(handles.paleo_data_plot_button,'Enable','off');
        
        %Extract demo values for ecc, obliq, and omega_bar
        %and set them to the corresponding text boxes/sliders
        handles.data.e = 0.6; 
        handles.data.obliquity = 45; 
        handles.data.precession = 225;
        
        set(handles.e, 'String', handles.data.e);
        set(handles.obliquity, 'String', handles.data.obliquity);
        set(handles.precession, 'String', handles.data.precession);
        
        %Update the sliders as well: 
        set(handles.e_slider, 'Value', handles.data.e);
        set(handles.obliquity_slider, 'Value', handles.data.obliquity);
        set(handles.precession_slider, 'Value', handles.data.precession);
        
    case 'Berger (1978)'
        handles.data.solutions_mode =1 ; %Berger
        set(handles.myear, 'Enable', 'on');
        set(handles.year_slider, 'Enable', 'on');
        set(handles.laskar_year_text, 'Enable', 'off');
        set(handles.laskar_year_slider, 'Enable', 'off');
        
        set(handles.e, 'Enable', 'off');
        set(handles.obliquity, 'Enable', 'off');
        set(handles.precession, 'Enable', 'off');
        set(handles.obliquity_slider, 'Enable', 'off');
        set(handles.precession_slider, 'Enable', 'off');
        set(handles.e_slider, 'Enable', 'off');
        
        %Set Berger year back to zero
        set(handles.myear, 'String', 0);
        set(handles.year_slider, 'Value', 0);
        
        %Set laskar year back to zero
        set(handles.laskar_year_text, 'String', 0);
        set(handles.laskar_year_slider, 'Value', 0);
        
        %Enable Time series options
        set(handles.timeSeries_button, 'Enable', 'on');
        set(handles.start_year_text, 'Enable', 'on');
        set(handles.end_year_text, 'Enable', 'on');
        set(handles.paleo_data_plot_button,'Enable','on');
        
        %Extract contemporary (year 0) values for ecc, obliq, and omega_bar
        %from Berger and set them to the corresponding text boxes - Adjust Milankovitch
        %text box parameters to the contemporary values based on the Berger solution
        %same algorithm as myear callback
        year_to_plot = 0;
        %Call Berger_orbpar with contemporary year
        [ecc, obliq, omega] = Berger_orbpar(year_to_plot+2000); %Fix to J2000.0 convention (add 2000)
        %Adjust the obliquity and omegar_bar values to the convention
        obliq = obliq*180/pi;
        omega_bar = omega*(180/pi)-180;
        if omega_bar <0
            omega_bar = 360+omega_bar;
        end
        %prec_for_orbit = 360-omega_bar*180/pi;
        
        set(handles.e, 'String', sprintf('%8.7f',ecc));
        set(handles.obliquity, 'String', sprintf('%6.4f',obliq));
        set(handles.precession, 'String', sprintf('%6.4f',omega_bar));
        
        %Update ecc, obliq & precession values in data structure
        handles.data.e = ecc;
        handles.data.precession = omega_bar;
        handles.data.obliquity = obliq;      
    case 'Laskar et al. (2004)'
        handles.data.solutions_mode = 2; %Laskar
        set(handles.myear, 'Enable', 'off');
        set(handles.year_slider, 'Enable', 'off');
        set(handles.laskar_year_text, 'Enable', 'on');
        set(handles.laskar_year_slider, 'Enable', 'on');
        
        %Update text boxes and sliders
        set(handles.e, 'Enable', 'off');
        set(handles.obliquity, 'Enable', 'off');
        set(handles.precession, 'Enable', 'off');
        set(handles.obliquity_slider, 'Enable', 'off');
        set(handles.precession_slider, 'Enable', 'off');
        set(handles.e_slider, 'Enable', 'off');
        
        %Set Laskar year back to zero
        set(handles.laskar_year_text, 'String', 0);
        set(handles.laskar_year_slider, 'Value', 0);
        
        %Set Berger year back to zero
        set(handles.myear, 'String', 0);
        set(handles.year_slider, 'Value', 0);
        
        %Enable Time series options
        set(handles.timeSeries_button, 'Enable', 'on');
        set(handles.start_year_text, 'Enable', 'on');
        set(handles.end_year_text, 'Enable', 'on');
        set(handles.paleo_data_plot_button,'Enable','on');
        
        %Adjust Milankovitch text box and slider parameters to the contemporary values
        %base on the Laskar solution
        year_to_plot = 0;
        %Call getLaskar with year 0
        [ecc, obliq, omega_bar] = getLaskar(year_to_plot, handles.data.laskar_neg, handles.data.laskar_pos);
        
        set(handles.e, 'String', sprintf('%8.7f',ecc));
        set(handles.obliquity, 'String', sprintf('%6.4f',obliq));
        set(handles.precession, 'String', sprintf('%6.4f',omega_bar));
        
        %set(handles.start_year_text, 'String', handles.data.start_year_text);
        %set(handles.end_year_text, 'String', handles.data.end_year_text);
        
        %Update ecc, oblilq & precession values in data structure
        handles.data.e = ecc;
        handles.data.precession = omega_bar;
        handles.data.obliquity = obliq;     
    otherwise
        error('');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function solution_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solution_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in help_button.
function help_button_Callback(hObject, eventdata, handles)
% hObject    handle to help_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wm = warndlg('The ReadMe.txt in the model root directory should open in ~6 s in your default editor. If not, open it manually to view the help.');

pause(6);
if ishandle(wm)
    close(wm);
end
edit('ReadMe.txt');

% --- Executes on button press in timeSeries_button. - Plots LASKAR's
% solution only - need to make a note of this in the GUI or help
function timeSeries_button_Callback(hObject, eventdata, handles)
% hObject    handle to timeSeries_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the start and end years from the corresonding time series text boxes

st_year = handles.data.start_year_text;
en_year = handles.data.end_year_text;

%figure out if Berger or Laskar solution should be used
switch handles.data.solutions_mode
    case 0                       %Demo solution chosen
        msgbox('Error - this functionality cannot be used in demo mode'); %This should never appear
    case 1                       %Use berger solution
        
        %Berger_data.dat = output data from Tiho's Berger78_driver function
        %output = year; ecc; obliq; omega
        Berger_solution = load('Berger78_data.dat');
        %Year values
        table_year = Berger_solution(:,1)/1000 - 2;
        %IMPORTANT!!! NEED to offset Berger years by 2,000 years because he has calendar years,
        %Laskar has years since J2000.0
        
        %eccentricity values
        ecc = Berger_solution(:,2);
        %obliquity values
        obliq = Berger_solution(:,3);
        %precession - need to adjust
        precess = Berger_solution(:,4)-180; %omega_bar
        
        pn = precess <0;
        precess(pn) = precess(pn) + 360;
        
    case 2                       %use Laskar solution
        %THIS CODE TAKEN AND AMENDED FROM GETLASKAR
        laskar_neg = handles.data.laskar_neg;
        laskar_pos = handles.data.laskar_pos;
        %Remove 0 row from negative file
        laskar_neg(1,:) = [];
        %Flip the negative file
        laskar_neg = flipud(laskar_neg);
        %Concatenate the files
        Laskar_solution = vertcat(laskar_neg, laskar_pos);
        %Extract each column and assign variables
        %Year values
        table_year = Laskar_solution(:,1);
        %Eccentricity values - no adjustment needed
        ecc = Laskar_solution(:,2);
        %obliquity values - no asjustment needed
        obliq = (180/pi)*Laskar_solution(:,3);
        %Need to adjust precession angle - find where it jumps from 360 to 0
        precess = (180/pi)*Laskar_solution(:,4); %This is omega_bar
        
    otherwise
        msgbox('Invalid solutions method')
end

%%%%%PREPARE TO COMPUTE INSOLATION TIME SERIES %%%%%%%%%%%
wm = warndlg('Please wait for the computation to complete.  This may take from a few moments to more than a minute, depending on the length of your time series and the speed of your computer...');
q = find(table_year<=st_year,1,'last');
q2 = find(table_year>=en_year,1,'first');

%Truncate table years and Milankovitch vectors to only the requested time interval
table_year = table_year(q:q2);
ecc = ecc(q:q2);
obliq = obliq(q:q2);
precess = precess(q:q2);

%Plot in steps of 1,000 years, otherwise it is too slow:

%Alternatively, or in additon, INSIDE insol_TS, if time series is too long,
%we can down sample more to up to 5,000 years steps.
if handles.data.solutions_mode==1 %Berger
    %Years are in units of thousands of years, select whole thousands from the tables only, with the end and beginning always selected
    %q3 = abs(rem(table_year,1)) < 1e-6; %Exact thousand years (integer)
    bind = 1:10:length(table_year); %get every 10th record from Berger (every 1000 years)
    table_year = [table_year(bind); table_year(end)]; %include last record for completeness
    ecc = [ecc(bind); ecc(end)];
    obliq = [obliq(bind); obliq(end)];
    precess = [precess(bind); precess(end)];
end

%The following snippet (4 code lines) is taken from the insolation snapshot button
%callback:
prec_for_orbit = 180-precess; %The one used internally by keplerian and orbit.m
if prec_for_orbit < 0 %If this is not implemented, prec_for_orbit becomes negatives and that
    %messes up some things within orbit.m, particularly the season length calculations
    prec_for_orbit = 360 + prec_for_orbit;
end

dayofyear = 1:365;
[sol, solstice_sol, annual_mean_sol] = insol_TS(handles.data.sm_axis,handles.data.AU,...
    handles.data.period,table_year,dayofyear, ecc,obliq,prec_for_orbit,handles.data.Fo,...
    handles.data.latitude,handles.data.month,handles.data.day);

if ishandle(wm)
    close(wm);
end

f3h = figure(3);
clf;
set(f3h,'Name','Milankovitch Parameters Time Series','NumberTitle','off','Units','normalized','OuterPosition',[.1 .1  .7 .85]);
subplot(3,1,1)
plot(table_year, ecc, 'b-')
title('Eccentricity')
axis([st_year en_year 0 0.06])
ylabel('dimensionless')
hold on
plot([0 0],[0 0.06],'m-')

%obliquity in green
subplot(3,1,2)
plot(table_year, obliq, 'g-')
title('Obliquity')
axis([st_year en_year 21.5 25])
ylabel('degrees')
hold on
plot([0 0],[21.5 25],'m-')
%precession in red

%Use precess (omega_bar, longitude of perihelion) to reconstruct omega,
%longitude of perigee, the quantity needed to reconstruct climatic precession:
omega = mod(precess + 180,360); %If greater than 360 degrees, get just remainder...

subplot(3,1,3)
[ax,h1,h2] = plotyy(table_year, precess,table_year,ecc.*sin(omega*(pi/180)));

title('Longitude of perihelion/Climatic Precession')
set(get(ax(1),'Ylabel'),'String','Long. of perihelion, \omega_t_i_l_d_e, deg.');
set(get(ax(2),'Ylabel'),'String','Climatic precession, e*sin(\omega)','Color','red');
set(h1,'LineStyle','-', 'Color',[.6 .6 .6])
set(h2,'LineStyle','-', 'Color','r')
axis(ax(1),[st_year en_year -5 365]);
axis(ax(2),[st_year en_year -0.051 0.051]);

set(ax(2),'YColor','red')
set(ax(1),'YColor','k', 'YTick',[0 90 180 270 360])
set(ax(2),'YColor','r', 'YTick',[-.05 -0.025 0 .025 .05])
%Plot some interesting paleoclimatological data here...
%Test for plotyy
box off
hold on
plot(ax(1),[0 0],[-5 365],'m-')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  PLOT TIME SERIES OF INSOLATION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mo_string = '';
switch handles.data.month
    case  1
        mo_string = 'January';
    case 2
        mo_string = 'February';
    case 3
        mo_string = 'March';
    case 4
        mo_string = 'April';
    case 5
        mo_string = 'May';
    case 6
        mo_string = 'June';
    case 7
        mo_string = 'July';
    case 8
        mo_string = 'August';
    case 9
        mo_string = 'September';
    case 10
        mo_string = 'October';
    case 11
        mo_string = 'November';
    case 12
        mo_string = 'December';
    otherwise
        error('');
end

f4h = figure(4);
clf;
set(f4h,'Name','Insolation Time Series at Specified Date and Latitude','NumberTitle','off','Units','normalized','OuterPosition',[.15 .3  .7 .6]); %Fig 3 is [.1 .1  .7 .85]);
[ax, h1, h2] = plotyy(table_year, solstice_sol,table_year,annual_mean_sol);
set(get(ax(1),'Ylabel'),'String',['Insolation for ', mo_string,' ',num2str(handles.data.day), ', W m^-^2']);
set(get(ax(2),'Ylabel'),'String','Annual mean of daily insolation, W m^-^2');
set(ax(2),'YColor','r');
set(h2,'Color','r');
xlabel('Thousand of years since J2000');
title(['Insolation at ',num2str(handles.data.latitude),'^o latitude, W m^-^2']);
rng = max(solstice_sol)-min(solstice_sol);
if rng<=0
    rng = 10; %arbitrarily assign
end
rng2 = max(annual_mean_sol)-min(annual_mean_sol);
if rng2<=0
    rng2 = 10; %arbitrarily assign
end
y1lim = [min(solstice_sol)-0.1*rng, max(solstice_sol)+0.1*rng];
y2lim = [min(annual_mean_sol)-0.1*rng2, max(annual_mean_sol)+0.1*rng2];
axis(ax(1),[st_year en_year y1lim]);
axis(ax(2),[st_year en_year y2lim]);
s = axis;
hold on
plot(ax(1),[0 0],[s(3) s(4)],'m-')
box off
set(ax(1),'YTickMode', 'auto');
set(ax(1),'YTickLabelMode', 'auto');
set(ax(2),'YTickMode', 'auto');
set(ax(2),'YTickLabelMode', 'auto');

f5h = figure(5);
clf;
set(f5h,'Name','3D Insolation Time Series','NumberTitle','off','Units','normalized','OuterPosition',[.2 .25  .5 .6])
surf(table_year,dayofyear,sol);
view(2);
axis([-Inf Inf 1 365]);
axis ij
shading interp
colorbar
xlabel('Thousands of years since J2000');
ylabel('Day of year');
title(['Daily insolation at ',num2str(handles.data.latitude),'^o latitude averaged over 24 hrs, W m^-^2']);

if handles.data.save_insol_data
    [ofile, opath] = uiputfile('orbit_output/*.dat', 'Save Insolation Time Series data as');
    if ofile
        tmp = [NaN dayofyear];
        tmp = tmp';
        out = [table_year'; sol];
        out = [tmp out];
        save('-ascii',[opath, ofile],'out');
    end
end

guidata(hObject, handles);


% --- Executes on button press in insolation_button.
function insolation_button_Callback(hObject, eventdata, handles)
% hObject    handle to insolation_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prec_for_orbit = 180-handles.data.precession; %The one used internally by keplerian and orbit.m
if prec_for_orbit < 0 %If this is not implemented, prec_for_orbit becomes negatives and that
    %messes up some things within orbit.m, particularly the season length calculations
    prec_for_orbit = 360 + prec_for_orbit;
end

wm = warndlg('Please wait for the computation to complete. This can take a few moments...');
dayofyear = [1:5:365, 365];
lats = 90:-5:-90;
[sol, ~] = insol_3d(handles.data.sm_axis,handles.data.AU,handles.data.period,handles.data.e,handles.data.obliquity,prec_for_orbit,handles.data.Fo,dayofyear,lats);
if ishandle(wm)
    close(wm);
end

f2h = figure(2);
clf;
set(f2h,'Name','3D Global Insolation in a Specified Year','NumberTitle','off','Units','normalized','OuterPosition',[0.05 0.25 .7 .75]);

subplot(2,1,1)
surf(dayofyear,lats,sol)
view(2)
axis([1 365 -90 90])
shading interp
colorbar
xlabel('Day of year')
ylabel('Latitude, degrees')
title('Daily insolation averaged over 24 hrs, W m^-^2')

contemporary_sol = load('contemp_insol.dat'); %file containing contemporary (J2000.0) insolation pattern...
sol_anomaly = sol - contemporary_sol;

subplot(2,1,2)
surf(dayofyear,lats,sol_anomaly)
view(2)
axis([1 365 -90 90])
shading interp
colorbar
xlabel('Day of year')
ylabel('Latitude, degrees')
title('Daily insolation anomaly, W m^-^2')

if handles.data.save_insol_data
    [ofile, opath] = uiputfile('orbit_output/*.dat', 'Save Insolation data as');
    if ofile
        out = [dayofyear; sol];
        tmp = [NaN lats];
        tmp = tmp';
        out = [tmp out];
        save('-ascii',[opath, ofile],'out');
    end
    
    [ofile, opath] = uiputfile('orbit_output/*.dat', 'Save Insolation Anomalies data as');
    if ofile
        out = [dayofyear; sol_anomaly];
        tmp = [NaN lats];
        tmp = tmp';
        out = [tmp out];
        save('-ascii',[opath,ofile],'out');
    end
end


function start_year_text_Callback(hObject, eventdata, handles)
% hObject    handle to start_year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_year_text as text
%        str2double(get(hObject,'String')) returns contents of start_year_text as a double
e = str2double(get(hObject, 'String'));
% if isnan(e) | ~(e>=get(handles.start_year_text,'Min') & e<= get(handles.obliquity_slider,'Max'))
%     e = handles.data.obliquity; %assume the old value
%     set(hObject, 'String', e);
%     errordlg(['Input must be a number b/n', num2str(get(handles.obliquity_slider,'Min')), ' and ' ,...
%         num2str(get(handles.obliquity_slider,'Max'))],'Error');
% end
% Save the new e value
handles.data.start_year_text = e;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function start_year_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function end_year_text_Callback(hObject, eventdata, handles)
% hObject    handle to end_year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_year_text as text
%        str2double(get(hObject,'String')) returns contents of end_year_text as a double
e = str2double(get(hObject, 'String'));

handles.data.end_year_text = e;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function end_year_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function laskar_year_slider_Callback(hObject, eventdata, handles)
% hObject    handle to laskar_year_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

a = get(hObject, 'Value');
%update value of corresponding field
set(handles.laskar_year_text, 'String', a);
%update value of e;
handles.data.laskar_year_text = a;


%update Milankovitch parameters
%calculate

switch handles.data.solutions_mode
    case 0                       %Demo solution chosen
        %msgbox('Demo method')
        %No need to do anything here...
        error('Calling year slider when demo is chosen is incorrect...')
        
    case 1                       %Berger solution chosen
        %msgbox('Berger method')
        error('Wrong solution slider chosen...')
        
    case 2                          %Laskar solution chosen
        %year_to_plot = handles.data.laskar_year_text;
        %         %msgbox('Laskar method')
        %         %Implement Laskar method
        %
        ecc = -999;     %*ones(size(t));
        obliq = ecc;
        omega_bar = ecc;
        
        %Call getLaskar method to obtain Milank values for the input year
        %Divide year by 1000 b/c Laskar gives his values every 1000 years
        [ecc, obliq, omega_bar] = getLaskar(handles.data.laskar_year_text, handles.data.laskar_neg, handles.data.laskar_pos);
        
        handles.data.e = ecc;
        handles.data.precession = omega_bar;
        handles.data.obliquity = obliq;
    otherwise
        msgbox('Invalid solutions method')
end

set(handles.e, 'String', sprintf('%8.7f',handles.data.e));
set(handles.obliquity, 'String', sprintf('%6.4f',handles.data.obliquity));
set(handles.precession, 'String', sprintf('%6.4f',handles.data.precession));

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function laskar_year_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laskar_year_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function laskar_year_text_Callback(hObject, eventdata, handles)
% hObject    handle to laskar_year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of laskar_year_text as text
%        str2double(get(hObject,'String')) returns contents of laskar_year_text as a double

e = str2double(get(hObject, 'String'));
if isnan(e) || ~(e>=get(handles.laskar_year_slider,'Min') && e<= get(handles.laskar_year_slider,'Max'))
    e = handles.data.laskar_year_text; %assume the old value
    set(hObject, 'String', e);
    errordlg(['Input must be a number b/n ', num2str(get(handles.laskar_year_slider,'Min')), ' and ' ,...
        num2str(get(handles.laskar_year_slider,'Max'))],'Error');
end
% Save the new e value
handles.data.laskar_year_text = e;
%update slider
set(handles.laskar_year_slider, 'Value',e);

switch handles.data.solutions_mode
    case 0                       %Demo solution chosen
        %msgbox('Demo method')
        %No need to do anything here...
        error('Calling year slider when demo is chosen is incorrect...')
        
    case 1                       %Berger solution chosen
        error('Wrong solution referenced...')
        
    case 2                          %Laskar solution chosen
        %year_to_plot = handles.data.laskar_year_text;
        %         %msgbox('Laskar method')
        %         %Implement Laskar method
        %
        ecc = -999;     %*ones(size(t));
        obliq = ecc;
        omega_bar = ecc;
        
        %Call getLaskar method to obtain Milank values for the input year
        %Divide year by 1000 b/c Laskar gives his values every 1000 years
        [ecc, obliq, omega_bar] = getLaskar(handles.data.laskar_year_text, handles.data.laskar_neg, handles.data.laskar_pos);
        
        handles.data.e = ecc;
        handles.data.precession = omega_bar;
        handles.data.obliquity = obliq;
        
    otherwise
        msgbox('Invalid solutions method')
end

set(handles.e, 'String', sprintf('%8.7f',handles.data.e));
set(handles.obliquity, 'String', sprintf('%6.4f',handles.data.obliquity));
set(handles.precession, 'String', sprintf('%6.4f',handles.data.precession));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function laskar_year_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laskar_year_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on help_button and none of its controls.
function help_button_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to help_button (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saving_data_checkbox.
function saving_data_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to saving_data_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')==get(hObject,'Max')
    handles.data.save_insol_data = 1;
else
    handles.data.save_insol_data = 0;
end
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of saving_data_checkbox


% --- Executes on button press in paleo_data_plot_button.
function paleo_data_plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to paleo_data_plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of paleo_data_plot_button
% if get(hObject,'Value')==get(hObject,'Max')
%     handles.data.plot_paleo_data = 1;
% else
%     handles.data.plot_paleo_data = 0;
% end
st_year = handles.data.start_year_text;
en_year = handles.data.end_year_text;

%Load EPICA CO2 data
co2 = load('EPICA_CO2.dat');
co2(:,1) = -co2(:,1);

%Load EPICA deuterium/temperature data
deut = load('EPICA_deuterium.dat');
deut(:,3) = -(deut(:,3)+50);

%Load the Lisiecki and Raymo 2005 benthic d-O-18 stack data: 
LR04 = load('LR04_dO18_benthic_stack.dat');

%Load the Zachos et al. (2001) benthic d-O-18 data set: 
Zachos = load('Zachos01_d18O.dat'); 

f6h = figure(6);
clf;
set(f6h,'Name','Paleoclimatological Data Time Series','NumberTitle','off','Units','normalized','OuterPosition',[.11 .11  .7 .85]);

subplot(2,1,1)
[ax,h1,h2] = plotyy(co2(:,1)/1000,co2(:,2),deut(:,3)/1000,deut(:,5));
axis(ax(1),[st_year en_year 180 295]);
axis(ax(2),[st_year en_year -11 5.5]);
title('EPICA CO_2 and temperature');
set(get(ax(1),'Ylabel'),'String','[CO_2], ppmv');
set(get(ax(2),'Ylabel'),'String','Temperature, \circC');
set(ax(1),'YTick',[180:20:300]);
set(ax(2),'YTick',[-10 -7.5 -5 -2.5 0 2.5 5])
xlabel('Thousands of years since J2000');
hold on
plot(ax(1),[0 0],[180 295],'m-');
set(ax(1),'box','off');

subplot(2,1,2)
q = -LR04(:,1)>=st_year & -LR04(:,1)<=en_year; 
q2 = -Zachos(:,1)>=st_year/1000 & -Zachos(:,1)<=en_year/1000;
plot(-LR04(q,1),LR04(q,2),'b-');
hold on
plot(-Zachos(q2,1)*1000, Zachos(q2,2),'r-');
ymin = 0.95*min([LR04(q,2); Zachos(q2,2)]);
ymax = 1.05*max([LR04(q,2); Zachos(q2,2)]);
axis([st_year en_year ymin ymax]);
set(gca,'YDir','reverse');
xlabel('Thousands of years since J2000');
ylabel(['Benthic Foraminiferal \delta^1^8O, ', char(8240)]); 
legend({'Lisiecki and Raymo 2005 \delta^1^8O','Zachos et al. (2001) \delta^1^8O'},'Location','NorthOutside','Orientation','Horizontal');
plot([0 0],[ymin ymax],'m-')

guidata(hObject,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over timeSeries_button.
function timeSeries_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to timeSeries_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over paleo_data_plot_button.
function paleo_data_plot_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to paleo_data_plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
