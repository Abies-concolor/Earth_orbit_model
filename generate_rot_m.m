function tilt_m = generate_rot_m(obliquity, precession)

%Construction of rotation matrix to rotate Earth around dir_vector by
%obliquity degrees. This is an affine transformation. 

%The code below uses ideas from Wikipedia's page on Rotation matrix:
%http://en.wikipedia.org/wiki/Rotation_matrix

%dir_vector = [cos(precession) sin(precession) 0 ]; 
% u,v,w components directly coded below
%Axis of rotation is in the orbital plane (w-component = 0), Earth to be
%rotated by an angle equal to obliquity

%This is a utility function expected to be called only internally by the
%orbit model, and no input checking/warnings are performed. 

%Dr. T.Kostadinov

cos_eps = cos(-obliquity);
sin_eps = sin(-obliquity);

onem_coseps = 1 - cos_eps;

u = cos(precession);
v = sin(precession);
w = 0;

tilt_m = [cos_eps+u^2*onem_coseps u*v*onem_coseps-w*sin_eps u*w*onem_coseps+v*sin_eps; ...
    u*v*onem_coseps+w*sin_eps cos_eps+v^2*onem_coseps v*w*onem_coseps-u*sin_eps; ...
    u*w*onem_coseps-v*sin_eps v*w*onem_coseps+u*sin_eps cos_eps+w^2*onem_coseps]';