function [r_vector] = rotate_vector(vector,angle)
%R Summary of this function goes here
%   Detailed explanation goes here

RotM = [...
    cos(angle), -sin(angle, 0); ...
    sin(angle), cos(angle, 0); ...
    0,  0,  1];

r_vector = RotM * vector;
end