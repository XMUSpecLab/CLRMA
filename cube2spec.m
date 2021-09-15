function spectrum = cube2spec(cube,xx,yy,zz)
% this function change the 3d data to one dimensional spectrum
% if xx yy are or one of them is range variable, the result return the average
% spectrum in the region
% if xx yy are number, the result return the single spectrum indicated by
% xx and yy
% zz is the selected range of the spectrum
if nargin == 3
    zz = 1:size(cube,3);
end

if length(xx) == 1 && length(yy) == 1
    spectrum = squeeze(cube(xx,yy,zz));
else 
    spectrum = mean(mean(cube(xx,yy,zz),1),2);
    spectrum = squeeze(spectrum);
end
end