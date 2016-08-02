function fwhm = rad2fwhm(rad,pow)

if nargin < 2
    pow = 7;
end

fwhm = 2*rad*acos((0.5^(1/pow)-0.5)/0.5)/pi;

return