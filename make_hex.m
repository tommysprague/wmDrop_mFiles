% make_hex
%
% requires an odd integer, makes a heaxagonal grid w/ that integer number
% of rows/that integer wide
%
%


function [x, y] = make_hex(n)


if mod(n,2)==0
    x = NaN;
    y = NaN;
    error('hexMap:invalidInput','Input must be an odd integer');
end

step_size = 2/(n-1);  % normalized so horizontal axis is longest, and reaches +/- 1

y_per_row = linspace(-1,1,n)*sqrt(3)/2;

points_per_row =  [(ceil(n/2):n) ((n-1):-1:ceil(n/2))];
x = NaN(sum(points_per_row),1);
y = NaN(sum(points_per_row),1);

startidx = 1;
for rr = 1:length(points_per_row)
    
        
        thisx = step_size*(1:points_per_row(rr));
        thisx = thisx - mean(thisx);
        thisy = y_per_row(rr)*ones(length(thisx),1);
        
        x(startidx:(startidx+points_per_row(rr)-1)) = thisx;
        y(startidx:(startidx+points_per_row(rr)-1)) = thisy;
    
    clear thisx thisy;
    startidx = startidx+points_per_row(rr);
    
end


return
