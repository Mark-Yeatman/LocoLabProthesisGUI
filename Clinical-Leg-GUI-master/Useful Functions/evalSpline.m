function [ q ] = evalSpline( x, xvals, yvals, kvals)
%EVALSPLINE Copied and modifeid from http://blog.ivank.net/interpolation-with-cubic-splines.html
%xvals, yvals, and kvals must be 1xn arrays
%is a 1xM array
%x inputs are sorted low to high
    % find indices of enclosing knots 
    indices = zeros(1,length(x));
    q = zeros(1,length(x));
    for j=1:length(x)
        if(x(j)==xvals(1))
            indices(j) = 2;
        else
            temp = find(xvals>=x(j));
            if(isempty(temp))
                error('x input value outside of range of xvals inputs')
            end
            indices(j) = temp(1); 
        end
    end

    %evaluate points
    for j=1:length(x)
    i= indices(j);
    t = ( x(j)- xvals(i-1) ) / ( xvals(i) - xvals(i-1) );
		
    a =  kvals(i-1) * ( xvals(i) - xvals(i-1) ) - ( yvals(i) - yvals(i-1) );
    b = -kvals(i)*(xvals(i)-xvals(i-1)) + (yvals(i)-yvals(i-1));
		
    q(j) = (1-t)*yvals(i-1) + t*yvals(i) + t*(1-t)*(a*(1-t)+b*t);
    end
end


