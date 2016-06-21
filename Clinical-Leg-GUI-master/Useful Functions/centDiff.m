function [ d ] = centDiff( points, xvals, yvals )
%CENTDIFF Uses the central difference method to estimate derivative of a function given discrete sets of inputs and outputs. 

%points, Nx1 array of input points to find derivative off
%xvals, Mx1 function input values
%yvals, Mx1 function output values
    if( any(size(xvals)~= size(yvals)) || size(xvals,2)~=1 )
        error('dimension function inputs and outputs arrays are incorrect');
    end
    d = zeros(1,length(points));
    for i=1:length(points)
        index = find(xvals == points(i),1);
        d(i) = ( yvals(index+1)-yvals(index-1) ) / (  xvals(index+1)-xvals(index-1) );
    end

end

