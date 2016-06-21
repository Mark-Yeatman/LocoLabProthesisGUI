function [p,t,info ] = grad7( d,deg,stop )
% GRAD7.M This function takes a given set of ordered data and returns
% the parameter values (nodes) and control points which determine the
% Bezier curve that fits the data in the total least squares sense.
% Instead of minimizing the vertical distance to the curve we minimize
% both vertical and horizontal distance. The central feature of this
% function is the use of the Gauss-Newton method to estimate the
% nearest point along the Bezier curve from each data point.

% Copied from BEZIER CURVE FITTING by Tim Andrew Pastva ,September 1998

% Input Arguments:
%   d       (i x 2) matrix of ordered data points
%   deg     degree of the curve the user wants to fit to the data
%   stop    stopping criteria number which in the algorithm becomes 10^(stop)

% Output Arguments:
%   P       control points for best fit Bezier curve
%   t       nodes for best fit Bezier curve
%   info    (2 x 1) vector which has the final norm of the residual 
%           and the number of iterations to convergence
% 
% Basic Algorithm:
% 
% 1. Determine and plot the best fit Bezier curve by
% solving a linear least squares problem using an
% initial set of nodes based on the
% 'spread' and 'bend' of the data.

% 2. Perform the following until the stopping
% criteria is met:
% 
%   a. Determine a new set of nodes
%   which estimates points on the current Bezier
%   points.
% 
%   b. Determine the control points for the Bezier
%   curve which will best fit the data using the
%   new set of nodes.

% ALGORITHM STEP 1.
%, Check the hold state so it can be returned to how the user had it.
    if (ishold)
        hold_was_off =0;
    else
        hold_was_off = 1;
    end
    
% 'i' is the number of data points and 'j' is the number of control
% points required for the degree of the curve specified by the user.
% aff_angle(d) returns the initial set of nodes, a vector of 'i'
% elements, based on the 'spread' and 'bend' of the data.

    i = size(d,1);
    j = deg+1;
    t = aff_angle(d);

% 'bez_mat' is a (i x j) matrix which is determined by the nodes
% and the degree of the curve desired. Since we would
% like to fit the data exactly, we should solve:
%           d = bez_mat * p

% for the desired (j x 2) matrix of control points 'p' . This is a
% linear least squares problem since the nodes are known. 'p'
% is determined by:
%        p = pinv(bez_mat) * d

% which is the same as using the matlab 'backslash' command.

    bez_mat = mxbern2(t,deg);
    p = bez_mat \ d;

% Once we have the control points 'p' , we can determine the Bezier
% curve and the points on the curve associated with the initial vector
% t. 'y'is the (i x 2) matrix of points on the Bezier curve
% associated with the initial vector t. 'tl' is a closely spaced
% set of parameter values which will produce the Bezier matrix,
% 'bez_mat_l', which will give enough points on the fitted Bezier curve
% for matlab to plot a smooth looking curve: 'yl' is the matrix of
% these points.

    y = bez_mat * p;
    t1 = (0:1/128:1)';
    bez_mat_l = mxbern2(t1,deg);
    y1 = bez_mat_l * p;

% Now we plot the results for the user. A legend is provided and the
% axes are made 'equal' to eliminate distortions.
% 'Pause' will keep the plot displayed and delay this function until
% the user presses any key on the keyboard.

    figure
    hold on
    plot(y1(:,1),y1(:,2),p(:,1),p(:,2),'*',d(:,1),d(:,2))
    %plot(p(:,1),p(:,2),'*')
    %plot(y(:,1),y(:,2),'o')
    %plot(d(:,1),d(:,2),'+')
    %axis('equal')
    legend('Fitted Curve','Control Points','Data Points')
    pause

% ALGORITHM STEP 2.
% We will perform steps 2a and 2b within a 'while' loop with stopping
% criteria to meet in order to end the loop. Our stopping criteria
% is based on the relative change of the residual. We also initialize
% the iteration counter to zero. The 'tic' command starts a clock so
% that the user will know how much computer time was required to
% meet the stopping criteria.

    iter = 0;
    resid_old = 0;
    resid_new = bez_mat * p - d;
    
    tic
    
    while norm(resid_new - resid_old)/max(1 , norm(resid_new)) > 10^(stop)
%   Algorithm step 2a. Each iteration of the 'while' loop produces a
%   new vector t by solving the nonlinear least squares problem
%           B(t) * p = d

%   where 'p' is the vector of control points, 'd' the matrix of data
%   points, and 'B(t)' is the Bernstein matrix with unknown parameter
%   values. Matrix 'B(t)' is nonlinear in terms of the parameter
%   values. The method used in this function is the Gauss-Newton method
%   where we let the residual be R(t) = B(t) * p - d and we let the
%   Jacobian matrix 'J' be such that J_i,j = dB(t_i)/dt_j . I.E., the
%   (i,j)"th element of 'J' is the slope along the Bezier curve at the
%   i'th parameter value with respect to the j~th parameter
%   value. The Gauss-Newton method says that the change in parameter
%   values which will minimize the residual is given by:
%           delta.t = -inv(J' * J) * J' * R

%   where 'J' and 'R' are evaluated at the current parameter values.

%   To compute the gradient at a point on the Bezier curve, we need the
%   forward difference of the control points. I.E., we need a
%   (j-1 x 2) matrix where the entries are p_i+l - p_i for i=l,..,j-l.
%   The slope, 'deriv', is then determined multiplying the Bernstein
%   matrix for a degree-minus-one curve by the forward difference matrix
%   of control points and then multiplying by the degree of the original
%   Bezier curve.

        deriv = deg * mxbern2(t,deg-1) * (p(2:j,:) - p(1:j-1,:));
        
%   y, Now we have what we need for the Gauss-Newton method. Since to use
%   this method the residual needs to be a vector, we simply take the
%   y-values of the residual and append them to the bottom of the
%   x-values. This is done using 'resid(:)'. Similarly, the Jacobian
%   matrix's elements must be for the new vector 'resid'. Each element
%   of the matrix 'J' is the x-value or y-value of the slope at each
%   parameter's point. Since the slope at any point on the curve
%   doesn't change with any parameter other than it's own, most of the
%   entries in 'J' are zero. 

%   't', the new nodes are given by 't = t - delta_t'
%   using the formula above. Because (J' * J) and J' have a form we
%   know in advance, we can form 'delta_t' using less computer time and flops.

        t = t - (deriv(:,1).*resid_new(:,1) + deriv(:,2).*resid_new(:,2)) ...
        ./ (deriv(:,1).^2 + deriv(:,2).^2) ;
    
%   Now we have a new vector t, but we want to make sure the values
%   are between 0 and 1. In most cases with well behaved data, the
%   following rescaling of the nodes also results in the endpoints
%   being associated with the nodes t_l=0 and t_m=l.

        t = -min(t)*ones(i,1) + t;
        t = t/max(t);
        
%   With ordered nodes we now want the new control points so
%   that we can reproduce the points 'tau' on the curve for the
%   next iteration of the 'while' loop or to be used in the final plot
%   if the stopping criteria is met. Note that if the condition is met
%   we can also use the below 'bez_mat' matrix for the final plot. We
%   also compute a new value for 'resid_new' and update the iteration
%   counter.       
        
        bez_mat = mxbern2(t,deg);
        p = bez_mat \ d;
        
        resid_old = resid_new;
        resid_new = bez_mat * p - d;
        
        iter = iter + 1; 

%   End 'while' loop and stop computer time clock to show user how long
%   it took to converge.
    end
    
    toc

% Now that we have the best fit vector t and the associated
% matrix 'p' of control points, we are ready to plot the final
% solution for the user, 'y' are the points on the Bezier curve
% associated with the vector t. 'yl' are the closely spaced
% points on the Bezier curve which matlab will use to plot a smooth
% looking curve.

    y = bez_mat * p;
    y1 = bez_mat_l * p;

% Plot final data 
    figure
    hold on
    plot(y1(:,1),y1(:,2),p(:,1),p(:,2),'*',d(:,1),d(:,2))
    %plot(p(:,1),p(:,2),'*')
    %plot(y(:,1),y(:,2),'o')
    %plot(d(:,1),d(:,2),'+')
    %axis('equal')
    legend('Fitted Curve','Control Points','Data Points')
    pause
% Return hold state to however the user had it.
    if (hold_was_off)
        hold off;
    end
    
    info = [norm(bez_mat*p-d), iter];
% End GRAD7.M    
end

function [B] = mxbern2(t,d)
% MXBERN2.M This function creates a Bernstein matrix of degree d
% using the values in the column vector t. A Bernstein matrix
% is a generalized Vandermonde matrix whose (i,j) entry is
% B_{j-l}~d(t_i) . The vector t must be a column vector with values
% between 0 and 1. Copyright (c) 3 December 1994 by Carlos F. Borges.
% All rights reserved. Modified by permission for this paper.
    [n,m] = size(t);
    
    if (m ~= 1)
        error('t must be a column vector.');
    end
% Check to see if nodes are within [0,1].
    if min(t) < 0 || max(t) > 1
        error('nodes are not within [0,1]')
    end
% Build the Bernstein matrix.
    ct = 1 - t;
    B = zeros(n,d+1);
    
    for i = 0 : d
        B(:,i+1) = (t.^i).*(ct.^(d-i));
    end
% If d > 22 we form the Bernstein matrix differently to
% avoid roundoff.
    if d < 23
        B = B*diag( [1 cumprod(d:-1:1)./cumprod(1:d)] );
    else
        B = B*diag(diag(fliplr(pascal(d+1))));
    end
% End MXBERM2.M
end

function [h] = aff_angle(X)
% AFF.ANGLE.M This function returns an affine invariant vector of
% nodes for a given set of ordered data X. It is
%assumed that the user sends this function the data arranged in a
% (n x 2) matrix ordered from row one to row n.
% Obtaining the covariance matrix A is from a function AFF_KNT.M which
% is copyrighted on 3 December 1994 by Carlos F. Borges and appears
% here with his approval. The affine invariant angle method of
% obtaining nodes is found in a paper by Thomas A. Foley
% and Gregory M. Nielson, "Knot Selection for Parametric Spline
% Interpolation", in the book "Mathematical Methods in Computer Aided
% Geometric Design", Academic Press, Inc., 1989. Notation from this
% paper is used to annotate the steps in this algorithm.
    n = size(X,1);
    Xbar = X - ones(size(X)) * diag(mean(X));
    Xcov = Xbar' * Xbar/n;
    A = inv(Xcov);
% Obtain the node spacing values using the metric
%   t_i = M[X](X_i,X_(i+l)).

    V = X(2:n,:) - X(1:n-1,:);
    t = diag(V * A * V') .^ (1/2);
% Obtain the values for
% M-2[X](X_(i-l),X(i+l))

    V_2 = X(3:n,:) - X(1:n-2,:);
    t_2 = diag(V_2 * A * V_2');

% Get theta_i values. This is what takes into account the bending
% of the data.

    theta = zeros(n-1,1);
    
    for j = 2:n-1
        theta(j) = min( pi - acos( (t(j-1)^2 + t(j)^2 - ...
            t_2(j-1))/(2*t(j)*t(j-1)) ) , pi/2);
    end
    
% Obtain the affine invariant angle node spacing values h_i.

    h = zeros(n-1,1);
    h(1) = t(1)*( 1 + (1.5*theta(2)*t(2))/(t(1) + t(2)) );
    
    for j = 2:n-2
        h(j) = t(j) * ( 1 + (1.5*theta(j)*t(j-1))/(t(j-1)+t(j)) + ...
            (1.5*theta(j+1)*t(j+1))/(t(j)+ t(j+1)) );
    end
    
    h(n-1) = t(n-1) * ( 1 + (1.5*theta(n-1)*t(n-2))/(t(n-2)+t(n-1)) );
    
% Now that we have the node spacing values, we want to normalize
% them so that they are within [0,1], with the first data point being
% associated with the value zero and the last data point with the
% value one.

    h = [0;h];
    h = cumsum(h);
    h = h / h(n);
    
% End AFF_ANGLE.M
end