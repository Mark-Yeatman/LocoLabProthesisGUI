function GUIDemoHermiteCubic(x , y )
%GUIDEMO A plot with a draggable point. 
%   Proof of concept of idea of using splines to create a plot based of some original input data that a user can then manipulate.  
%   Two points, one line

    close all
    clc
    
    global Index Grabbed SourceData KnotData;
    Index = [];
    KnotData = [];
    Grabbed = 0;
    SourceData = [x,y];
    
    H = figure('Position', [150 100 1100 800],'NumberTitle','off');
    SaveButton = uicontrol(H,'Style','pushbutton','String','Save',...
                'Position',[50 20 60 40],'CallBack', {@SaveClick});
            
    set(H, 'MenuBar','none');
    set(H,'ToolBar','none');
    set(H, 'Renderer','OpenGL');
    hold on

    set(H,'WindowButtonMotionFcn',@MouseMove);
    set(H,'WindowButtonDownFcn',@MouseClick);
    %set(H,'WindowScrollWheelFcn',@MouseScroll);
    %set(H,'KeyPressFcn',@KeyPress );
    
    H.MenuBar = 'figure';
    SetupCurve();
end

%Event functions 
function MouseMove(~,~)
    global Grabbed Index KnotData; %Grabbed MoveAlong Snap
    if Grabbed
        CurrentPoints = get(gca,'CurrentPoint');
        KnotData(Index,1:2) = CurrentPoints(1,1:2);
        UpdateCurve();
    end
end

function MouseClick(~,~)
    global Grabbed Index; %SelectMode Grabbed Snap Points Weights Xi Degree
    %Left Click
    if strcmpi(get(gcf,'SelectionType'), 'Normal')         
        Grabbed = true;
        Index = FindClosest();
        set(gcf,'Pointer','Cross')
    %Right Click
    elseif strcmpi(get(gcf,'SelectionType'), 'Alt')
        Grabbed = false;
    end
    UpdateCurve();
end

function SaveClick(~,~)
    global KnotData;
    X = linspace(KnotData(1,1),KnotData(end,1),100);
    Y = evalSpline(X,KnotData(:,1),KnotData(:,2),KnotData(:,3));
    temp = [X;Y]';
    save('SaveData','temp');
end

%Sub functions 
function [maxima] = findMaxs(x)
%this function finds the indices of the maximums of a periodic function, where x is
%a data set representing one period of the function. 
    x = x(:);
    % Identify whether signal is rising or falling
    upordown = sign(diff(x));
    % Find points where signal is rising before, falling after
    maxflags = [upordown(1)<0; diff(upordown)<0; upordown(end)>0];
    maxima_1   = find(maxflags);
    maxima=maxima_1;
    if length(maxima_1) > 1
        shift = maxima_1(2);
        x = circshift(x,shift);
        % Identify whether signal is rising or falling
        upordown = sign(diff(x));
        % Find points where signal is rising before, falling after
        maxflags = [upordown(1)<0; diff(upordown)<0; upordown(end)>0];
        maxima_2   = find(maxflags);
        maxima_2 = maxima_2 -shift;
        for i=1:length(maxima_2)
            if maxima_2(i)<=0
                maxima_2(i)=length(x)+maxima_2(i);
            end
        end   
        maxima = intersect(maxima_1,maxima_2);
    end
end

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

function Index = FindClosest()
%Calculates the point in the plot data set closest to where the user
%clicked. 
    global KnotData;
    
    CurrentPoint = get(gca,'CurrentPoint');
    Pt = CurrentPoint(1,1:2);

    Distances = zeros(1,length(KnotData));

    for P = 1:length(Distances)
        Distances(P) = norm( KnotData(P,1:2)-Pt, 2);
    end
    [~,Index] = min( Distances );
end

function RetData = CalcKnotData(Data)
%This function descritizes input data into multiple splines, determines
%starting and ending points of those splines, and estimates values of
%the derivatives of the data at those points. 
    KneeMaxs = findMaxs(Data(:,2));
    KneeMins = findMaxs(-Data(:,2));
    temp = sort([KneeMins;KneeMaxs]);
    extrapoints = zeros(length(temp)-1,1);
    %extra points are the half way points between local mins and maxs, so
    %that a local min or max ends up as roughly the mid point of a spline.
    %Ends up creating a spline for every local min or max.
    for i=1:length(extrapoints)
       extrapoints(i) = floor( ( temp(i) + temp(i+1) ) / 2 );
    end

    temp = sort( [ temp ; extrapoints ] );
    
    %Create Splines
    kVals = centDiff(Data(temp,1),Data(:,1),Data(:,2)); %derivative values
    xVals = Data(temp,1);
    yVals = Data(temp,2);
    RetData = [xVals,yVals,kVals'];
end

function [ q ] = evalSpline( x, xvals, yvals, kvals)
%EVALSPLINE Acts as a cubic spline function, takes x values and returns
%function evaluations. The function is defined a set of knots and
%derivative values at those knots. 

%xvals, yvals, and kvals must be 1xn arrays
%is a 1xM array
%x inputs are sorted low to high
%Copied and modified from http://blog.ivank.net/interpolation-with-cubic-splines.html
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

function SetupCurve()
    %set starting points
    set(gcf,'Pointer','FullCrossHair')
    
    global SourceData KnotData;
    KnotData = CalcKnotData(SourceData);    
    UpdateCurve();
end

function UpdateCurve()
%Plots points onto graph. 
    global SourceData KnotData;
    cla;
    %Control Points
    plot( KnotData(:,1), KnotData(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0] );
    
    %Source
    plot( SourceData(:,1) ,SourceData(:,2));
    
    %Spline
    X = linspace(KnotData(1,1),KnotData(end,1),100);
    Y = evalSpline(X,KnotData(:,1),KnotData(:,2),KnotData(:,3));
    plot( X, Y, 'Color', [255 128 0]/255, 'LineWidth', 3 );
end