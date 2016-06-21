function GUIDemoHermiteCubic(x , y )
%GUIDEMO A plot with a draggable point. 
%   Proof of concept of idea of using splines to create a plot based of some original input data that a user can then manipulate.  
%   Two points, one line

    close all
    clc
    
    %global variable declaration 
    global Index Grabbed SourceData KnotData;
    Index = []; 
    KnotData = []; %Nx3
    Grabbed = 0; %boolean 
    SourceData = [x,y]; %Nx2
    
    %Set up GUI buttons and figures
    H = figure('Position', [150 100 1100 800],'NumberTitle','off');
    
    % Create buttons and Register event handler functions
    SaveButton = uicontrol(H,'Style','pushbutton','String','Save',...
                'Position',[50 20 60 40],'CallBack', {@SaveClick});
    
    StatsButton = uicontrol(H,'Style','pushbutton','String','Stats',...
                'Position',[50 80 60 40],'CallBack', {@StatsClick});
            
    set(H, 'MenuBar','none');
    set(H,'ToolBar','none');
    set(H, 'Renderer','OpenGL');
    hold on

    %Register event handlers for graph figure
    set(H,'WindowButtonMotionFcn',@MouseMove);
    set(H,'WindowButtonDownFcn',@MouseClick);
    %set(H,'WindowScrollWheelFcn',@MouseScroll);
    %set(H,'KeyPressFcn',@KeyPress );
    
    H.MenuBar = 'figure';
    SetupCurve();
end

function SetupCurve()
    %set starting points
    %set(gcf,'Pointer','FullCrossHair')
    
    global SourceData KnotData;
    KnotData = CalcKnotData(SourceData);    
    UpdateCurve();
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
        if Grabbed == false
            Grabbed = true;
            Index = FindClosest();
        else
           Grabbed = false; 
        end
        %set(gcf,'Pointer','Cross')
%     Right Click
%     elseif strcmpi(get(gcf,'SelectionType'), 'Alt')
%         Grabbed = false;
    end
    UpdateCurve();
end

function SaveClick(~,~)
    global KnotData;
    X = linspace(KnotData(1,1),KnotData(end,1),100);
    Y = evalSpline(X,KnotData(:,1),KnotData(:,2),KnotData(:,3));
    coordinates = [X;Y]';    
    disp('Saving')
    save('SaveData','coordinates','KnotData');
end

function StatsClick(~,~)
    disp(Stats());
end

%Sub functions 
function UpdateCurve()
%Plots points onto graph. 
    global SourceData KnotData;
    cla;
    %Control Points
    plot( KnotData(:,1), KnotData(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 0] );
    
    %Source
    plot( SourceData(:,1) ,SourceData(:,2));
    
    %Spline
    X = SourceData(:,1); %linspace(KnotData(1,1),KnotData(end,1),100);
    Y = evalSpline(X,KnotData(:,1),KnotData(:,2),KnotData(:,3));
    plot( X, Y, 'Color', [255 128 0]/255);
end

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

function [ d ] = Diff( points, xvals, yvals )
%DIFF Uses the central difference method to estimate derivative of a function given discrete sets of inputs and outputs. 
% Uses right and left difference for boundary points.

%points, Nx1 array of input points to find derivative off
%xvals, Mx1 function input values
%yvals, Mx1 function output values
    if( any(size(xvals)~= size(yvals)) || size(xvals,2)~=1 )
        error('dimension function inputs and outputs arrays are incorrect');
    end
    d = zeros(1,length(points));
    for i=1:length(points)
        index = find(xvals == points(i),1);
        if index>1 && index < length(yvals)
            d(i) = ( yvals(index+1)-yvals(index-1) ) / (  xvals(index+1)-xvals(index-1) );
        elseif index == 1
            d(i) = ( yvals(2)-yvals(1) ) / (  xvals(2)-xvals(1) );
        elseif index == length(yvals)
            d(i) = ( yvals(index)-yvals(index-1) ) / (  xvals(index)-xvals(index-1) );
        end
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

%spline stuff functions 
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
    temp = [1;temp;length(Data(:,2))];
    %Create Splines
    kVals = Diff(Data(temp,1),Data(:,1),Data(:,2)); %derivative values
    xVals = Data(temp,1);
    yVals = Data(temp,2);
    
    %xVals = [ xVals ; Data(end,1) + xVals(1) - Data(1,1) ]; %slightly cheating here, really shouldn't use global variables like this
    %yVals = [ yVals ; yVals(1)]; %slightly cheating here, really shouldn't use global variables like this
    %kVals = [ kVals , kVals(1) ]; 
    
    %Create Loop
    
    RetData = [xVals,yVals,kVals'];
end

function [ output ] = evalSpline( x, xvals, yvals, kvals)
%EVALSPLINE Acts as a periodic cubic spline function, takes x values and returns
%function evaluations. The function is defined a set of knots and
%derivative values at those knots. 

%xvals, yvals, and kvals must be 1xn arrays, they are knot information, and
    %functions assumes they are sorted
%is a 1xM array
%x inputs are sorted low to high
%Copied and modified from http://blog.ivank.net/interpolation-with-cubic-splines.html
    %find indices of enclosing knots 
    
%    global SourceData %slightly cheating here, really shouldn't use global variables like this
    indices = zeros(1,length(x));
    output = zeros(1,length(x));
       
    x = mod(x, xvals(end) ); %adjust frame to (0,max) knot
%     xvals = [ xvals ; SourceData(end,1) + xvals(1) - SourceData(1,1) ]; %slightly cheating here, really shouldn't use global variables like this
%     yvals = [ yvals ; yvals(1)]; %slightly cheating here, really shouldn't use global variables like this
%     kvals = [ kvals ; kvals(1) ]; 
    
    for j=1:length(x)      
        modvalue = mod(x(j)-xvals(1), xvals(end)-xvals(1) ); %adjust frame to (0,max) knot
        temp = find(xvals - xvals(1) > modvalue);
        indices(j) = temp(1);     
    end

    %evaluate points
    for j=1:length(x)
        i = indices(j);
%         if i~= 8
            t = ( x(j)- xvals(i-1) ) / ( xvals(i) - xvals(i-1) );

            a =  kvals(i-1) * ( xvals(i) - xvals(i-1) ) - ( yvals(i) - yvals(i-1) );
            b = -kvals(i)*(xvals(i)-xvals(i-1)) + (yvals(i)-yvals(i-1));

            output(j) = (1-t)*yvals(i-1) + t*yvals(i) + t*(1-t)*(a*(1-t)+b*t);
       
%         else     
%             x(j) = x(j) + SourceData(end,1) - xvals(7);
%             t = (x(j)) / ( xvals(i) - xvals(i-1) );
% 
%             a =  kvals(i-1) * ( xvals(i) - xvals(i-1) ) - ( yvals(i) - yvals(i-1) );
%             b = -kvals(i)*(xvals(i)-xvals(i-1)) + (yvals(i)-yvals(i-1));
% 
%             output(j) = (1-t)*yvals(i-1) + t*yvals(i) + t*(1-t)*(a*(1-t)+b*t);
%         end
    end
    
end

%Stat functions 
function statdata = Stats()
%Calculates goodness of fit statistics for current state of curve vs
%original data
%Returns a tab delimited string containing data
    statdata ='no data calculated';
    R2 = RSquared();
    disp('R squared value');
    R2
end

function output = RSquared()
    output = 1 - ( SSresid() ./ SStotal() );
end

function output = SSresid()
    global SourceData KnotData;
    residual = SourceData(:,2) - evalSpline(SourceData(:,1),KnotData(:,1),KnotData(:,2),KnotData(:,3))';
    save residual;
    output = sum(residual.^2);
end

function output = SStotal()
    global SourceData;
    output = (length(SourceData(:,2))-1) * var(SourceData(:,2));
end