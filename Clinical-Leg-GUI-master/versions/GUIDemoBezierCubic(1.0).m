% Author: Mark Yeatman
% Written: Summer, 2016

%TO DO: Fix the dragging of controls forming spline start/end.
function GUIDemoBezierCubic(x , y )
%GUIDEMO A plot with a draggable point. 
%   Proof of concept of idea of using splines to create a plot based of some original input data that a user can then manipulate.  
%   Two points, one line

    close all
    clc
    
    %global variable declaration 
    global Index Grabbed SourceData BezierData BezierCurve;

    BezierCurve = BezierSpline(3);
    Index = [];  %indexs of flattened weight array in BezierSpline
    Grabbed = 0; %boolean 
    SourceData = [x,y]; %Nx2
    BezierData = zeros(length(x),2);
    
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
    
    global SourceData BezierCurve BezierData I;   
    
    %Bezier Setup
    MxAllowSqD=0.1; % Max. allowed Square Distance between original and fitted data    
    [x,y,I]=BezierCurve.GenerateCurve(SourceData(:,1),SourceData(:,2),MxAllowSqD);
    BezierData = [x,y]; 
    UpdateCurve();
    
end


%Event functions 
function MouseMove(~,~)
    global Grabbed Index BezierCurve;
    if Grabbed
        KnotData = cell2mat(BezierCurve.Weights);
        CurrentPoints = get(gca,'CurrentPoint');
        if mod(Index,4) == 0 %Potential Double Knot Point
            if Index == 1 || Index == length(KnotData) % First or last point in entire set
                KnotData(Index,1:2) = CurrentPoints(1,1:2);
            else    %is a double knot point
                KnotData(Index,1:2) = CurrentPoints(1,1:2);
                KnotData(Index+1,1:2) = CurrentPoints(1,1:2);
            end
        else
            KnotData(Index,1:2) = CurrentPoints(1,1:2);
        end
        
        BezierCurve.Weights = mat2cell(KnotData,ones(length(BezierCurve.Weights),1)*4,2);
        UpdateCurve();
    end
end

function MouseClick(~,~)
    global Grabbed Index; 
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
    global BezierData;
    save('SaveData','BezierData');
end

function StatsClick(~,~)
    disp(Stats());
end


%Event helper functions
function UpdateCurve()
%Plots points onto graph. 
    global SourceData BezierData BezierCurve I;
    cla;
    
    colors = [ 0.2,1,1
               0.5,0.5,1;
               1,0,1;
               0,1,1;
               1,0,0;
               0,1,0;
               0,0,1;
               1,0.5,0.5;
               0.3,0.8,0.4;
               0.8,0.3,0.7];
    
    %Control Points
    
    for i=1:BezierCurve.NumOSplines
        W = BezierCurve.Weights{i};
        for j= 1:length(W)
            if mod(j,4)==1
                plot( W(j,1),W(j,2), '*', 'color' , 'black' );
            else
                plot( W(j,1),W(j,2), '*', 'color' , colors(i,:) );
            end           
        end         
    end
    [x,y,I]=BezierCurve.GenerateCurve(I);
    BezierData = [x,y];
    %Spline and Source Data
    plot( SourceData(:,1) ,SourceData(:,2), BezierData(:,1) ,BezierData(:,2));
    
end

function Index = FindClosest()
%Calculates the point in the plot data set closest to where the user
%clicked. 
    global BezierCurve;
    
    KnotData = cell2mat(BezierCurve.Weights);
    CurrentPoint = get(gca,'CurrentPoint');
    Pt = CurrentPoint(1,1:2);

    Distances = zeros(1,length(KnotData));

    for P = 1:length(Distances)
        Distances(P) = norm( KnotData(P,1:2)-Pt, 2);
    end
    [~,Index] = min( Distances );
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
    global SourceData BezierData;
    residual = zeros( length(SourceData(:,1)) ,2);
    residual(:,2) = SourceData(:,2) - BezierData(:,2);
    save residual;
    output = sum(residual.^2);
end

function output = SStotal()
    global SourceData;
    output(1) = (length(SourceData(:,2))-1) * var(SourceData(:,2));
end