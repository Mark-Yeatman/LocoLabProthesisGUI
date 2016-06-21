% Author: Mark Yeatman
% Written: Summer, 2016

function GUIDemoBezierCubic(x , y )
%GUIDEMO A plot with a draggable point. 
%   Proof of concept of idea of using splines to create a plot based of some original input data that a user can then manipulate.  
%   Two points, one line

    close all
    clc
    
    %global variable declaration 
    global Index Grabbed SourceData BezierData BezierCurve StartPoint EndPoint;

    BezierCurve = BezierSpline(3);
    Index = [];  %indexs of flattened weight array in BezierSpline
    Grabbed = 0; %boolean 
    SourceData = [x,y]; %Nx2
    BezierData = [];%zeros(length(x),2);
    StartPoint = [];
    EndPoint = [];
    
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
    
    %global SourceData BezierCurve BezierData I;   
    
    %Bezier Setup
    %MxAllowSqD=0.1; % Max. allowed Square Distance between original and fitted data    
    %[x,y,I]=BezierCurve.GenerateCurve(SourceData(:,1),SourceData(:,2),MxAllowSqD);
    %BezierData = [x,y]; 
    UpdateCurve();
    
end


%Event functions 
function MouseMove(~,~)
    global Grabbed Index BezierCurve;
    if Grabbed
        KnotData = cell2mat(BezierCurve.Weights);
        CurrentPoints = get(gca,'CurrentPoint');
%         if mod(Index,4) == 0 %Potential Double Knot Point
%             if Index == 1 || Index == length(KnotData) % First or last point in entire set
%                 KnotData(Index,1:2) = CurrentPoints(1,1:2);
%             else    %is a double knot point
%                 KnotData(Index,1:2) = CurrentPoints(1,1:2);
%                 KnotData(Index+1,1:2) = CurrentPoints(1,1:2);
%             end
%         else
%             KnotData(Index,1:2) = CurrentPoints(1,1:2);
%         end
        KnotData(Index,1:2) = CurrentPoints(1,1:2);  
        BezierCurve.Weights = mat2cell(KnotData,ones(length(BezierCurve.Weights),1)*length(BezierCurve.Weights{1}),2);
        UpdateCurve();
    end
end

function MouseClick(~,~)
    
    global Grabbed Index StartPoint EndPoint BezierData BezierCurve I SourceData; 
    
    %Left Click
    if strcmpi(get(gcf,'SelectionType'), 'Normal')     
        
        %Select Starting Point
        if isempty(StartPoint)
            StartPoint = FindClickedPoint();
            
        %Select Ending Point
        elseif isempty(EndPoint)
            EndPoint = FindClickedPoint();
            MxAllowSqD=0.1; % Max. allowed Square Distance between original and fitted data    
            [~,~,I]=BezierCurve.GenerateCurve(SourceData(StartPoint:EndPoint,1),SourceData(StartPoint:EndPoint,2),MxAllowSqD);
            BezierCurve.RaiseTo(5);
            [x,y,I]=BezierCurve.GenerateCurve(I);
            save BezierCurve
            BezierData = [x,y]; 
        %Select Control Point
        else                      
            if Grabbed == false
                
                Index = FindClosest();
                if Index == 0
                    Grabbed = false;
                else
                    Grabbed = true;
                end
            else
               Grabbed = false; 
            end
        end
    end
    UpdateCurve();
end

function SaveClick(~,~)
    global BezierData BezierCurve;
    save('SaveData','BezierData','BezierCurve');
end

function StatsClick(~,~)
    disp(Stats());
end


%Event helper functions
function UpdateCurve()
%Plots points onto graph. 
    global SourceData BezierData BezierCurve I StartPoint EndPoint;
    cla;
    
    %Array of normalized RGB values. 
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
        %for j= 1:length(W)
            %if mod(j,length(W))==1
             %   plot( W(j,1),W(j,2), '*', 'color' , 'black' );
            %else
                plot( W(3:4,1),W(3:4,2), '*', 'color' , colors(i,:) );
            %end           
        %end         
    end
    [x,y,I]=BezierCurve.GenerateCurve(I);
    BezierData = [x,y];
    %Spline and Source Data
    
    if ~isempty(StartPoint)
        plot( SourceData(StartPoint,1),SourceData(StartPoint,2),'*','color','red');
    end
    
    if ~isempty(EndPoint)
        plot( SourceData(EndPoint,1),SourceData(EndPoint,2),'*','color','red');
    end    
    
    if ~isempty(BezierData)
       plot( BezierData(:,1),BezierData(:,2),'linewidth',2,'color','green'); 
    end
    
    if ~isempty(SourceData) && ~isempty(EndPoint)
        plot( SourceData(1:StartPoint,1), SourceData(1:StartPoint,2),'color','blue');
        plot( SourceData(StartPoint+1:EndPoint-1,1),SourceData(StartPoint+1:EndPoint-1,2),'--','color','blue');
        plot( SourceData(EndPoint:end,1), SourceData(EndPoint:end,2),'color','blue');
    else
        plot( SourceData(:,1),SourceData(:,2),'color','blue');
    end
end

function Index = FindClosest()
%Calculates the ctrlpoint closest to where the user clicked. 
    global BezierCurve SourceData;
    msf = 0.005 ; %magic scale factor
    maxAllowDistX =   msf*( max(SourceData(:,1)) - min(SourceData(:,1)) );
    maxAllowDistY =   msf*( max(SourceData(:,2)) - min(SourceData(:,2)) );
       
    KnotData = cell2mat(BezierCurve.Weights);
    CurrentPoint = get(gca,'CurrentPoint');
    Pt = CurrentPoint(1,1:2);

    Distances = zeros(1,length(KnotData));

    for P = 1:length(Distances)
        diffVector = KnotData(P,1:2)-Pt;
        %Scale vector so that calculation matches up with visual scale. 
        diffVector(1) = diffVector(1)/( max(SourceData(:,1)) - min(SourceData(:,1))  );
        diffVector(2) = diffVector(2)/( max(SourceData(:,2)) - min(SourceData(:,2)) );
        
        Distances(P) = norm(diffVector,2);
    end
    [~,Index] = min( Distances );
    
    %see selection area gaphically
    %th = 0:pi/50:2*pi;
    %xunit = maxAllowDistX * cos(th) + KnotData(3,1);
    %yunit = maxAllowDistY * sin(th) + KnotData(3,2);
    
    %plot(xunit, yunit,Pt(1,1),Pt(1,2),'+');
    
    if (Index ~= 3 && Index ~= 4) || abs(KnotData(Index,1)-Pt(1,1)) > maxAllowDistX || abs(KnotData(Index,2)-Pt(1,2)) > maxAllowDistY  %cheat and only report on the middle ctrl points to hide end/intentional hidden points
        Index = 0;
    end
    
end

function Index = FindClickedPoint()
%Calculates the point in the plot data set closest to where the user
%clicked. 
    global SourceData;
    
    CurrentPoint = get(gca,'CurrentPoint');
    Pt = CurrentPoint(1,1:2);

    Distances = zeros(1,length(SourceData));

    for P = 1:length(Distances)
        Distances(P) = norm( SourceData(P,1:2)-Pt, 2);
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