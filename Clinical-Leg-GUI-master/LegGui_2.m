% Author: Mark Yeatman and Frank Moreno
% Written: Summer, 2016

%TO DO: Fix the dragging of controls forming spline start/end.
function LegGui_2
%GUIDEMO A plot with a draggable point. 
%   Proof of concept of idea of using splines to create a plot based of some original input data that a user can then manipulate.  
%   Two points, one line

    close all
    clc
    load Winter_NWvsPV.mat
    x = Winter_NWvsPV(:,1);
    y = Winter_NWvsPV(:,2);
    %global variable declaration 
    global Index Grabbed SourceData BezierData BezierCurve StartPoint EndPoint...
        n h1 h2 h3 h4 h5;

    BezierCurve = BezierSpline(3);
    Index = [];  %indexs of flattened weight array in BezierSpline
    Grabbed = 0; %boolean 
    SourceData = [x,y];
    BezierData = [];%zeros(length(x),2);
    StartPoint = [];
    EndPoint = [];
    n = 1;                          %Number after the saved lab files
    h1 = '  ';
    h2 = '  ';
    h3 = '  ';
    h4 = '  ';
    h5 = '  ';
    %Set up GUI buttons and figures
    H = figure('Visible','on','Position',[150 50 1100 775]);     %setting size of figure
    hp = uipanel(H,'Title','Controls','FontSize',12,'BackgroundColor',[0.70 0.7 0.7],...
                'Units','pixels','Position',[715 70 380 650]);
    graph_axes = axes('Units','pixels','Position',[50,70,650,650]);
    % Create buttons and Register event handler functions
    SaveButton = uicontrol('Style','pushbutton','String','Save',...
                'Position',[750 600 130 70],...
                'Callback',{@SaveClick});  
    moreinfo =    uicontrol('Style','popupmenu',...
                'String',{'More Info...','Heel Strike','Foot Flat','Mid-Stance','Heel-Off','Toe-Off','Mid-Swing'},...
                'Position',[930 400 150 70],...
                'Callback',{@moreinfo_Callback});
    graphchoice = uicontrol('Style','listbox','String',{h1,h2,h3,h4,h5},...
                'Max',2,'Min',0,'Value',[],...
                'Position',[930 600 130 70],...
                'Callback',{@graphchoice_Callback});
    choicetitle = uicontrol('Style','text','String',('Saved Files'),...
                'Position',[960 675 70 20]);
    StatsButton = uicontrol(H,'Style','pushbutton','String','Stats',...
                'Position',[750 400 130 70],'CallBack', {@StatsClick});
    ManInput = uicontrol('Style','pushbutton','String','Manual Input',...
                'Position',[750 250 130 70],...
                'Callback',{@ManInput}); 
    title('Gait Cycle');  
    xlabel('Percent of Gait Cycle');
    ylabel('Angle of Knee');
    set(H,'Name','Gait Cycle'); 
    set(H, 'Color',[0.70,0.7,0.7]);
    set(H, 'MenuBar','none');
    set(H,'ToolBar','none');
    set(H, 'Renderer','OpenGL');
    hold on

    %Register event handlers for graph figure
    set(H,'WindowButtonMotionFcn',@MouseMove);
    set(H,'WindowButtonDownFcn',@MouseClick);
    %set(H,'WindowScrollWheelFcn',@MouseScroll);
    set(H,'KeyPressFcn',@KeyPress );
    
    H.MenuBar = 'figure';
    SetupCurve();
    graphBoundaries();
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
            [x,y,I]=BezierCurve.GenerateCurve(SourceData(StartPoint:EndPoint,1),SourceData(StartPoint:EndPoint,2),MxAllowSqD);
            BezierData = [x,y]; 
        %Select Control Point
        else                      
            if Grabbed == false
                Grabbed = true;
                Index = FindClosest();
            else
               Grabbed = false; 
            end
        end
    end    
    UpdateCurve();
end

function SaveClick(~,~)
    global BezierData SourceData StartPoint EndPoint n h1 h2 h3 h4 h5;
    X = [SourceData(1:StartPoint,1);BezierData(:,1);SourceData(EndPoint:end,1)];
    Y = [SourceData(1:StartPoint,end);BezierData(:,end);SourceData(EndPoint:end,end)];
    coordinates = [X,Y];
    file_input = inputdlg('Enter the name of the file.',...
                     'Input',1);
    file_name = strcat(file_input{1},'.mat');
    save(file_name,'coordinates');  
    if(n == 1)
            h1 = file_name;
        elseif (n == 2)
            h2 = file_name;
        elseif (n == 3)
            h3 = file_name;
        elseif (n == 4)
            h4 = file_name;
        elseif (n == 5)
            h5 = file_name;
     end
        graphchoice = uicontrol('Style','listbox','String',{h1,h2,h3,h4,h5},...
                'Max',2,'Min',0,'Value',[],...
                'Position',[930 600 130 70],...
                'Callback',{@graphchoice_Callback});
         n = n + 1;
         RemovingGraph();
end

function StatsClick(~,~)
   Stats();
end

function KeyPress(~,~)
if strcmpi(get(gcf,'CurrentCharacter'), 'x')
    RemovingGraph();
end    
end

function moreinfo_Callback(source, ~)
    global msgdsp shape;
        msgdsp = uicontrol('Style','text','String','','Position',[760 120 300 70]);
        str = source.String;
        val = source.Value;
        switch str{val};
            case 'More Info...'
                ClearObjects();
            case 'Heel Strike'
                message = sprintf('%s\n%s\n%s', ...
                '- First 0% to 2% of gait cycle',...
                '- Ideally first contact for foot; begins when heel makes contact',...
                '- Body weight is immediately transferred to leg');
%                 borderx = [0.00 0.00 0.02 0.02];
%                 bordery = [-70 10 10 -70 ];
%                 shape = fill(borderx,bordery,'blue');
%                 set(shape,'FaceAlpha',0.1); % or alpha(h,0.1); would also work
%                 set(shape,'EdgeColor','None');
                msgdsp = uicontrol('Style','text','String',message,'Position',[760 120 300 70]);
            case 'Foot Flat'
               message = sprintf('%s\n%s\n%s', ...
                '- Approximately from 2% to 12% of the gait cycles',...
                '- Following heel strike, the rest of the foot makes contact with floor',...
                '- Absorbs shock and provides stability for other leg');
%                 borderx = [0.02 0.02 0.12 0.12];
%                 bordery = [-70 10 10 -70 ];
%                 shape = fill(borderx,bordery,'blue');
%                 set(shape,'FaceAlpha',0.1); % or alpha(h,0.1); would also work
%                 set(shape,'EdgeColor','None');
                msgdsp = uicontrol('Style','text','String',message,'Position',[760 120 300 70]);
            case 'Mid-Stance'
                message = sprintf('%s\n%s\n%s', ...
                '- Approximately from 12% to 31% of the gait cycles',...
                '- Continues from opposite leg being lifted and making contact again',...
                '- During this stage, reference foot must support all of the weight');
%                 borderx = [0.12 0.12 0.31 0.31];
%                 bordery = [-70 10 10 -70 ];
%                 shape = fill(borderx,bordery,'blue');
%                 set(shape,'FaceAlpha',0.1); % or alpha(h,0.1); would also work
%                 set(shape,'EdgeColor','None');
                msgdsp = uicontrol('Style','text','String',message,'Position',[760 120 300 70]);
            case 'Heel-Off'
                message = sprintf('%s\n%s\n%s', ...
                '- Approximately from 31% to 50% of the gait cycles',...
                '- Begins as reference heel begins to rise until opposite foot makes contact',...
                '- Load of body weight is trasnferred during this stage');
%                 borderx = [0.31 0.31 0.50 0.50];
%                 bordery = [-70 10 10 -70 ];
%                 shape = fill(borderx,bordery,'blue');
%                 set(shape,'FaceAlpha',0.1); % or alpha(h,0.1); would also work
%                 set(shape,'EdgeColor','None');
                msgdsp = uicontrol('Style','text','String',message,'Position',[760 120 300 70]);
            case 'Toe-Off'
                message = sprintf('%s\n%s\n%s', ...
                '- Approximately from 50% to 75% of the gait cycles',...
                '- Weight is then transferred to reference foot as limb pushes off ground',...
                '- Foot is lifted from the floor at 62% and swings until it is adjacent to other leg');
%                 borderx = [0.50 0.50 0.75 0.75];
%                 bordery = [-70 10 10 -70 ];
%                 shape = fill(borderx,bordery,'blue');
%                 set(shape,'FaceAlpha',0.1); % or alpha(h,0.1); would also work
%                 set(shape,'EdgeColor','None');
                msgdsp = uicontrol('Style','text','String',message,'Position',[760 120 300 70]);
            case 'Mid-Swing'
                message = sprintf('%s\n%s\n%s', ...
                '- Approximately from 75% to 87% of the gait cycles',...
                '- Middle third of the overall swing stage',...
                '- Continues from leg being opposite to adjacant foot until it is straight');
%                 borderx = [0.75 0.75 0.87 0.87];
%                 bordery = [-70 10 10 -70 ];
%                 shape = fill(borderx,bordery,'blue');
%                 set(shape,'FaceAlpha',0.1); % or alpha(h,0.1); would also work
%                 set(shape,'EdgeColor','None');
                msgdsp = uicontrol('Style','text','String',message,'Position',[760 120 300 70]);
        end
end

function graphchoice_Callback(source, ~)
        items = get(source, 'String');
        index_selected = get(source,'Value');
        item_selected = items{index_selected};
        final = load(item_selected);
        x = final.coordinates(:,1);
        y = final.coordinates(:,2); 
        plot(x,y);
        legend(item_selected);
    end
    
function ManInput(~,~)
global StartPoint EndPoint BezierCurve BezierData SourceData;
prompt = {'Enter the starting point:','Enter the end point'};
default = {num2str(StartPoint),num2str(EndPoint)};
RemovingGraph();
answer = inputdlg(prompt,'input',1,default);
StartPoint = str2num(answer{1});
EndPoint =  str2num(answer{2});
MxAllowSqD=0.1; % Max. allowed Square Distance between original and fitted data    
[x,y,I]=BezierCurve.GenerateCurve(SourceData(StartPoint:EndPoint,1),SourceData(StartPoint:EndPoint,2),MxAllowSqD);
BezierData = [x,y]; 
UpdateCurve();
end

%Event helper functions
function UpdateCurve()
%Plots points onto graph. 
    global SourceData BezierData BezierCurve I StartPoint EndPoint Red1 Red2 GrnLine BlckX ClrX;
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
        for j= 1:length(W)
            if mod(j,4)==1
                BlckX = plot( W(j,1),W(j,2), '*', 'color' , 'black' );
            else
                ClrX = plot( W(j,1),W(j,2), '*', 'color' , colors(i,:) );
            end           
        end         
    end
    [x,y,I]=BezierCurve.GenerateCurve(I);
    BezierData = [x,y];
    %Spline and Source Data
    
    if ~isempty(StartPoint)
        Red1 = plot( SourceData(StartPoint,1),SourceData(StartPoint,2),'*','color','red');
    end
    
    if ~isempty(EndPoint)
        Red2 = plot( SourceData(EndPoint,1),SourceData(EndPoint,2),'*','color','red');
    end    
    
    if ~isempty(BezierData)
       GrnLine = plot( BezierData(:,1),BezierData(:,2),'linewidth',2,'color','green'); 
    end
    
    plot( SourceData(:,1),SourceData(:,2),'color','blue');
    graphBoundaries();
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

function Index = FindClickedPoint()
    global SourceData;

    CurrentPoint = get(gca,'CurrentPoint');
    Pt = CurrentPoint(1,1:2);

    Distances = zeros(1,length(SourceData));

    for P = 1:length(Distances)
        Distances(P) = norm( SourceData(P,1:2)-Pt, 2);
    end
    [~,Index] = min( Distances );
end

function RemovingGraph()
    global StartPoint EndPoint Red1 Red2 GrnLine BlckX ClrX BezierCurve;
    StartPoint = [];
    EndPoint = [];
    delete(Red1);
    delete(Red2);
    delete(GrnLine);
    delete(BlckX);
    delete(ClrX);
    delete(BezierCurve);
    BezierCurve = BezierSpline(3);
end

%Stat functions 
function statdata = Stats()
%Calculates goodness of fit statistics for current state of curve vs
%original data
%Returns a tab delimited string containing data
    statdata ='no data calculated';
    R2 = RSquared();
    %disp('R squared value');
    R = num2str(R2);
    message = sprintf('%s\n%s%s',statdata,'R squared value: ',R);
    msgbox(message);
    
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

function graphBoundaries()
    %w = [0.62,0.62];
    %z = [-70,10];
    borderx = [0.62 0.62 1 1];
    bordery = [-70 10 10 -70 ];
    shape = fill(borderx,bordery,[0.4 0.4 0.4]);
    set(shape,'FaceAlpha',0.1); % or alpha(h,0.1); would also work
    set(shape,'EdgeColor','None');
    %plot(w,z,'--k');
end

function ClearObjects()
    clear msgdsp;
    clear shape
end