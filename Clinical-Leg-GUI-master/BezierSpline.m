%Author: Mark Yeatman
%With code and references from: 
%   cubic Bezier least square fitting by Dr. Murtaza Khan, 27 Jan 2016
%   https://pomax.github.io/bezierinfo/ 
%   http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-elev.html

classdef BezierSpline < handle
    %BEZIER Summary of this class goes here
    %   Detailed explanation goes here
    
    %   Code assumptions:
    %       All splines have the same degree
    
    %   How to use:
    %       Initialize a BezierSpline object to degree 3
    %       Call GenerateCurve with input data and degree 1 tolerance. 
    %       Call RecoverEquations to get list of symbolic equations for
    %       each spline. 
    
    properties
        
        Degree;         %largest x exponent in the equation / order of the curve        
        NumOSplines;    %number of splines that make up curve       
        Weights;        %cell array of x,y pairs of control points ,{i} of (Nx2)
                                
        BinomialLUT = [ 1,0,0,0,0,0,0;           
                        1,1,0,0,0,0,0;         
                        1,2,1,0,0,0,0;         
                        1,3,3,1,0,0,0;       
                        1,4,6,4,1,0,0;       
                        1,5,10,10,5,1,0;    
                        1,6,15,20,15,6,1; ];
        
    end
    
    methods
        
        function obj = BezierSpline(n,W,S)
            % Input:
            %   n = degree of curve
            %   W = list of weights
            %   S = Number of Splines 
            
            % add path for helper functions in fitting algorithm
            path(path,strcat(pwd,'\cubicbezierlsufit'));
            
            % add path for toLatex function
            path(path,strcat(pwd,'\Useful Functions'));
            
            if nargin == 1
                obj.Degree = n;
                obj.Weights = [];
                obj.NumOSplines = [];
            end
            
            if nargin == 2 && size(W,1) == n+1 && size(W,2) == 2
                obj.Degree = n;
                obj.Weights = W;
                obj.NumOSplines = 1;
            end
            
            if nargin == 3
                obj.Degree = n;
                
                if iscell(W)
                    obj.Weights = W;
                else
                    obj.Weights = cell{S};
                    for i=1:S
                       obj.Weights{i} = W(i,:,:);
                    end
                end
                
                obj.NumOSplines = S;
                
            end
        end
        
        function coefficient = Binomial(Bezier,n,k)
        %  Calculates binomial integers, basically the same as n choose k.  
           n = n+1;
           while(n > size(Bezier.BinomialLUT,1))
                Bezier.BinomialLUT = [ Bezier.BinomialLUT , zeros( size(Bezier.BinomialLUT, 1), 1) ];
                s = size(Bezier.BinomialLUT,2);
                prev = s-1;
                nextRow = zeros( 1, size(Bezier.BinomialLUT,1) );
                nextRow(1) = 1;
                for i=2:prev
                  nextRow(i) = Bezier.BinomialLUT(prev,i-1) + Bezier.BinomialLUT(prev,i);
                end
                nextRow(s) = 1;
                Bezier.BinomialLUT = [Bezier.BinomialLUT ; nextRow] ; 
           end
           
           coefficient = Bezier.BinomialLUT(n,k);
           
        end
        
        function [x,y] = Crunch(Bezier,t,j)
        % Calculates x,y coordinate pairs based off parametric t value and current weights/degree 
        % Inputs: 
        %   t = parametric variable
        %   j = spline #
        
        % Output:
        %   2 Nx1 matrices of x,y coordinates
            W = Bezier.Weights{j};
            n = Bezier.Degree ;
            x=0;
            y=0;
            for k = 1:n+1
                xterm = W(k,1) .* Bezier.Binomial(n,k) .* (1-t).^(n-k+1) .* t.^(k-1) ;
                yterm = W(k,2) .* Bezier.Binomial(n,k) .* (1-t).^(n-k+1) .* t.^(k-1) ;

                x = x + xterm;
                y = y + yterm;
            end                          
        end
        
        function [x,y,BreakIndices] = GenerateCurve(Bezier,DataX,DataY,MxAllowSqD,varargin)
            
            if nargin == 4              
                BreakIndices = Bezier.Fit(DataX,DataY,MxAllowSqD);

                x = [];
                y = [];

                for j=1:Bezier.NumOSplines;

                    %takes care of point overlap between splines, insures same
                    %number of data points are output
                    if j ~= Bezier.NumOSplines
                        t = linspace(0,1,BreakIndices(j+1)-BreakIndices(j));
                    else
                        t = linspace(0,1,BreakIndices(j+1)-BreakIndices(j)+1); 
                    end

                    [tempX,tempY] = Bezier.Crunch(t, j);
                    x = [x,tempX];
                    y = [y,tempY];
                end
                x = x';
                y = y';
            end
            if nargin == 2
                BreakIndices = DataX;

                x = [];
                y = [];

                for j=1:Bezier.NumOSplines;

                    %takes care of point overlap between splines, insures same
                    %number of data points are output
                    if j ~= Bezier.NumOSplines
                        t = linspace(0,1,BreakIndices(j+1)-BreakIndices(j));
                    else
                        t = linspace(0,1,BreakIndices(j+1)-BreakIndices(j)+1); 
                    end

                    [tempX,tempY] = Bezier.Crunch(t, j);
                    x = [x,tempX];
                    y = [y,tempY];
                end
                x = x';
                y = y';
            end
        end
            
        function [ EqList ] = RecoverEquations(Bezier,S,varargin)
            % Returns Nx2 Cell array of symbolic parametric equations for
            % each spline. 
            % First Column is x(t), 2nd is y(t)
            % Enter 'show' for S arguement to get latex display of
            % equations. 
            EqList = cell(Bezier.NumOSplines,2);  
            
            for i = 1:Bezier.NumOSplines
               [EqList{i,1},EqList{i,2}] = Bezier.RecoverEquation(i);
               if nargin ==2 
                   if strcmp(S,'show')
                       toLatex(EqList{i,1});
                       toLatex(EqList{i,2});
                   end
               end
            end
        end
        
        function [X,Y] = RecoverEquation(Bezier,SplineNum)
        % Creates a list of paired parametric symbolic equations for each
        % spline based on weight values. 
            X = 0;
            Y = 0;
            syms t;
            n = Bezier.Degree ;
            W = Bezier.Weights{SplineNum};
            
            for k = 1:n+1
                xterm = W(k,1) * Bezier.Binomial(n,k) * (1-t)^(n-k+1) * t^(k-1) ;
                yterm = W(k,2) * Bezier.Binomial(n,k) * (1-t)^(n-k+1) * t^(k-1) ;
                
                X = X + xterm;
                Y = Y + yterm;
            end
            
        end
        
        function fbi = Fit(Bezier,DataX, DataY, MxAllowSqD)
        % Algorithms from Dr. Murtaza Khan  
        % Only works for cubic splines
        % Returns indices of spline break points. 
        % Requires many many functions from the 'cubicbezierlsufit' folder
            if Bezier.Degree ~= 3
               error('Curve must be degree 3'); 
            end
            
            if size(DataX,1) ==1
                DataX = DataX';
            end
            
            if size(DataY,1) ==1
                DataY = DataY';
            end
            
            Mat=[DataX,DataY];
            ei= length(DataX);
            ibi=[1;ei]; %first and last point are taken as initial break points
            
            [p0mat,p1mat,p2mat,p3mat,fbi]=bzapproxu(Mat,MxAllowSqD,ibi);    %fitting function, external to class
            
            Bezier.NumOSplines = size(p0mat,1);
            Bezier.Weights = cell(Bezier.NumOSplines,1);
            
            %convert Khan's output format in objet data format
            for i=1:Bezier.NumOSplines
               Bezier.Weights{i} = [p0mat(i,:);p1mat(i,:);p2mat(i,:);p3mat(i,:)]; 
            end           
        end
        
        function RaiseTo(Bezier,n)
        %Algorithm taken directly from http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-elev.html  
        
            %must run (n-original degree curve) times
            for p = Bezier.Degree+1:n
                %for all spline functions
                for j=1:Bezier.NumOSplines
                    
                    
                    W = [Bezier.Weights{j}]; 
                    
                    NewW = [Bezier.Weights{j}(j,:);
                            zeros(p-1,2);
                            Bezier.Weights{j}(end,:)]; %add another row for new control point
                        
                    %execute n -> n+1 raise procedure
                    for i=2:p
                        NewW(i,:) = W(i-1,:).*(i-1)/(p) + ( 1 - (i-1)/p).*W(i,:) ;
                    end
                    
                    Bezier.Weights{j} = NewW;
                    
                end
            end
            Bezier.Degree = n;
        end
    end
    
end

