clear all ;
load Winter_NWvsPV;
B = BezierSpline(3);
[x,y]=B.GenerateCurve(Winter_NWvsPV(:,1),Winter_NWvsPV(:,2),1);

colors = [ 0,0,0
           0.5,0.5,1;
           1,0,1;
           0,1,1;
           1,0,0;
           0,1,0;
           0,0,1;
           1,0.5,0.5];
       
plot(x,y);
hold on;

for i=1:B.NumOSplines
    W = B.Weights{i};
    plot( W(:,1),W(:,2), '*', 'color' , colors(i,:) );   
end

hold off;