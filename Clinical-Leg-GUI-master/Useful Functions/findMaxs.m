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

