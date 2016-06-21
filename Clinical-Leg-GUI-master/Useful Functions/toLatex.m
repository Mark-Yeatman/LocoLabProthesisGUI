function [] = toLatex(fx)

    persistent count;
    persistent h;
    
    
    if ismatrix(fx)
       close all  
       
       h = figure; 
       
       set(h, 'color', 'white'); axis off;
       count =0;
       
       for i=1:size(fx,1)
           for j=1:size(fx,2)
               count = count + 1;
               string = latex(fx(i,j)); 
               my_text = strcat( '$$' ,string ,'$$');
               text('units','pixels', 'position', [2,50*count], 'fontsize', 10, 'color', 'k', ...
               'string', strcat( 'Eq #: ', strcat(num2str(i) ,',',num2str(j) ) ) );
               text('units','pixels', 'position', [60,50*count], 'fontsize', 14, 'color', 'k', ...
               'Interpreter', 'latex', 'string', my_text);
           end
       end
       return
    end
    
    string = latex(fx); 
    my_text = strcat( '$$' ,string ,'$$');
    
    if isempty(count)
        count =0;
    end
    
    if isempty(h)
       h=99; %takes care of first function call
    end
    
    if ~ishandle(h)
        count = 0;
        h = figure; %takes care of function call after closing window
    end
    
    if ishandle(h)
        set(h, 'color', 'white'); axis off;       
    end
            
    count = count + 1;
    text('units','pixels', 'position', [2,50*count], 'fontsize', 10, 'color', 'k', ...
    'string', strcat( 'Eq #: ', num2str(count)) );
    text('units','pixels', 'position', [50,50*count], 'fontsize', 14, 'color', 'k', ...
    'Interpreter', 'latex', 'string', my_text);

end

