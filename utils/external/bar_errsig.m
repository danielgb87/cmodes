function [hbar herr] = bar_errsig(varargin)
%     [hbar herr] = bar_errsig(varargin)
%     Bar plot with error bars and significance "stars"
%     Uses the plotwitherr function (by Martina F. Callaghan, available from Mathworks)
%     FORMS
%         [hbar herr] = bar_errsig([],Xstd,X)
%         [hbar herr] = bar_errsig(h,Xstd,X)
%         [hbar herr] = bar_errsig(h,Xstd,barvars)
%         [hbar herr] = bar_errsig(h,Xstd,bar_vars,h_options)
%     INTPUTS
%         h - optional - vector or matrix of hypotheses size(h,1) must
%           equal size(X,1); 
%         Xstd - vector for error bars
%         X - vector or matrix of data
%         bar_vars - variables as passed to conventional bar plot command
%         h_options - additional name value pairs listed below are
%                     available to bar_errsig
%     H_OPTIONS
%         'mysymbol' - symbol to use for denoting significance
%         'myfontsize' - Font size of symbol (default 25)
%         'star_shift_factor' - amount by which to shift the significance "star" 
%     OUTPUTS
%         [HL] - plot handle
% 
%     EXAMPLES
%     Example 1 - 1D barplot, 1d h matrix
%     y = randn(1,5);         
%     errY = 0.1.*y;          % 10% error
%     h = [0 1 0 0 1];
%     h = bar_errsig(h,errY, y,'g');% Plot with errorbars
% 
%     Example 2 - 2D barplot, 1d h matrix, with some optional arguments
%     y = randn(3,4);        
%     errY = 0.1.*y;          % 10% error
%     h = [0 1 0];
%     h = bar_errsig(h,errY, y,'mysymbol','#','star_shift_factor',1.2);% Plot with errorbars
% 
%     Example 3 - 2D barplot, 2d h matrix
%     y = randn(3,4);        
%     errY = 0.1.*y;          % 10% error
%     h = [0 1 0 1; 0 1 1 1; 1 0 1 0];
%     out = bar_errsig(h,errY, y);% Plot with errorbars
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    [mysymbol, varargin] = parse_pair(varargin,'mysymbol','*');
    [star_shift_factor, varargin] = parse_pair(varargin,'star_shift_factor',1.1);
    [myfontsize, varargin] = parse_pair(varargin,'myfontsize',25);

    % Same inputs as barwitherr, except 1st argument is 
    % an array of h values (1 or 0) indicating whether
    % to put a star
    
    % Example usage
    % [h1 h2 ] = bar_errsig(any(h,1),Xste,Xmu);
    
    h=varargin{1};
    Xstd = varargin{2};
    [hbar herr ] =barwitherr(varargin{2:end});
    
    Xmu = [];
    for i = 1:length(hbar)
        Xmu = [Xmu; get(hbar(i),'YData')];
    end
    Xmu = Xmu';
    
    if isempty(h); return; end  % If hypothesis vector is empty, return
    
    if isvector(h); h = h(:); end
    
    xdata=get(get(hbar(end),'children'),'xdata');
    xlocs = mean(xdata,1);
%     xlocs2D = unique(xdata,'rows');   % This method for calculating xlocs2D is wrong
%     xlocs2D = xlocs2D';
    
    xlocs2D = [];
    for i = 1:length(hbar)
        xdata=get(get(hbar(i),'children'),'xdata');
        xlocs2D = [xlocs2D; mean(xdata,1)];
    end
    xlocs2D=xlocs2D';
        
%     This should be fine, but it was causing a lot of false errors, so I'm commmenting it out
%     if ~isvector(Xmu) && (length(xlocs) ~= size(h,1))
%         fprintf('Error, size(X,3) of input data X must equal the number of hypotheses (h) provided \n.');
%         return
%     end
    
    % 1 d significance vector
    if isvector(h)
        for i = 1:length(h)
            % old code - should really specify a 2D h vector
            if ~isvector(Xmu)
                if h(i)
                    text(mean(xlocs2D(i,:)),star_shift_factor*(max(Xmu(i,:)+Xstd(i,:))),mysymbol,'FontSize',myfontsize,'HorizontalAlignment','center');
                end
            else
                if h(i)
                    text(xlocs(i),star_shift_factor*(Xmu(i)+Xstd(i)),mysymbol,'FontSize',myfontsize,'HorizontalAlignment','center');
                end

            end
        end
    else
        % 2d significance vector
        sz1 =size(h); sz2 = size(Xmu);
        if sz1 == sz2
            for i = 1:sz1(1)
                for j = 1:sz1(2)
                    if h(i,j)
                        text(xlocs2D(i,j),star_shift_factor*max(Xmu(i,:)+Xstd(i,:)),mysymbol,'FontSize',myfontsize,'HorizontalAlignment','center');
                    end
                end
            end
        end
    end
    
    
end