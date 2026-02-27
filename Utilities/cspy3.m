function [] = cspy3(matrix,varargin)
%CSPY Visualize sparsity pattern.
%   CSPY(S) plots the sparsity pattern of the matrix S with L levels 
%
%   CSPY(S,'Marker', '*') use the this marker to plot the matrix S.
%
%   CSPY(S, 'Marker', {'*', '+'}) use each marker per level in the same
%   order. The size of the marker cell array must be equal to the number of
%   levels, in otherwise is used the first marker for all levels.
%
%   CSPY(S,'MarkerSize', M) use the M markersize to plot the matrix S.
%
%   CSPY(S, 'Levels', L) use L levels to show different colors in the
%   matrix S. If the levels is a vector and has the same size to channels
%   is used one value for channel. Where level hasn't the same size is used
%   the mean value.
%   
%   CSPY(S, 'XDir', 'Reverse') change the X direction (Reverse or Normal).
%   Default Normal
%
%   CSPY(S, 'YDir', 'Reverse') change the Y direction (Reverse or Normal). 
%   Default Reverse
%   
%
%   NOTES:
%       - set up an odd number of levels, for a small even number the color 
%       map will be represented wrongly
%
%    the original code of :
%   Hugo Gualdron - gualdron@usp.br
%   Jose Rodrigues - junio@icmc.usp.br
%   University of São Paulo
%   https://www.mathworks.com/matlabcentral/fileexchange/46551-cspy-m
%   
% modified by Artem Mavliutov, 21/07/2023
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default values
    islevels = 0;
    levels = 0;
    MAP = [[1 0 0];[0.90 0.90 0.90];[0 0 1]];
    markerSize = 8;
    marker = '.';
    xdir = 'normal';
    ydir = 'reverse';
    % input
    for i=1:2:nargin-1
        switch(lower(char(varargin(i))))
            case 'marker'
                marker = varargin(i+1);
                if iscell(marker)
                    marker = marker{:};
                end
            case 'markersize'
                markerSize = cell2mat(varargin(i+1));
            case {'levels','lvl','lvls','lev','level','l'}
                levels = cell2mat(varargin(i+1));
                islevels = 1;
            case 'colormap'
                MAP = cell2mat(varargin(i+1));
            case 'ydir'
                ydir = lower(char(varargin(i+1)));
            case 'xdir'
                xdir = lower(char(varargin(i+1)));
            otherwise
                warning('Unexpected input argument')
        end
    end
    %------------------------------------------------------------------------- 
    figure;

    % Set the figure position to full screen with fixed ratio
    screenSize = get(0, 'Screensize');
    set(gcf, 'Position', [screenSize(3)/2-screenSize(4)/2,screenSize(4)/2, screenSize(4), screenSize(4)]);

    [y, x, z] = find(matrix);
    minvalue = min(z);
    maxvalue = max(z);
    
    %number of levels if not defined
    if ~islevels 
        levels = 3;
        % levels = max(ceil(max(abs(maxvalue),abs(minvalue)) ...
        %     / min(abs(maxvalue),abs(minvalue)) ),2);
    % take the odd number of levels
    elseif mod(levels,2) == 0
        levels = levels + 1;
    end
    
    step = abs(maxvalue - minvalue)/levels;
    % step = max(abs(minvalue),abs(maxvalue)) / levels;
    colors = flipud(redgreenblue(levels,MAP));

    % centralize the values around zero (usefull when the min and max are
    % distanced differently from zero), see also clim
    minvalue = 0-step/2 - (levels-1) / 2 * step; % rem. pseudo minvalue
    limit = max(abs(minvalue),abs(maxvalue));
    
    colormap(colors);
    % if size(matrix,1) == size(matrix, 2)
    %     axis square;
    % end
    if levels > 1
        cb = colorbar;
        % caxis([minvalue maxvalue]);
        clim([sign(minvalue)*limit sign(maxvalue)*limit]);
        if levels <= 10
            cb.Ticks = minvalue:step:limit;
        end
    end
    xlim([1, size(matrix, 2)]);
    ylim([1, size(matrix, 1)]);
    title(['nz = ' num2str(nnz(matrix))  '    scal = 1e' num2str(log10(limit),1)])
    hold on;

    set(gca,'XDir', xdir);
    set(gca,'YDir', ydir);
    set(gca,'FontSize',12)
    % set(gca, 'XAxisLocation', 'top');
    
    % permutation from ceter out is needed in order to plot the zero values
    % at first and then add the other values increasingly.
    permlevels = middleout(levels);
    for k=1:levels
        i = permlevels(k);
        step_init = minvalue + (i-1)*step;
        step_end = minvalue + i*step;
        if i == levels
            ids = find(z>=step_init);
        elseif i == 1
            ids = find(z<step_end);
        else
            ids = find(z>=step_init & z<step_end);
        end

        % markerSize = rescale(z(ids),1,100); % too slow

        scatter3(x(ids), y(ids), z(ids), markerSize, colors(i,:),"filled");

        % curtick = get(gca, 'XTick');
        % set(gca, 'XTickLabel', cellstr(num2str(curtick(:)))); 
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function arr = middleout(n)
    % find a permutation array that access a general array of size n 
    % from center out
    % Artem Mavliutov, 21/07/2023
    arr = zeros(n,1);
    i = ceil(n / 2 );
    j = i + 1;
    ind = 1;
    while (j <= n)
        if (mod(n,2) == 0)
            arr(ind) = j; ind = ind + 1;
            arr(ind) = i; ind = ind + 1;
        else
            arr(ind) = i; ind = ind + 1;
            arr(ind) = j; ind = ind + 1;
        end
        j = j + 1;
        i = i - 1;
    end
    if (i == 1), arr(ind) = i; end
end

function c = redgreenblue(M,varargin)
    % given the input of 3 colors define the intermediate shades
    % returns the color map 'c' of size Mx3
    % the default color values are red green and blue
    % Artem Mavliutov, 21/07/2023
    if nargin < 1, M = size(get(gcf,'colormap'),1); end
    
    if nargin < 2
        % default rgb color values
        r = [1 0 0];
        g = [1 1 1];
        b = [0 0 1];
    else
        % custom rgb color values
        r = varargin{1}(1,:);
        g = varargin{1}(2,:);
        b = varargin{1}(3,:);
    end
    
    if (mod(M,2) == 0)
        M1 = M*0.5;
        if M == 2
            c1 = [r(1);g(1)];
            c2 = [r(2);g(2)];
            c3 = [r(3);g(3)];
        else
            c1 = [linspace(r(1),g(1),M1)';linspace(g(1),b(1),M1)'];
            c2 = [linspace(r(2),g(2),M1)';linspace(g(2),b(2),M1)'];
            c3 = [linspace(r(3),g(3),M1)';linspace(g(3),b(3),M1)'];
        end
    else
        M1 = floor(M*0.5);
        if M == 1
            c1 = g(1);
            c2 = g(2);
            c3 = g(3);
        else
            c1 = linspace(g(1),b(1),M1+1)';
            c2 = linspace(g(2),b(2),M1+1)';
            c3 = linspace(g(3),b(3),M1+1)';
    
            c1 = [linspace(r(1),g(1),M1+1)';c1(2:M1+1)];
            c2 = [linspace(r(2),g(2),M1+1)';c2(2:M1+1)];
            c3 = [linspace(r(3),g(3),M1+1)';c3(2:M1+1)];
        end
    end
    c = [c1 c2 c3];
end