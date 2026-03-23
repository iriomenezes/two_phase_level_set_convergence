function plot_contour_radial(C,plots,linestyle,filename)
% Plots the contours described in C.
%#ok<*AGROW>
%#ok<*NASGU>

% Default contour for testing.
if nargin < 1
    C = [0 1 0.5 0 ; 3 0 0.8 1];
end

% Default axis handle.
if nargin < 2
    plots = 1;
end

% Default linestyle.
if nargin < 3
    linestyle = 'b-';
end

% Default filename.
if nargin < 4
    filename = [];
end

% Check for quarter domain.
q = max(abs(C(:,2) - C(:,1+C(2,1)))) > 1e-15;

% Setup figure.
if plots ~= 0
    if plots == 1
        figure, h1 = gca;
    else
        h1 = plots;
    end
    hold(h1,'on')
    xlabel(h1,'x')
    ylabel(h1,'y')
    box(h1,'on')
end

K = size(C,2);
k = 0;
c = 0;
while k < K
    k = k + 1;
    kk = C(2,k);
    x = C(1,k+1:k+kk);
    y = C(2,k+1:k+kk);
    % Check for consistent ordering.
    if q == 1
        if y(1) ~= 0
            x = x(end:-1:1);
            y = y(end:-1:1);
        end
    end
    % Add symmetric contour.
    X = x;
    Y = y;
    if q == 1
        X = [X -x(end-1:-1:1)];
        Y = [Y y(end-1:-1:1)];
        X = [X X(end-1:-1:1)];
        Y = [Y -Y(end-1:-1:1)];
    end
    % Plot contour.
    if plots ~= 0
        plot(h1,X,Y,linestyle)
    end
    % Save contour.
    if ~isempty(filename)
        A = [X' Y'];
        str = [filename sprintf('_%.0f.txt',c)];
        save(str,'A','-ascii')
    end
    k = k + kk;
    c = c + 1;
end





























