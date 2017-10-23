function plotSO3curve(X, delta_tau, save_plots)
% Render a discrete curve on SO(3) as a movie.
%
% X is a matrix of size 3x3xM such that each slice X(:, :, k) is a rotation
% matrix (orthogonal, determinant +1).
%
% delta_tau > 0 is the time separation between two slices of X.
%
% By default, save_plots = []. If it is a string, then frames are saved to
% frames/frame_"save_plots"_#.png, where # is the frame number and
% "save_plots" is the given string. Code should be adapted to synchronize
% videos with different delta_tau's and different M's.
%
% Nicolas Boumal, Oct. 2017

    assert(size(X, 1) == 3 && size(X, 2) == 3, ...
           'Can only show rotations in 3D: X must be 3x3xM');
       
    M = size(X, 3);

    if ~exist('save_plots', 'var') || ~ischar(save_plots)
        save_plots = [];
    end

    [Xs, Ys, Zs] = ellipsoid(0, 0, 0, 1, .5, .5, 30);

    % Alternative to ellipsoid: a cube.
    % vert = [1 1 1; 1 2 1; 2 2 1; 2 1 1 ; ... 
    %         1 1 2; 1 2 2; 2 2 2; 2 1 2]; 
    % fac = [1 2 3 4; ... 
    %        2 6 7 3; ... 
    %        4 3 7 8; ... 
    %        1 5 8 4; ... 
    %        1 2 6 5; ... 
    %        5 6 7 8];

    % SO(3) geodesic between A at time t = 0 and B at time t = 1.
    f = @(A, B, t) real(A*expm(t*logm(A.'*B)));

    frameno = 1;

    % Desired number of frames: this could be adapted to synchronize with
    % respect to delta_tau.
    framecount = 2*(M-1)+1;
    
    tau = delta_tau*(0:(M-1));

    for t = linspace(tau(1), tau(end), framecount)

        i = find(tau > t, 1, 'first');
        if isempty(i)
            i = M;
        end

        C = f(X(:, :, i-1), X(:, :, i), (t-tau(i-1))/delta_tau);

        S = C*[ Xs(:).' ; Ys(:).' ; Zs(:).' ];
        Xr = reshape(S(1,:), size(Xs));
        Yr = reshape(S(2,:), size(Ys));
        Zr = reshape(S(3,:), size(Zs));

        % If working with cube:
        % vertrotated = (C*vert.').';

        clf;
        hold on;
        
        % Plot ellipsoid
        surf(Xr, Yr, Zr, 'FaceColor', [203 220 237]/255);
        
        % If working with cube:
        % patch('Faces', fac, 'Vertices', vertrotated, 'FaceColor', [.8 .4 0]);
        
        % Plot axes
        r = 1.1;
        plot3(r*[1 -1], [0 0], [0 0], '-', 'LineWidth', 2, 'Color', [42 83 125]/255);
        plot3([0 0], r*[1 -1], [0 0], '-', 'LineWidth', 2, 'Color', [42 83 125]/255);
        plot3([0 0], [0 0], r*[1 -1], '-', 'LineWidth', 2, 'Color', [42 83 125]/255);
        
        view(25, 29);
        camzoom(2);
        
        hold off;

        set(gcf, 'color', 'w');
        
        axis equal;
        axis off;
        
        drawnow;

        if ~isempty(save_plots)
            fname = sprintf('frames/frame_%s_%05d', save_plots, frameno);
            print('-dpng', '-opengl', '-r120', [fname '.png']);
        end
        frameno = frameno+1;

    end
