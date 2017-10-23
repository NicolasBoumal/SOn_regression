function p = plotrotations(R)
% Display a set of rotations as a sequence of rotated copies of a 3D object.
%
% R is a 3D matrix of size 3x3xM such that each slice R(:, :, k) is a
% rotation matrix (orthogonal, determinant +1).
%
% In the currently active figure, this function will render a 3D object
% M times on a horizontal line, such that the furthest left has been
% rotated by R(:, :, 1) and the furthest right by R(:, :, M).
%
% Nicolas Boumal, Oct. 2017.
    
    % Read a 3D model in STL format
    try
        [model.v, model.f, model.n, model.c] = stlread('gorilla.stl');
    catch
        error(sprintf(['Please add the STL directory to your ' ...
               'Matlab path:\naddpath(''STL'');'])); %#ok<SPERR>
    end
    
    % Shift center of mass to origin
    model = center_model(model);
    
    % Horizontal spacing between two copies of the object.
    D = 90;

    % Create a larger model, which comprises size(R, 3) copies of the
    % model, rotated and translated.
    full_model = rotate_model(model, R(:, :, 1));
    for k = 2 : size(R, 3)
        new_model = rotate_model(model, R(:, :, k));
        new_model = translate_model(new_model, (k-1)*[D, 0, 0]);
        full_model = add_models(full_model, new_model);
    end
    
    % Render
    p = show_model(full_model);
    
    % Change color and lighting
    p.LineStyle = 'none';
    p.FaceLighting = 'flat';
    p.FaceColor = [1 1 1];

    % Place a light at some point (this is not ideally defined...)
    view(-38, 6);
    light;
    view(0,0);
    
    % Try to force Matlab to use all the available space for rendering.
    set(gca,'Position', [0 0 1 1]);
    
    axis equal;
    axis off;
    set(gcf, 'color', 'k');
    
end
