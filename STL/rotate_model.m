function model = rotate_model(model, R)
% Apply rotation matrix R to all vertices of a model.
%
% Current code does not update face normals. This is inconsequential for
% rendering since Matlab's patch function computes normals itself.
%
% Nicolas Boumal, Oct. 2017

    model.v = model.v * R';
    
end
