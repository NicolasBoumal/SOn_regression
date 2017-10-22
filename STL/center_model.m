function model = center_model(model)
% Center a model so its center of mass is at the origin.
% Nicolas Boumal, Oct. 2017

    center = mean(model.v, 1);
    n = size(model.v, 1);
    model.v = model.v - repmat(center, [n, 1]);
    
end
