function model = translate_model(model, t)
% Translate a model by t (a vector of length three, either row or col.)
% Nicolas Boumal, Oct. 2017

    model.v = model.v + repmat(t(:)', size(model.v, 1), 1);

end
