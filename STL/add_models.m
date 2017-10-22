function model = add_models(model1, model2)
% Add two models so they can be rendered together.
% Nicolas Boumal, Oct. 2017

    model.v = [ model1.v ; model2.v ];
    model.f = [ model1.f ; model2.f + size(model1.v, 1) ];
    model.n = [ model1.n ; model2.n ];
    model.c = [ model1.c ; model2.c ];

end
