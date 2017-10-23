function [X, info, optim_problem] = digress(problem, X0, options)
% DIGRESS algorithm for discrete regression on SO(n).
%
% See the paper (PDF included on the GIT repository):
%
% Interpolation and Regression of Rotation Matrices,
% Nicolas Boumal, 2013
% In: Nielsen F., Barbaresco F. (eds) Geometric Science of Information.
% Lecture Notes in Computer Science, vol 8085. Springer, Berlin, Heidelberg
% https://link.springer.com/chapter/10.1007/978-3-642-40020-9_37
% 
% Code was adapted from the original in October 2017.
%
% Input:
%
%   problem
%       Structure describing the regression problem. Fields are as defined
%       in the example script main.m.
%
%   X0
%       Initial guess for the discrete curve: a 3D matrix
%       of size n x n x Nd whose Nd slices are rotation matrices. If this
%       is omitted or empty ([]), one is created automatically.
%
%   options
%       Structure passed to the Manopt solver trustregions. This is
%       optional. Values passed in this structure overwrite all defaults.
%       Parameters of particular interest to fix include:
%           options.tolgradnorm (stopping criterion for gradient norm)
%           options.verbosity (set to 0 to silence solver)
%           options.maxinner (to limit Hessian calls per outer iteration)
%       See Manopt.org for more information.
%
%
% The requires Manopt, freely available at http://www.manopt.org.
%
% Nicolas Boumal, Oct. 2017.

    n = problem.n;
    Nd = problem.Nd;
    
    % If no initial guess is provided, create one.
    if ~exist('X0', 'var') || isempty(X0)
        X0 = initguess(problem);
    end
    
    % Options for the optimization algorithm are optional.
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
   
    % Obtain Manopt description of the manifold SO(n)^Nd.
    manifold = rotationsfactory(n, Nd);

    % Define the cost function, at X.
    function [E, store] = mycost(X, store)
        if ~isfield(store, 'E')
            [store, E] = cost(store, problem, X);
            store.E = E;
        end
        E = store.E;
    end

    % Define the gradient of the cost, at X.
    function [gradE, store] = mygrad(X, store)
        if ~isfield(store, 'gradE')
            [store, ~, gradE] = cost(store, problem, X);
            store.gradE = manifold.egrad2rgrad(X, gradE);
        end
        gradE = store.gradE;
    end

    % Define the Hessian of the cost, at X, along XO.
    function [hessE, store] = myhess(X, O, store)
        [store, ~, gradE, hessE] = cost(store, problem, X, O);
        hessE = manifold.ehess2rhess(X, gradE, hessE, O);
    end
    
    % Combine the manifold, the cost and its derivatives into a problem
    % structure for Manopt.
    optim_problem.M = manifold;
    optim_problem.cost = @mycost;
    optim_problem.grad = @mygrad;
    optim_problem.hess = @myhess;
    
    % If you decide to change the cost function and want to check that the
    % gradient and Hessian are still correct, see checkgradient and
    % checkhessian tools in Manopt.
    
    % Options for the optimization algorithm. These will be overwritten by
    % user supplied options, if any.
    options_default.maxinner = 3*manifold.dim();
    options_default.tolgradnorm = 1e-5;
    options = mergeOptions(options_default, options);
    
    % Call the solver.
    [X, ~, info] = trustregions(optim_problem, X0, options);

end
