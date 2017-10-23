function [problem, Y] = refine(problem, X)
% Refines a discrete curve on SO(n) by adding one point between every two
% points, sampled as the midpoint of the geodesic joining them.
%
% Thus, an input curve with Nd points is returned with 2*Nd - 1 points.
%
% Nicolas Boumal, Oct. 2017.

    n = problem.n;
    Nd = problem.Nd;
    
    assert(Nd >= 2, 'To refine, need Nd >= 2.');
    assert(ndims(X) == 3 && all(size(X) == [n, n, Nd]), ...
                                     'X must be n x n x Nd with Nd >= 2.');

    Y = zeros(n, n, 2*Nd-1);
    
    % Every other point is simply copied from X
    Y(:, :, 1:2:end) = X;
    
    % The remaining points are obtained along the geodesics of SO(n).
    for k = 1 : (Nd-1)
        prev = k;
        next = k+1;
        L = real(logm(X(:, :, prev)'*X(:, :, next)));
        Y(:, :, 2*k) = X(:, :, prev)*expm(.5*L);
    end
    
    % Update problem structure
    problem.Nd = 2*problem.Nd-1;
    problem.s = 2*problem.s-1;
    problem.delta_tau = problem.delta_tau/2;

end
