function X = initguess(problem)
% Produce an initial curve for discrete regression on SO(n).
%
% The initial curve X is based on piecewise geodesic interpolation.
%
% Input 'problem' is a structure defining a regression problem on SO(n):
% see digress.
%
% X is a 3D matrix of size n x n n Nd, as defined by problem.n and
% problem.Nd. Each slice X(:, :, k) is in SO(n) (orthogonal, determinant
% +1).
%
% This function assumes s (as defined in problem.s) satisfies: s(1) = 1,
% s(end) = Nd (as defined in problem.Nd) and s(k+1) > s(k). If this is not
% the case, the returned curve is still valid, but may be a poor
% initializer for the optimization phase in digress.
% 
% By construction, X(:, :, s(k)) == p(k), where p is defined in problem.p
% and k = 1:N, with N as defined in problem.N.
%
% See also: digress
%
% Nicolas Boumal, Oct. 2017.

    n = problem.n;
    N = problem.N;
    p = problem.p;
    Nd = problem.Nd;
    s = problem.s;
    
    % Obtain a Manopt representation of the rotation group
    SOn = rotationsfactory(n);
    
    X = zeros(n, n, Nd);
    
    % Make sure all matrices are in SO(n) to start.
    for k = 1 : Nd
        X(:, :, k) = eye(n);
    end
    
    for k = 1 : N-1
        
        R = p(:, :, k);
        H = SOn.log(p(:, :, k), p(:, :, k+1));
        
        X(:, :, s(k)) = p(:, :, k);
        
        % Sample geodesic from p(:, :, k) to p(: ,:, k+1)
        steps = s(k+1) - s(k);
        for j = 1 : steps-1
            X(:, :, s(k)+j) = SOn.exp(R, H, j/steps);
        end
        
    end
    
    X(:, :, s(end)) = p(:, :, end);

end
