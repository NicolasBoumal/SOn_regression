function [speed, acc] = compute_profiles(problem, X)
% Compute speed and acceleration along discrete curve on SO(n).
%
% Speed and acceleration are defined based on finite difference
% approximations on the manifold SO(n) (the rotation group.)
%
% X is a 3D matrix of size n x n x Nd.
% problem is a problem structure for a regression problem.
% speed and acc are vectors of length Nd such that speed(k) and acc(k)
% represent the speed and acceleration of X at X(:, :, k) on SO(n).
% Acceleration is NaN at the first and last points.
%
% Nicolas Boumal, Oct. 2017.

    n = problem.n;
    Nd = problem.Nd;
    dtau = problem.delta_tau;
    
    assert(all(size(X) == [n, n, Nd]), ...
      'X must have size n x n x Nd. Check consistency with problem struct.');
    
    % Get a Manopt representation of SO(n) for easy access to geometric
    % operations such as exponential (geodesic), distance and logarithm.
    SOn = rotationsfactory(n);

    % Simple forward difference formula for speed.
    speed = zeros(1, Nd);
    for k = 1 : Nd-1
        fw = SOn.log(X(:, :, k), X(:, :, k+1));
        v = fw/dtau;
        speed(k) = SOn.norm(X(:, :, k), v);
    end
    % Backward difference for last point.
    speed(end) = speed(end-1);
    
    % Acceleration is NaN at first and last point. For all the others,
    % using a symmetric difference formula.
    acc = NaN(1, Nd);
    for k = 2 : Nd-1
        fw = SOn.log(X(:, :, k), X(:, :, k+1));
        bw = SOn.log(X(:, :, k), X(:, :, k-1));
        a = ( fw + bw ) / ( dtau^2 );
        acc(k) = SOn.norm(X(:, :, k), a);
    end

end
