%% Discrete regression curves on rotation group SO(n)
%
% Example script for code accompanying the paper:
%
% https://link.springer.com/chapter/10.1007/978-3-642-40020-9_37
%
% SO(n) is the set of orthogonal matrices of size n x n and determinant +1.
%
% This code requires Manopt, freely available at http://www.manopt.org.
%
% Nicolas Boumal, Oct. 2017

clear;
close all;
clc;

%% Define a regression problem by defining N control points on SO(n)

% Example 1: explicit construction
% n = 3;
% N = 4;
% p = zeros(n, n, 4);
% p(:, :, 1) = eye(n);
% p(:, :, 2) = [0 1 0 ; -1 0 0 ;  0  0  1];
% p(:, :, 3) = [0 0 1 ; -1 0 0 ;  0 -1  0];
% p(:, :, 4) = [0 0 1 ;  0 1 0 ; -1  0  0];

% Example 2: load from mat file
data = load('controlpoints.mat');
n = data.n;
N = data.N;
p = data.p;

% For each control point, pick a weight (positive number). A larger value
% means the regression curve will pass closer to that control point.
w = ones(N, 1);

%% Define parameters of the discrete regression curve

% The curve has Nd points on SO(n)
Nd = 97;

% Each control point attracts one particular point of the regression curve.
% Specifically, control point k (in 1:N) attracts curve point s(k).
% The vector s of length N usually satsifies:
% s(1) = 1, s(end) = Nd and s(k+1) > s(k).
s = round(linspace(1, Nd, N));

% Time interval between two discretization points of the regression curve.
% This is only used to fix a scaling. It is useful in particular so that
% other parameter values such as w, lambda and mu (see below) have the same
% sense even when the discretization parameter Nd is changed.
delta_tau = 1/(Nd-1);

% Weight of the velocity regularization term (nonnegative). The larger it
% is, the more velocity along the discrete curve is penalized. A large
% value usually results in a shorter curve.
lambda = 0;

% Weight of the acceleration regularization term (nonnegative). The larger
% it is, the more acceleration along the discrete curve is penalized. A
% large value usually results is a 'straighter' curve (closer to a
% geodesic.)
mu = 1e-2;

%% Pack all data defining the regression problem in a problem structure.
problem.n = n;
problem.N = N;
problem.Nd = Nd;
problem.p = p;
problem.s = s;
problem.w = w;
problem.delta_tau = delta_tau;
problem.lambda = lambda;
problem.mu = mu;

%% Call the optimization procedure to compute the regression curve.

% Compute an initial guess for the curve. If this step is omitted, digress
% (below) will compute one itself. X0 is a 3D matrix of size n x n x Nd,
% such that each slice X0(:, :, k) is a rotation matrix.
%
X0 = initguess(problem);

% Run the optimization procedure to compute X1, the discrete regression
% curve. X1 is a 3D matrix of size n x n x Nd with each slice a rotation
% matrix. The second output, info, is a struct-array containing information
% returned by the optimization algorithm. The third output, optim_problem,
% is the Manopt optimization problem structure used to produce X1. It can
% be used to run another algorithm, e.g., for research purposes.
%
[X1, info, optim_problem] = digress(problem, X0);

%% Compare Initial guess and optimized curve

% We only plot a subset of the regression curve. The numbers are such that
% the control points attract points that are represented. Specifically, the
% value in s appear in 1:8:97 as rotation 1, 5, 9, 13. In the case of X0,
% these rotations are the control points, by construction.

figure(1);
plotrotations(X0(:, :, 1:8:Nd));
view(0, 0);

figure(2);
plotrotations(X1(:, :, 1:8:Nd));
view(0, 0);

%% Plot optimization information

figure(3);
semilogy([info.time], [info.gradnorm], 'k.-');
title('Gradient norm');
xlabel('Computation time [s]');
pbaspect([1.6, 1, 1]);

% Can also display information about Hessian at X1:
% plot_hessian_condition_number(optim_problem, X1);

%% Plot speed and acceleration of X0 and X1

[speed0, acc0] = compute_profiles(problem, X0);
[speed1, acc1] = compute_profiles(problem, X1);

% Passage time of each point on the discrete curves.
time = problem.delta_tau*( 0 : (problem.Nd-1) );

figure(4);

subplot(1, 2, 1);
plot(time, speed0, time, speed1);
title('Speed of initial curve and optimized curve');
xlabel('Time');
ylabel('Speed');
legend('Initial curve', 'Optimized curve', 'Location', 'SouthEast');
pbaspect([1.6, 1, 1]);

subplot(1, 2, 2);
plot(time, acc0, time, acc1);
title('Acceleration of initial curve and optimized curve');
xlabel('Time');
ylabel('Acceleration');
legend('Initial curve', 'Optimized curve', 'Location', 'NorthWest');
pbaspect([1.6, 1, 1]);

ylim([0, 20]);

%% Refine a regression curve

% Especially for large values of mu, the optimization problem can be
% poorly conditioned, and as a result take a long time to solve. One way to
% alleviate this is to proceed in stages: compute a regression curve with
% small Nd, then refine that curve through piecewise geodesic interpolation
% and reoptimize. This procedure can be iterated.

compute_refinement = false;
if compute_refinement

    [problem_refined, X1r] = refine(problem, X1);

    X2 = digress(problem_refined, X1r);

end

%% Display the whole regression curves as movies

plot_movies = false;
if plot_movies
    
    plotSO3curve(X0, problem.delta_tau);
    pause;
    
    plotSO3curve(X1, problem.delta_tau);
    pause;
    
    % The refined curve appears slower because the movie must produce more
    % frames, but the frames can be synchronized when outputing a video
    % file.
    if compute_refinement
        plotSO3curve(X2, problem_refined.delta_tau);
    end
    
end
