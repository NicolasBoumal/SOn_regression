function hcn = plot_hessian_condition_number(optim_problem, X)
% Evaluate the condition number of the Hessian of the cost function at a
% given point X. This is typically interesting when X is a critical point
% the optimization algorithm reached. A large number indicates a poorly
% conditioned Hessian, which typically results in slow first-order methods.
%
% This can take a while to run if optim_problem.M.dim() is more than a
% couple hundreds
%
% Nicolas Boumal, Oct. 2017.
    
    evs = hessianspectrum(optim_problem, X);
    
    hcn = evs(1)/evs(end);
    
    figure;
    stairs(log10(evs));
    xlabel('Index k');
    ylabel('log_{10}(\lambda_k)');
    title(sprintf('Condition number: %g', hcn));
    drawnow;

end
