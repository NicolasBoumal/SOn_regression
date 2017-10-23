function [store, E, gradE, hessE] = cost(store, problem, X, O)
% Cost function with gradient and Hessian for DIGRESS --- see paper.
%
% store is a structure used for caching. problem is regression problem
% structure (/not/ a Manopt problem structure.)
%
% Nicolas Boumal, Oct. 2017

    needhess = (nargin == 4);

    gradE = zeros(size(X));
    if needhess
        hessE = zeros(size(X));
    end
    
    %% Data fidelity term

    s = problem.s;
    w = problem.w;
    p = problem.p;
    
    if ~isfield(store, 'memory_df')
        store.memory_df = struct();
    end
    memory_df = store.memory_df;
    
    if needhess
        [memory_df, df_vals, df_grads, ~, df_hess] = ...
                                funf(memory_df, X(:, :, s), p, O(:, :, s));
        hessE(:, :, s) = hessE(:, :, s) + multiscale(w, df_hess);
    else
        [memory_df, df_vals, df_grads] = funf(memory_df, X(:, :, s), p);
    end
    store.memory_df = memory_df;
    
    data_fidelity = w' * df_vals;
    
    gradE(:, :, s) = gradE(:, :, s) + multiscale(w, df_grads);
    
    
    %% Velocity regularity term
    lambda = problem.lambda;
    Nd = problem.Nd;
    delta_tau = problem.delta_tau;
    
    if lambda > 0

        I = 1:Nd-1;
        velocity_coeff = lambda/delta_tau;

        if ~isfield(store, 'memory_vr')
            store.memory_vr = struct();
        end
        memory_vr = store.memory_vr;

        if needhess
            [memory_vr, vr_vals, vr_grads_A, vr_grads_B, vr_hess_A, vr_hess_B] = ...
                      funf(memory_vr, X(:, :, I), X(:, :, I+1), O(:, :, I), O(:, :, I+1));
            hessE(:, :, I)   = hessE(:, :, I)   + velocity_coeff*vr_hess_A;
            hessE(:, :, I+1) = hessE(:, :, I+1) + velocity_coeff*vr_hess_B;

        else
            [memory_vr, vr_vals, vr_grads_A, vr_grads_B] = funf(memory_vr, X(:, :, I), X(:, :, I+1));
        end
        store.memory_vr = memory_vr;

        velocity_reg = velocity_coeff * sum(vr_vals);

        gradE(:, :, I)   = gradE(:, :, I)   + velocity_coeff*vr_grads_A;
        gradE(:, :, I+1) = gradE(:, :, I+1) + velocity_coeff*vr_grads_B;
        
    else
        
        velocity_reg = 0;
        
    end
    
    
    %% Acceleration regularity term
    
    mu = problem.mu;
    
    if mu > 0
    
        I = 2:Nd-1;
        accel_coeff = mu/delta_tau^3;

        if ~isfield(store, 'memory_ar')
            store.memory_ar = struct();
        end
        memory_ar = store.memory_ar;

        if needhess
            [memory_ar, ar_vals, ar_grads_A, ar_grads_B, ar_grads_C, ...
                      ar_hess_A,  ar_hess_B,  ar_hess_C ] = ...
                               fung(memory_ar, X(:, :, I), X(:, :, I+1), X(:, :, I-1), ...
                                    O(:, :, I), O(:, :, I+1), O(:, :, I-1));
            hessE(:, :, I)   = hessE(:, :, I)   + accel_coeff*ar_hess_A;
            hessE(:, :, I+1) = hessE(:, :, I+1) + accel_coeff*ar_hess_B;
            hessE(:, :, I-1) = hessE(:, :, I-1) + accel_coeff*ar_hess_C;
        else
            [memory_ar, ar_vals, ar_grads_A, ar_grads_B, ar_grads_C] = ...
                                  fung(memory_ar, X(:, :, I), X(:, :, I+1), X(:, :, I-1));
        end
        store.memory_ar = memory_ar;

        accel_reg = accel_coeff * sum(ar_vals);

        gradE(:, :, I)   = gradE(:, :, I)   + accel_coeff*ar_grads_A;
        gradE(:, :, I+1) = gradE(:, :, I+1) + accel_coeff*ar_grads_B;
        gradE(:, :, I-1) = gradE(:, :, I-1) + accel_coeff*ar_grads_C;
        
    else
        
        accel_reg = 0;
        
    end

    E = data_fidelity + velocity_reg + accel_reg;
    
end
