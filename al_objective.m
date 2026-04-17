function L = al_objective(x, p_struct, lambda, p_weight)
    % Calculate base weight
    W = beam_weight(x, p_struct);
    
    % Get constraint violations
    [g, ~] = beam_constraints(x, p_struct);
    
    % Add Augmented Lagrangian terms (only for violated constraints)
    penalty_sum = 0;
    for i = 1:length(g)
        viol = max(0, g(i)); 
        penalty_sum = penalty_sum + lambda(i)*viol + p_weight*(viol^2);
    end
    
    L = W + penalty_sum;
end