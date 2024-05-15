% Clear Command Window
clc;


% Definition of the maximization optimization
% problem for the estimation of the upper bound
% for the CHSH Inequality (F_2_1)
f_2_1_chsh_inequality_local_bound_optimization_prob = ...
    optimproblem('ObjectiveSense', 'max');

% Definition of the optimization variables/unknowns
% to be optimized in the maximization optimization
% problem for the estimation of the upper bound
% for the CHSH Inequality
a_1_1 = optimvar('a_1_1', 1, 1, 'Type', 'integer');
a_1_2 = optimvar('a_1_2', 1, 1, 'Type', 'integer');
a_2_1 = optimvar('a_2_1', 1, 1, 'Type', 'integer');
a_2_2 = optimvar('a_2_2', 1, 1, 'Type', 'integer');

% Definition of the objective function of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality (F_2_1)
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Objective = a_1_1(1) * a_2_1(1) + a_1_1(1) * a_2_2(1) + ... 
                 a_1_2(1) * a_2_1(1) - a_1_2(1) * a_2_2(1);

% Definition of the constraints
% for the a_1_1 value E [-1,+1]
constraint_1 = a_1_1(1) >= -1;
constraint_2 = a_1_1(1) <= 1;
constraint_3 = a_1_1(1)^2 <= 1;
constraint_4 = a_1_1(1)^2 >= 1;

% Definition of the constraints
% for the a_1_2 value E [-1,+1]
constraint_5 = a_1_2(1) >= -1;
constraint_6 = a_1_2(1) <= 1;
constraint_7 = a_1_2(1)^2 <= 1;
constraint_8 = a_1_2(1)^2 >= 1;

% Definition of the constraints
% for the a_2_1 value E [-1,+1]
constraint_9 = a_2_1(1) >= -1;
constraint_10 = a_2_1(1) <= 1;
constraint_11 = a_2_1(1)^2 <= 1;
constraint_12 = a_2_1(1)^2 >= 1;

% Definition of the constraints
% for the a_2_2 value E [-1,+1]
constraint_13 = a_2_2(1) >= -1;
constraint_14 = a_2_2(1) <= 1;
constraint_15 = a_2_2(1)^2 <= 1;
constraint_16 = a_2_2(1)^2 >= 1;


% Setup of all the constraints
% for the estimation values
% a_1_1, a_1_2, a_2_1, and a_2_2,
% for the maximization optimization
% problem for the estimation of
% the upper bound for the CHSH Inequality (F_2_1)
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons1 = constraint_1;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons2 = constraint_2;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons3 = constraint_3;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons4 = constraint_4;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons5 = constraint_5;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons6 = constraint_6;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons7 = constraint_7;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons8 = constraint_8;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons9 = constraint_9;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons10 = constraint_10;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons11 = constraint_11;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons12 = constraint_12;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons13 = constraint_13;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons14 = constraint_14;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons15 = constraint_15;
f_2_1_chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons16 = constraint_16;

% Show the mathematical formulation of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality (F_2_1)
show(f_2_1_chsh_inequality_local_bound_optimization_prob);

% Solve the maximization optimization problem
% for the estimation of the (classical) local bound
% for the CHSH Inequality (F_2_1)
f_2_1_chsh_inequality_local_bound_optimization_prob_sol = ...
    solve(f_2_1_chsh_inequality_local_bound_optimization_prob);

% Compute the a_1_1, a_1_2, a_2_1, and a_2_2 values
% for the CHSH Inequality (F_2_1) that allow to
% estimate the (classical) local upper bound L
f_2_1_chsh_inequality_local_bound_a_1_1_value = ...
    f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_1_1;
f_2_1_chsh_inequality_local_bound_a_1_2_value = ...
    f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_1_2;
f_2_1_chsh_inequality_local_bound_a_2_1_value = ...
    f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_2_1;
f_2_1_chsh_inequality_local_bound_a_2_2_value = ...
    f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_2_2;

% Compute the individual expectation values
% for the CHSH Inequality (F_2_1) that allow to estimate
% the (classical) local upper bound L
a_1_1 = f_2_1_chsh_inequality_local_bound_a_1_1_value(1);
a_1_2 = f_2_1_chsh_inequality_local_bound_a_1_2_value(1);
a_2_1 = f_2_1_chsh_inequality_local_bound_a_2_1_value(1);
a_2_2 = f_2_1_chsh_inequality_local_bound_a_2_2_value(1);

% Compute the (classical) local
% upper bound L for the CHSH Inequality (F_2_1)
f_2_1_chsh_inequality_local_upper_bound_L = ...
    a_1_1 * a_2_1 + a_1_1 * a_2_2 + ...
    a_1_2 * a_2_1 - a_1_2 * a_2_2;

% Print a blank line
fprintf('\n');

% Print of the mathematical form
% for the CHSH Inequality (F_2_1),
% in the deterministic strategy form
fprintf('CHSH Inequality:\n');
fprintf('  a_1_1 x a_2_1 + a_1_1 x a_2_1 + a_1_2 x a_2_1 - a_1_2 x a_2_2 <= L\n');

% Print a blank line
fprintf('\n');

% Print of the constraints of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality (F_2_1)
fprintf('Such that:\n');
fprintf('  a_1_1 E {-1, +1}\n');
fprintf('  a_1_2 E {-1, +1}\n');
fprintf('  a_2_1 E {-1, +1}\n');
fprintf('  a_2_2 E {-1, +1}\n');


% Print a blank line
fprintf('\n');

% Print a blank line
fprintf('\n');

% Print of the complete solution of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality (F_2_1)
fprintf(['(Classical) Local Upper Bound L ' ...
         'for CHSH Inequality:\n']);
fprintf('  a_1_1 x a_2_1 + a_1_1 x a_2_2 + a_1_2 x a_2_1 - a_1_2 x a_2_2 =\n');
fprintf(['         = %d x %d + %d x %d + %d x %d - %d x %d =\n' ...
         '         = %d + %d + %d - %d <= L^(C) = %d\n'], ...
         a_1_1, a_2_1, a_1_1, a_2_2, a_1_2, a_2_1, a_1_2, a_2_2, ...
         a_1_1 * a_2_1, a_1_1 * a_2_2, a_1_2 * a_2_1, a_1_2 * a_2_2, ...
         f_2_1_chsh_inequality_local_upper_bound_L);

% Print a blank line
fprintf('\n');

% Print of the variables/unknowns' solution
% of the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality (F_2_1)
fprintf('Such that:\n');
fprintf('  a_1_1 = %d\n', a_1_1);
fprintf('  a_1_2 = %d\n', a_1_2);
fprintf('  a_2_1 = %d\n', a_2_1);
fprintf('  a_2_2 = %d\n', a_2_2);


% Print a blank line
fprintf('\n');


num_final_total_observers = 4;


a = compute_f_n_1_chsh_inequality_local_upper_bound_L(num_final_total_observers);

disp(a);




% F_n = 1/2 (a_n + a'_n)xF_(n-1) + 1/2 (a_n - a'_n)F'_(n-1)


function f_2_1_chsh_inequality_local_upper_bound_L = ...
        compute_f_2_1_chsh_inequality_local_upper_bound_L()
        
    % Definition of the maximization optimization
    % problem for the estimation of the upper bound
    % for the CHSH Inequality (F_2_1)
    f_2_1_chsh_inequality_local_bound_optimization_prob = ...
        optimproblem('ObjectiveSense', 'max');
    
    % Definition of the optimization variables/unknowns
    % to be optimized in the maximization optimization
    % problem for the estimation of the upper bound
    % for the CHSH Inequality
    a_1_1 = optimvar('a_1_1', 1, 1, 'Type', 'integer');
    a_1_2 = optimvar('a_1_2', 1, 1, 'Type', 'integer');
    a_2_1 = optimvar('a_2_1', 1, 1, 'Type', 'integer');
    a_2_2 = optimvar('a_2_2', 1, 1, 'Type', 'integer');
    
    % Definition of the objective function of
    % the maximization optimization problem
    % for the estimation of the upper bound
    % for the CHSH Inequality (F_2_1)
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Objective = a_1_1(1) * a_2_1(1) + a_1_1(1) * a_2_2(1) + ... 
                     a_1_2(1) * a_2_1(1) - a_1_2(1) * a_2_2(1);
    
    
    % Definition of the constraints
    % for the a_1_1 value E [-1,+1]
    constraint_1 = a_1_1(1) >= -1;
    constraint_2 = a_1_1(1) <= 1;
    constraint_3 = a_1_1(1)^2 <= 1;
    constraint_4 = a_1_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_1_2 value E [-1,+1]
    constraint_5 = a_1_2(1) >= -1;
    constraint_6 = a_1_2(1) <= 1;
    constraint_7 = a_1_2(1)^2 <= 1;
    constraint_8 = a_1_2(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_2_1 value E [-1,+1]
    constraint_9 = a_2_1(1) >= -1;
    constraint_10 = a_2_1(1) <= 1;
    constraint_11 = a_2_1(1)^2 <= 1;
    constraint_12 = a_2_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_2_2 value E [-1,+1]
    constraint_13 = a_2_2(1) >= -1;
    constraint_14 = a_2_2(1) <= 1;
    constraint_15 = a_2_2(1)^2 <= 1;
    constraint_16 = a_2_2(1)^2 >= 1;
    
    
    % Setup of all the constraints
    % for the estimation values
    % a_1_1, a_1_2, a_2_1, and a_2_2,
    % for the maximization optimization
    % problem for the estimation of
    % the upper bound for the CHSH Inequality (F_2_1)
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons1 = constraint_1;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons2 = constraint_2;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons3 = constraint_3;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons4 = constraint_4;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons5 = constraint_5;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons6 = constraint_6;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons7 = constraint_7;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons8 = constraint_8;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons9 = constraint_9;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons10 = constraint_10;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons11 = constraint_11;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons12 = constraint_12;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons13 = constraint_13;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons14 = constraint_14;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons15 = constraint_15;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons16 = constraint_16;


    % Solve the maximization optimization problem
    % for the estimation of the (classical) local bound
    % for the CHSH Inequality (F_2_1)
    f_2_1_chsh_inequality_local_bound_optimization_prob_sol = ...
        solve(f_2_1_chsh_inequality_local_bound_optimization_prob);
    
    % Compute the a_1_1, a_1_2, a_2_1, and a_2_2 values
    % for the CHSH Inequality (F_2_1) that allow to
    % estimate the (classical) local upper bound L
    f_2_1_chsh_inequality_local_bound_a_1_1_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_1_1;
    f_2_1_chsh_inequality_local_bound_a_1_2_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_1_2;
    f_2_1_chsh_inequality_local_bound_a_2_1_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_2_1;
    f_2_1_chsh_inequality_local_bound_a_2_2_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_2_2;
    
    % Compute the individual expectation values
    % for the CHSH Inequality (F_2_1) that allow to estimate
    % the (classical) local upper bound L
    a_1_1 = f_2_1_chsh_inequality_local_bound_a_1_1_value(1);
    a_1_2 = f_2_1_chsh_inequality_local_bound_a_1_2_value(1);
    a_2_1 = f_2_1_chsh_inequality_local_bound_a_2_1_value(1);
    a_2_2 = f_2_1_chsh_inequality_local_bound_a_2_2_value(1);
    
    % Compute the (classical) local
    % upper bound L for the CHSH Inequality (F_2_1)
    f_2_1_chsh_inequality_local_upper_bound_L = ...
        a_1_1 * a_2_1 + a_1_1 * a_2_2 + ...
        a_1_2 * a_2_1 - a_1_2 * a_2_2;

end


function f_2_2_chsh_inequality_local_upper_bound_L = ...
        compute_f_2_2_chsh_inequality_local_upper_bound_L()
        
    % Definition of the maximization optimization
    % problem for the estimation of the upper bound
    % for the CHSH Inequality (F_2_2)
    f_2_2_chsh_inequality_local_bound_optimization_prob = ...
        optimproblem('ObjectiveSense', 'max');
    
    % Definition of the optimization variables/unknowns
    % to be optimized in the maximization optimization
    % problem for the estimation of the upper bound
    % for the CHSH Inequality (F_2_2)
    a_1_1 = optimvar('a_1_1', 1, 1, 'Type', 'integer');
    a_1_2 = optimvar('a_1_2', 1, 1, 'Type', 'integer');
    a_2_1 = optimvar('a_2_1', 1, 1, 'Type', 'integer');
    a_2_2 = optimvar('a_2_2', 1, 1, 'Type', 'integer');
    
    % Definition of the objective function of
    % the maximization optimization problem
    % for the estimation of the upper bound
    % for the CHSH Inequality (F_2_2)
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Objective = a_1_2(1) * a_2_2(1) + a_1_2(1) * a_2_1(1) + ... 
                     a_1_1(1) * a_2_2(1) - a_1_1(1) * a_2_1(1);
    
    % Definition of the constraints
    % for the a_1_1 value E [-1,+1]
    constraint_1 = a_1_1(1) >= -1;
    constraint_2 = a_1_1(1) <= 1;
    constraint_3 = a_1_1(1)^2 <= 1;
    constraint_4 = a_1_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_1_2 value E [-1,+1]
    constraint_5 = a_1_2(1) >= -1;
    constraint_6 = a_1_2(1) <= 1;
    constraint_7 = a_1_2(1)^2 <= 1;
    constraint_8 = a_1_2(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_2_1 value E [-1,+1]
    constraint_9 = a_2_1(1) >= -1;
    constraint_10 = a_2_1(1) <= 1;
    constraint_11 = a_2_1(1)^2 <= 1;
    constraint_12 = a_2_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_2_2 value E [-1,+1]
    constraint_13 = a_2_2(1) >= -1;
    constraint_14 = a_2_2(1) <= 1;
    constraint_15 = a_2_2(1)^2 <= 1;
    constraint_16 = a_2_2(1)^2 >= 1;
    
    
    % Setup of all the constraints
    % for the estimation values
    % a_1_1, a_1_2, a_2_1, and a_2_2,
    % for the maximization optimization
    % problem for the estimation of
    % the upper bound for the CHSH Inequality (F_2_2)
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons1 = constraint_1;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons2 = constraint_2;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons3 = constraint_3;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons4 = constraint_4;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons5 = constraint_5;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons6 = constraint_6;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons7 = constraint_7;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons8 = constraint_8;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons9 = constraint_9;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons10 = constraint_10;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons11 = constraint_11;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons12 = constraint_12;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons13 = constraint_13;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons14 = constraint_14;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons15 = constraint_15;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons16 = constraint_16;


    % Solve the maximization optimization problem
    % for the estimation of the (classical) local bound
    % for the CHSH Inequality (F_2_2)
    f_2_2_chsh_inequality_local_bound_optimization_prob_sol = ...
        solve(f_2_2_chsh_inequality_local_bound_optimization_prob);
    
    % Compute the a_1_1, a_1_2, a_2_1, and a_2_2 values
    % for the CHSH Inequality (F_2_2) that allow to
    % estimate the (classical) local upper bound L
    f_2_2_chsh_inequality_local_bound_a_1_1_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_1_1;
    f_2_2_chsh_inequality_local_bound_a_1_2_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_1_2;
    f_2_2_chsh_inequality_local_bound_a_2_1_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_2_1;
    f_2_2_chsh_inequality_local_bound_a_2_2_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_2_2;
    
    % Compute the individual expectation values
    % for the CHSH Inequality (F_2_2) that allow to estimate
    % the (classical) local upper bound L
    a_1_1 = f_2_2_chsh_inequality_local_bound_a_1_1_value(1);
    a_1_2 = f_2_2_chsh_inequality_local_bound_a_1_2_value(1);
    a_2_1 = f_2_2_chsh_inequality_local_bound_a_2_1_value(1);
    a_2_2 = f_2_2_chsh_inequality_local_bound_a_2_2_value(1);


    f_2_2_chsh_inequality_local_upper_bound_L = ...
        a_1_2 * a_2_2 + a_1_2 * a_2_1 + ...
        a_1_1 * a_2_2 - a_1_1 * a_2_1;
    
    
end


function f_n_1_chsh_inequality_local_upper_bound_L = ...
        compute_f_n_1_chsh_inequality_local_upper_bound_L(n)
    
    if n == 2
        
        f_n_1_chsh_inequality_local_upper_bound_L = ...
            compute_f_2_1_chsh_inequality_local_upper_bound_L();

    end


    if n > 2
       
        % Definition of the maximization optimization
        % problem for the estimation of the upper bound
        % for the CHSH Inequality (F_n_1)
        f_n_1_chsh_inequality_local_bound_optimization_prob = ...
            optimproblem('ObjectiveSense', 'max');
        
        % Definition of the optimization variables/unknowns
        % to be optimized in the maximization optimization
        % problem for the estimation of the upper bound
        % for the CHSH Inequality
        a_n_1 = optimvar('a_n_1', 1, 1, 'Type', 'integer');
        a_n_2 = optimvar('a_n_2', 1, 1, 'Type', 'integer');
        
        % Definition of the objective function of
        % the maximization optimization problem
        % for the estimation of the upper bound
        % for the CHSH Inequality (F_n_1)
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Objective = (1 / 2) * (a_n_1 + a_n_2) * ...
                            compute_f_n_1_chsh_inequality_local_upper_bound_L(n - 1) + ...
                         (1 / 2) * (a_n_1 - a_n_2) * ...
                            compute_f_n_2_chsh_inequality_local_upper_bound_L(n - 1);
        
        % Definition of the constraints
        % for the a_n_1 value E [-1,+1]
        constraint_1 = a_n_1(1) >= -1;
        constraint_2 = a_n_1(1) <= 1;
        constraint_3 = a_n_1(1)^2 <= 1;
        constraint_4 = a_n_1(1)^2 >= 1;
        
        % Definition of the constraints
        % for the a_n_2 value E [-1,+1]
        constraint_5 = a_n_2(1) >= -1;
        constraint_6 = a_n_2(1) <= 1;
        constraint_7 = a_n_2(1)^2 <= 1;
        constraint_8 = a_n_2(1)^2 >= 1;
        
        % Setup of all the constraints
        % for the estimation values a_n_1, and a_n_2,
        % for the maximization optimization
        % problem for the estimation of
        % the upper bound for the CHSH Inequality (F_2_1)
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons1 = constraint_1;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons2 = constraint_2;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons3 = constraint_3;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons4 = constraint_4;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons5 = constraint_5;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons6 = constraint_6;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons7 = constraint_7;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons8 = constraint_8;
        
        
        % Solve the maximization optimization problem
        % for the estimation of the (classical) local bound
        % for the CHSH Inequality (F_n_1)
        f_n_1_chsh_inequality_local_bound_optimization_prob_sol = ...
            solve(f_n_1_chsh_inequality_local_bound_optimization_prob);
        
        % Compute the a_n_1, and a_n_2 values
        % for the CHSH Inequality (F_n_1) that allow to
        % estimate the (classical) local upper bound L
        f_n_1_chsh_inequality_local_bound_a_n_1_value = ...
            f_n_1_chsh_inequality_local_bound_optimization_prob_sol.a_n_1;
        f_n_1_chsh_inequality_local_bound_a_n_2_value = ...
            f_n_1_chsh_inequality_local_bound_optimization_prob_sol.a_n_2;
        
        % Compute the individual expectation values
        % for the CHSH Inequality (F_n_1) that allow to estimate
        % the (classical) local upper bound L
        a_n_1 = f_n_1_chsh_inequality_local_bound_a_n_1_value(1);
        a_n_2 = f_n_1_chsh_inequality_local_bound_a_n_2_value(1);
        
        % Compute the (classical) local
        % upper bound L for the CHSH Inequality (F_n_1)
        f_n_1_chsh_inequality_local_upper_bound_L = ...
            (1 / 2) * (a_n_1 + a_n_2) * ...
                compute_f_n_1_chsh_inequality_local_upper_bound_L(n - 1) + ...
            (1 / 2) * (a_n_1 - a_n_2) * ...
                compute_f_n_2_chsh_inequality_local_upper_bound_L(n - 1);
    
    end

end


function f_n_2_chsh_inequality_local_upper_bound_L = ...
        compute_f_n_2_chsh_inequality_local_upper_bound_L(n)
    
    if n == 2
        
        f_n_2_chsh_inequality_local_upper_bound_L = ...
            compute_f_2_2_chsh_inequality_local_upper_bound_L();
        
    end


    if n > 2
       
        % Definition of the maximization optimization
        % problem for the estimation of the upper bound
        % for the CHSH Inequality (F_n_2)
        f_n_2_chsh_inequality_local_bound_optimization_prob = ...
            optimproblem('ObjectiveSense', 'max');
        
        % Definition of the optimization variables/unknowns
        % to be optimized in the maximization optimization
        % problem for the estimation of the upper bound
        % for the CHSH Inequality
        a_n_1 = optimvar('a_n_1', 1, 1, 'Type', 'integer');
        a_n_2 = optimvar('a_n_2', 1, 1, 'Type', 'integer');
        
        % Definition of the objective function of
        % the maximization optimization problem
        % for the estimation of the upper bound
        % for the CHSH Inequality (F_n_1)
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Objective = (1 / 2) * (a_n_2 + a_n_1) * ...
                            compute_f_n_1_chsh_inequality_local_upper_bound_L(n - 1) + ...
                         (1 / 2) * (a_n_2 - a_n_1) * ...
                            compute_f_n_2_chsh_inequality_local_upper_bound_L(n - 1);
        
        % Definition of the constraints
        % for the a_n_1 value E [-1,+1]
        constraint_1 = a_n_1(1) >= -1;
        constraint_2 = a_n_1(1) <= 1;
        constraint_3 = a_n_1(1)^2 <= 1;
        constraint_4 = a_n_1(1)^2 >= 1;
        
        % Definition of the constraints
        % for the a_n_2 value E [-1,+1]
        constraint_5 = a_n_2(1) >= -1;
        constraint_6 = a_n_2(1) <= 1;
        constraint_7 = a_n_2(1)^2 <= 1;
        constraint_8 = a_n_2(1)^2 >= 1;
        
        % Setup of all the constraints
        % for the estimation values a_n_1, and a_n_2,
        % for the maximization optimization
        % problem for the estimation of
        % the upper bound for the CHSH Inequality (F_n_2)
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons1 = constraint_1;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons2 = constraint_2;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons3 = constraint_3;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons4 = constraint_4;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons5 = constraint_5;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons6 = constraint_6;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons7 = constraint_7;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons8 = constraint_8;
        
        
        % Solve the maximization optimization problem
        % for the estimation of the (classical) local bound
        % for the CHSH Inequality (F_n_2)
        f_n_2_chsh_inequality_local_bound_optimization_prob_sol = ...
            solve(f_n_2_chsh_inequality_local_bound_optimization_prob);
        
        % Compute the a_n_1, and a_n_2 values
        % for the CHSH Inequality (F_n_2) that allow to
        % estimate the (classical) local upper bound L
        f_n_2_chsh_inequality_local_bound_a_n_1_value = ...
            f_n_2_chsh_inequality_local_bound_optimization_prob_sol.a_n_1;
        f_n_2_chsh_inequality_local_bound_a_n_2_value = ...
            f_n_2_chsh_inequality_local_bound_optimization_prob_sol.a_n_2;
        
        % Compute the individual expectation values
        % for the CHSH Inequality (F_n_2) that allow to estimate
        % the (classical) local upper bound L
        a_n_1 = f_n_2_chsh_inequality_local_bound_a_n_1_value(1);
        a_n_2 = f_n_2_chsh_inequality_local_bound_a_n_2_value(1);
        
        % Compute the (classical) local
        % upper bound L for the CHSH Inequality (F_n_2)
        f_n_2_chsh_inequality_local_upper_bound_L = ...
            (1 / 2) * (a_n_2 - a_n_1) * ...
                compute_f_n_1_chsh_inequality_local_upper_bound_L(n - 1) + ...
            (1 / 2) * (a_n_1 + a_n_2) * ...
                compute_f_n_2_chsh_inequality_local_upper_bound_L(n - 1);

         
    end

end