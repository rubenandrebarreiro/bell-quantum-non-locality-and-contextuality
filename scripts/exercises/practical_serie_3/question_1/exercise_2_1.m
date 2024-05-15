% Clear Command Window
clc;


% Definition of the maximization optimization
% problem for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
f_2_1_mermin_inequality_algebraic_max_optimization_prob = ...
    optimproblem('ObjectiveSense', 'max');

% Definition of the optimization variables/unknowns
% to be optimized in the maximization optimization
% problem for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
a_xy = optimvar('a_xy', 4, 1, 'Type', 'integer');

% Definition of the objective function of
% the maximization optimization problem
% for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Objective = a_xy(1) + a_xy(2) + a_xy(3) - a_xy(4);

% Definition of the constraints
% for the estimation value a_{11} value E {-1,+1}
constraint_1 = a_xy(1) >= -1;
constraint_2 = a_xy(1) <= 1;
constraint_3 = a_xy(1)^2 <= 1;
constraint_4 = a_xy(1)^2 >= 1;

% Definition of the constraints
% for the estimation value a_{12} value E {-1,+1}
constraint_5 = a_xy(2) >= -1;
constraint_6 = a_xy(2) <= 1;
constraint_7 = a_xy(2)^2 <= 1;
constraint_8 = a_xy(2)^2 >= 1;

% Definition of the constraints
% for the estimation value a_{21} value E {-1,+1}
constraint_9 = a_xy(3) >= -1;
constraint_10 = a_xy(3) <= 1;
constraint_11 = a_xy(3)^2 <= 1;
constraint_12 = a_xy(3)^2 >= 1;

% Definition of the constraints
% for the estimation value a_{22} value E {-1,+1}
constraint_13 = a_xy(4) >= -1;
constraint_14 = a_xy(4) <= 1;
constraint_15 = a_xy(4)^2 <= 1;
constraint_16 = a_xy(4)^2 >= 1;


% Setup of all the constraints
% for the estimation values
% a_{11}, a_{12}, a_{21}, and a_{22},
% for the maximization optimization
% problem for the estimation of
% the upper bound for the Mermin Inequality (F_2_1)
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons1 = constraint_1;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons2 = constraint_2;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons3 = constraint_3;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons4 = constraint_4;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons5 = constraint_5;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons6 = constraint_6;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons7 = constraint_7;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons8 = constraint_8;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons9 = constraint_9;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons10 = constraint_10;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons11 = constraint_11;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons12 = constraint_12;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons13 = constraint_13;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons14 = constraint_14;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons15 = constraint_15;
f_2_1_mermin_inequality_algebraic_max_optimization_prob...
    .Constraints.cons16 = constraint_16;

% Show the mathematical formulation of
% the maximization optimization problem
% for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
show(f_2_1_mermin_inequality_algebraic_max_optimization_prob);

% Solve the maximization optimization problem
% for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
f_2_1_mermin_inequality_algebraic_max_optimization_prob_sol = ...
    solve(f_2_1_mermin_inequality_algebraic_max_optimization_prob);

% Compute all the expectation values
% for the Mermin Inequality (F_2_1) that allow to
% estimate the maximum algebraic upper bound L
f_2_1_mermin_inequality_algebraic_max_exp_values = ...
    f_2_1_mermin_inequality_algebraic_max_optimization_prob_sol.a_xy;

% Compute the individual expectation values
% for the Mermin Inequality (F_2_1) that allow to
% estimate the maximum algebraic upper bound L
a_11 = f_2_1_mermin_inequality_algebraic_max_exp_values(1);
a_12 = f_2_1_mermin_inequality_algebraic_max_exp_values(2);
a_21 = f_2_1_mermin_inequality_algebraic_max_exp_values(3);
a_22 = f_2_1_mermin_inequality_algebraic_max_exp_values(4);

% Compute the maximum algebraic
% upper bound L for the Mermin Inequality (F_2_1)
f_2_1_mermin_inequality_algebraic_max_upper_bound_L = ...
    a_11 + a_12 + a_21 - a_22;

% Print a blank line
fprintf('\n');

% Print of the mathematical form
% for the Mermin Inequality (F_2_1),
% in the correlation form
fprintf('Mermin (F_2_1) Inequality:\n');
fprintf('  a_{11} + a_{12} + a_{12} - a_{22} <= L\n');

% Print a blank line
fprintf('\n');

% Print of the constraints of
% the maximization optimization problem
% for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
fprintf('Such that:\n');
fprintf('  a_{11} E {-1, +1}\n');
fprintf('  a_{12} E {-1, +1}\n');
fprintf('  a_{21} E {-1, +1}\n');
fprintf('  a_{22} E {-1, +1}\n');


% Print a blank line
fprintf('\n');

% Print a blank line
fprintf('\n');

% Print of the complete solution of
% the maximization optimization problem
% for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
fprintf(['Maximum Algebraic Upper Bound L ' ...
         'for Mermin Inequality (F_2_1):\n']);
fprintf('  a_11 + a_12 + a_21 - a_22 =\n');
fprintf('         = %d + %d + %d - %d <= L^(C) = %d\n', ...
        a_11, a_12, a_21, a_22, ...
        f_2_1_mermin_inequality_algebraic_max_upper_bound_L);

% Print a blank line
fprintf('\n');

% Print of the variables/unknowns' solution
% of the maximization optimization problem
% for the estimation of the upper bound
% for the Mermin Inequality (F_2_1)
fprintf('Such that:\n');
fprintf('  a_{11} = %d\n', a_11);
fprintf('  a_{12} = %d\n', a_12);
fprintf('  a_{21} = %d\n', a_21);
fprintf('  a_{22} = %d\n', a_22);


% Print a blank line
fprintf('\n');



num_final_total_observers = 2;

a = compute_f_n_2_mermin_inequality_algebraic_max_upper_bound_L...
    (num_final_total_observers);

disp(a);


function f_2_1_mermin_inequality_algebraic_max_upper_bound_L = ...
         compute_f_2_1_mermin_inequality_algebraic_max_upper_bound_L()
        
    % Definition of the maximization optimization
    % problem for the estimation of the upper bound
    % for the Mermin Inequality (F_2_1)
    f_2_1_mermin_inequality_algebraic_max_optimization_prob = ...
        optimproblem('ObjectiveSense', 'max');

    % Definition of the optimization variables/unknowns
    % to be optimized in the maximization optimization
    % problem for the estimation of the upper bound
    % for the Mermin Inequality (F_2_1)
    a_xy = optimvar('a_xy', 4, 1, 'Type', 'integer');

    % Definition of the objective function of
    % the maximization optimization problem
    % for the estimation of the upper bound
    % for the Mermin Inequality (F_2_1)
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Objective = a_xy(1) + a_xy(2) + a_xy(3) - a_xy(4);


    % Definition of the constraints
    % for the estimation value a_{11} value E {-1,+1}
    constraint_1 = a_xy(1) >= -1;
    constraint_2 = a_xy(1) <= 1;
    constraint_3 = a_xy(1)^2 <= 1;
    constraint_4 = a_xy(1)^2 >= 1;
    
    % Definition of the constraints
    % for the estimation value a_{12} value E {-1,+1}
    constraint_5 = a_xy(2) >= -1;
    constraint_6 = a_xy(2) <= 1;
    constraint_7 = a_xy(2)^2 <= 1;
    constraint_8 = a_xy(2)^2 >= 1;
    
    % Definition of the constraints
    % for the estimation value a_{21} value E {-1,+1}
    constraint_9 = a_xy(3) >= -1;
    constraint_10 = a_xy(3) <= 1;
    constraint_11 = a_xy(3)^2 <= 1;
    constraint_12 = a_xy(3)^2 >= 1;
    
    % Definition of the constraints
    % for the estimation value a_{22} value E {-1,+1}
    constraint_13 = a_xy(4) >= -1;
    constraint_14 = a_xy(4) <= 1;
    constraint_15 = a_xy(4)^2 <= 1;
    constraint_16 = a_xy(4)^2 >= 1;
    
    
    % Setup of all the constraints
    % for the estimation values
    % a_{11}, a_{12}, a_{21}, and a_{22},
    % for the maximization optimization
    % problem for the estimation of
    % the upper bound for the Mermin Inequality (F_2_1)
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons1 = constraint_1;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons2 = constraint_2;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons3 = constraint_3;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons4 = constraint_4;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons5 = constraint_5;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons6 = constraint_6;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons7 = constraint_7;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons8 = constraint_8;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons9 = constraint_9;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons10 = constraint_10;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons11 = constraint_11;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons12 = constraint_12;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons13 = constraint_13;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons14 = constraint_14;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons15 = constraint_15;
    f_2_1_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons16 = constraint_16;
    
    % Solve the maximization optimization problem
    % for the estimation of the upper bound
    % for the Mermin Inequality (F_2_1)
    f_2_1_mermin_inequality_algebraic_max_optimization_prob_sol = ...
        solve(f_2_1_mermin_inequality_algebraic_max_optimization_prob);

    % Compute all the expectation values
    % for the Mermin Inequality (F_2_1) that allow to
    % estimate the maximum algebraic upper bound L
    f_2_1_mermin_inequality_algebraic_max_exp_values = ...
        f_2_1_mermin_inequality_algebraic_max_optimization_prob_sol.a_xy;

    % Compute the individual expectation values
    % for the Mermin Inequality (F_2_1) that allow to
    % estimate the maximum algebraic upper bound L
    a_11 = f_2_1_mermin_inequality_algebraic_max_exp_values(1);
    a_12 = f_2_1_mermin_inequality_algebraic_max_exp_values(2);
    a_21 = f_2_1_mermin_inequality_algebraic_max_exp_values(3);
    a_22 = f_2_1_mermin_inequality_algebraic_max_exp_values(4);

    % Compute the maximum algebraic
    % upper bound L for the Mermin Inequality (F_2_1)
    f_2_1_mermin_inequality_algebraic_max_upper_bound_L = ...
        a_11 + a_12 + a_21 - a_22;

end


function f_2_2_mermin_inequality_algebraic_max_upper_bound_L = ...
         compute_f_2_2_mermin_inequality_algebraic_max_upper_bound_L()
        
    % Definition of the maximization optimization
    % problem for the estimation of the upper bound
    % for the Mermin Inequality (F_2_2)
    f_2_2_mermin_inequality_algebraic_max_optimization_prob = ...
        optimproblem('ObjectiveSense', 'max');

    % Definition of the optimization variables/unknowns
    % to be optimized in the maximization optimization
    % problem for the estimation of the upper bound
    % for the Mermin Inequality (F_2_2)
    a_xy = optimvar('a_xy', 4, 1, 'Type', 'integer');

    % Definition of the objective function of
    % the maximization optimization problem
    % for the estimation of the upper bound
    % for the Mermin Inequality (F_2_2)
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Objective = a_xy(4) + a_xy(3) + a_xy(2) - a_xy(1);


    % Definition of the constraints
    % for the estimation value a_{11} value E {-1,+1}
    constraint_1 = a_xy(1) >= -1;
    constraint_2 = a_xy(1) <= 1;
    constraint_3 = a_xy(1)^2 <= 1;
    constraint_4 = a_xy(1)^2 >= 1;
    
    % Definition of the constraints
    % for the estimation value a_{12} value E {-1,+1}
    constraint_5 = a_xy(2) >= -1;
    constraint_6 = a_xy(2) <= 1;
    constraint_7 = a_xy(2)^2 <= 1;
    constraint_8 = a_xy(2)^2 >= 1;
    
    % Definition of the constraints
    % for the estimation value a_{21} value E {-1,+1}
    constraint_9 = a_xy(3) >= -1;
    constraint_10 = a_xy(3) <= 1;
    constraint_11 = a_xy(3)^2 <= 1;
    constraint_12 = a_xy(3)^2 >= 1;
    
    % Definition of the constraints
    % for the estimation value a_{22} value E {-1,+1}
    constraint_13 = a_xy(4) >= -1;
    constraint_14 = a_xy(4) <= 1;
    constraint_15 = a_xy(4)^2 <= 1;
    constraint_16 = a_xy(4)^2 >= 1;
    
    
    % Setup of all the constraints
    % for the estimation values
    % a_{11}, a_{12}, a_{21}, and a_{22},
    % for the maximization optimization
    % problem for the estimation of
    % the upper bound for the Mermin Inequality (F_2_2)
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons1 = constraint_1;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons2 = constraint_2;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons3 = constraint_3;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons4 = constraint_4;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons5 = constraint_5;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons6 = constraint_6;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons7 = constraint_7;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons8 = constraint_8;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons9 = constraint_9;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons10 = constraint_10;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons11 = constraint_11;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons12 = constraint_12;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons13 = constraint_13;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons14 = constraint_14;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons15 = constraint_15;
    f_2_2_mermin_inequality_algebraic_max_optimization_prob...
        .Constraints.cons16 = constraint_16;
    
    % Solve the maximization optimization problem
    % for the estimation of the upper bound
    % for the Mermin Inequality (F_2_2)
    f_2_2_mermin_inequality_algebraic_max_optimization_prob_sol = ...
        solve(f_2_2_mermin_inequality_algebraic_max_optimization_prob);

    % Compute all the expectation values
    % for the Mermin Inequality (F_2_2) that allow to
    % estimate the maximum algebraic upper bound L
    f_2_2_mermin_inequality_algebraic_max_exp_values = ...
        f_2_2_mermin_inequality_algebraic_max_optimization_prob_sol.a_xy;

    % Compute the individual expectation values
    % for the Mermin Inequality (F_2_2) that allow to
    % estimate the maximum algebraic upper bound L
    a_11 = f_2_2_mermin_inequality_algebraic_max_exp_values(1);
    a_12 = f_2_2_mermin_inequality_algebraic_max_exp_values(2);
    a_21 = f_2_2_mermin_inequality_algebraic_max_exp_values(3);
    a_22 = f_2_2_mermin_inequality_algebraic_max_exp_values(4);

    % Compute the maximum algebraic
    % upper bound L for the Mermin Inequality (F_2_2)
    f_2_2_mermin_inequality_algebraic_max_upper_bound_L = ...
        a_11 + a_12 + a_21 - a_22;

end

function f_n_1_mermin_inequality_algebraic_max_upper_bound_L = ...
         compute_f_n_1_mermin_inequality_algebraic_max_upper_bound_L(n)
    
    if n == 2
        
        f_n_1_mermin_inequality_algebraic_max_upper_bound_L = ...
            compute_f_2_1_mermin_inequality_algebraic_max_upper_bound_L();

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


function f_n_2_mermin_inequality_algebraic_max_upper_bound_L = ...
         compute_f_n_2_mermin_inequality_algebraic_max_upper_bound_L(n)
    
    if n == 2
        
        f_n_2_mermin_inequality_algebraic_max_upper_bound_L = ...
            compute_f_2_2_mermin_inequality_algebraic_max_upper_bound_L();
        
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