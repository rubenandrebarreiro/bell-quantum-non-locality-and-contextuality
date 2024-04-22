% Clear Command Window
clc;


% Definition of the maximization optimization
% problem for the estimation of the upper bound
% for the CHSH Inequality
chsh_inequality_local_bound_optimization_prob = ...
    optimproblem('ObjectiveSense', 'max');

% Definition of the optimization variables/unknowns
% to be optimized in the maximization optimization
% problem for the estimation of the upper bound
% for the CHSH Inequality
a_0 = optimvar('a_0', 1, 1, 'Type', 'integer');
a_1 = optimvar('a_1', 1, 1, 'Type', 'integer');
b_0 = optimvar('b_0', 1, 1, 'Type', 'integer');
b_1 = optimvar('b_1', 1, 1, 'Type', 'integer');

% Definition of the objective function of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality
chsh_inequality_local_bound_optimization_prob...
    .Objective = a_0(1) * b_0(1) + a_0 * b_1 + a_1*b_0 - a_1 * b_1;

% Definition of the constraints
% for the a_0 value E [-1,+1]
constraint_1 = a_0(1) >= -1;
constraint_2 = a_0(1) <= 1;
constraint_3 = a_0(1)^2 <= 1;
constraint_4 = a_0(1)^2 >= 1;

% Definition of the constraints
% for the a_1 value E [-1,+1]
constraint_5 = a_1(1) >= -1;
constraint_6 = a_1(1) <= 1;
constraint_7 = a_1(1)^2 <= 1;
constraint_8 = a_1(1)^2 >= 1;

% Definition of the constraints
% for the b_0 value E [-1,+1]
constraint_9 = b_0(1) >= -1;
constraint_10 = b_0(1) <= 1;
constraint_11 = b_0(1)^2 <= 1;
constraint_12 = b_0(1)^2 >= 1;

% Definition of the constraints
% for the b_1 value E [-1,+1]
constraint_13 = b_1(1) >= -1;
constraint_14 = b_1(1) <= 1;
constraint_15 = b_1(1)^2 <= 1;
constraint_16 = b_1(1)^2 >= 1;


% Setup of all the constraints
% for the estimation values
% a_0, a_1, b_0, and b_1,
% for the maximization optimization
% problem for the estimation of
% the upper bound for the CHSH Inequality
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons1 = constraint_1;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons2 = constraint_2;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons3 = constraint_3;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons4 = constraint_4;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons5 = constraint_5;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons6 = constraint_6;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons7 = constraint_7;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons8 = constraint_8;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons9 = constraint_9;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons10 = constraint_10;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons11 = constraint_11;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons12 = constraint_12;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons13 = constraint_13;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons14 = constraint_14;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons15 = constraint_15;
chsh_inequality_local_bound_optimization_prob...
    .Constraints.cons16 = constraint_16;

% Show the mathematical formulation of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality
show(chsh_inequality_local_bound_optimization_prob);

% Solve the maximization optimization problem
% for the estimation of the (classical) local bound
% for the CHSH Inequality 
chsh_inequality_local_bound_optimization_prob_sol = ...
    solve(chsh_inequality_local_bound_optimization_prob);

% Compute the a_0, a_1, b_0, and b_1 values
% for the CHSH Inequality that allow to
% estimate the (classical) local upper bound L
chsh_inequality_local_bound_a_0_value = ...
    chsh_inequality_local_bound_optimization_prob_sol.a_0;
chsh_inequality_local_bound_a_1_value = ...
    chsh_inequality_local_bound_optimization_prob_sol.a_1;
chsh_inequality_local_bound_b_0_value = ...
    chsh_inequality_local_bound_optimization_prob_sol.b_0;
chsh_inequality_local_bound_b_1_value = ...
    chsh_inequality_local_bound_optimization_prob_sol.b_1;

% Compute the individual expectation values
% for the CHSH Inequality that allow to estimate
% the (classical) local upper bound L
a_0 = chsh_inequality_local_bound_a_0_value(1);
a_1 = chsh_inequality_local_bound_a_1_value(1);
b_0 = chsh_inequality_local_bound_b_0_value(1);
b_1 = chsh_inequality_local_bound_b_1_value(1);

% Compute the (classical) local
% upper bound L for the CHSH Inequality
chsh_inequality_local_upper_bound_L = ...
    a_0 * b_0 + a_0 * b_0 + a_1 * b_0 - a_1 * b_1;

% Print a blank line
fprintf('\n');

% Print of the mathematical form
% for the CHSH Inequality,
% in the deterministic strategy form
fprintf('CHSH Inequality:\n');
fprintf('  a_0 x b_0 + a_0 x b_0 + a_1 x b_0 - a_1 x b_1 <= L\n');

% Print a blank line
fprintf('\n');

% Print of the constraints of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality
fprintf('Such that:\n');
fprintf('  a_0 E {-1, +1}\n');
fprintf('  a_1 E {-1, +1}\n');
fprintf('  b_0 E {-1, +1}\n');
fprintf('  b_1 E {-1, +1}\n');


% Print a blank line
fprintf('\n');

% Print a blank line
fprintf('\n');

% Print of the complete solution of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality
fprintf(['(Classical) Local Upper Bound L ' ...
         'for CHSH Inequality:\n'])
fprintf('  a_0 x b_0 + a_0 x b_0 + a_1 x b_0 - a_1 x b_1 =\n');
fprintf('         = %d + %d + %d - %d <= L^(C) = %d\n', ...
         a_0, a_1, b_0, b_1, ...
         chsh_inequality_local_upper_bound_L);

% Print a blank line
fprintf('\n');

% Print of the variables/unknowns' solution
% of the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality
fprintf('Such that:\n');
fprintf('  a_0 = %d\n', a_0);
fprintf('  a_1 = %d\n', a_1);
fprintf('  b_0 = %d\n', b_0);
fprintf('  b_1 = %d\n', b_1);


% Print a blank line
fprintf('\n');