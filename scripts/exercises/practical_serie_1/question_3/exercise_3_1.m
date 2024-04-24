% Clear Command Window
clc;


% Compute the CHSH Coefficients Matrix
% for the CHSH Inequality with
% the notation a_0 x b_0 + a_0 x b_0 + 
%              a_1 x b_0 - a_1 x b_1 <= L,
% which results on the following matrix:
%      B_0  B_1  B_2
% A_0   0    0    0
% A_1   0    1    1
% A_2   0    1   -1
chsh_coefficients_matrix = ...
    [0  0  0; ...
     0  1  1; ...
     0  1 -1];

% Create the number of inputs and outputs
% for Alice and Bob, regarding the setup of
% the CHSH Inequality
chsh_num_inputs_alice = 2;
chsh_num_inputs_bob = 2;
chsh_num_outputs_alice = 2;
chsh_num_outputs_bob = 2;

% Compute the vector of 4 elements
% with the number of outputs of Alice,
% number of outputs of Bob, number of inputs of Alice,
% and number of inputs of Bob, in this order
chsh_description = [ chsh_num_outputs_alice chsh_num_outputs_bob ...
                     chsh_num_inputs_alice chsh_num_inputs_bob ];

% Compute the maximum algebraic
% upper bound for the CHSH Inequality
chsh_inequality_algebraic_max_upper_bound_L = ...
    BellInequalityMax( chsh_coefficients_matrix, ...
                       chsh_description, ...
                       'fc', 'nosignal');

% Print of the estimation of the maximum algebraic
% upper bound for the CHSH Inequality
fprintf(['Al Upper Bound L ' ...
         'for CHSH Inequality:\n'])
fprintf('  <A1B1> + <A1B2> + <A2B1> - <A2B2> <= L^(A) = %d\n', ...
             chsh_inequality_algebraic_max_upper_bound_L);

% Print a blank line
fprintf('\n');


% Print of the constraints of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality
fprintf('Such that:\n');
fprintf('  <A1B1> E {-1, +1}\n');
fprintf('  <A1B2> E {-1, +1}\n');
fprintf('  <A2B1> E {-1, +1}\n');
fprintf('  <A2B2> E {-1, +1}\n');


% Print a blank line
fprintf('\n');
