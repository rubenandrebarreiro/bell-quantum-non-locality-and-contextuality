% Clear Command Window
clc;


% Create the ket vector for the two-dimensional
% quantum Einstein-Podolski-Rosen (EPR) pair denoted
% as |{phi}^{+}> = ( 1 / sqrt(2) ) * ([1, 0, 0, 1])^T
ket_epr_pair_phi_plus = ( 1 / sqrt(2) ) * [ 1, 0, 0, 1 ].';

% Create the bra vector for the two-dimensional
% quantum Einstein-Podolski-Rosen (EPR) pair denoted
% as <{phi}^{+}| = ( 1 / sqrt(2) ) * ([1, 0, 0, 1])
bra_epr_pair_phi_plus = ket_epr_pair_phi_plus.';

% Compute the quantum upper bound
% for the CHSH Inequality
chsh_inequality_quantum_upper_bound_L = 2 * sqrt(2);

% Compute the error tolerance for
% the estimation of the several quantum
% upper bounds for the CHSH Inequality
error_tol = 1e+10;


% Create the Pauli I (sigma_i) Matrix
pauli_sigma_i = full(Pauli('I'));

% Create the Pauli X (sigma_x) Matrix
pauli_sigma_x = full(Pauli('X'));

% Create the Pauli Y (sigma_y) Matrix
pauli_sigma_y = full(Pauli('Y'));

% Create the Pauli Z (sigma_z) Matrix
pauli_sigma_z = full(Pauli('Z'));

% Create vector of the Pauli X, Y, and Z
% (sigma_x, sigma_y, sigma_z) Matrices
pauli_sigmas_vec = { pauli_sigma_x, ...
                     pauli_sigma_y, ...
                     pauli_sigma_z };


% Compute the density
% matrix rho = |psi><psi| for
% the two-dimensional quantum state psi
rho_density_matrix = ket_epr_pair_phi_plus * bra_epr_pair_phi_plus;


% Create the measurement choices
% for Alice as (sigma_x, sigma_z)
alice_measurements = { pauli_sigma_x, ...
                       pauli_sigma_z };

% Create the measurement choices for Bob as
% diagonal and anti-diagonal measurements, i.e.,
% ( ( 1 / sqrt(2) ) x ( sigma_x + sigma_z),
%   ( 1 / sqrt(2) ) x ( sigma_x - sigma_z) )
bob_measurements = { ( 1 / sqrt(2) ) * ( pauli_sigma_x + pauli_sigma_z ), ...
                     ( 1 / sqrt(2) ) * ( pauli_sigma_x - pauli_sigma_z ) };

% Compute the eigenvectors matrix for
% the (sigma_x) measurement choice of Alice
[eigenvectors_matrix, ~] = ...
    eig(alice_measurements{1});

% Compute the ket vector for the 1st basis
% for the (sigma_x) measurement choice of Alice
alice_meas_0_ket_basis_0 = [ eigenvectors_matrix(1:2) ].';

% Compute the ket vector for the 2nd basis
% for the (sigma_x) measurement choice of Alice
alice_meas_0_ket_basis_1 = [ eigenvectors_matrix(3:4) ].';

% Compute the bra vector for the 1st basis
% for the (sigma_x) measurement choice of Alice
alice_meas_0_bra_basis_0 = alice_meas_0_ket_basis_0.';

% Compute the bra vector for the 2nd basis
% for the (sigma_x) measurement choice of Alice
alice_meas_0_bra_basis_1 = alice_meas_0_ket_basis_1.';


% Compute the density matrices for the bases of
% the (sigma_x) measurement choice of Alice
alice_meas_0 = { alice_meas_0_ket_basis_0 * alice_meas_0_bra_basis_0, ...
                 alice_meas_0_ket_basis_1 * alice_meas_0_bra_basis_1 };


% Compute the eigenvectors matrix for
% the (sigma_z) measurement choice of Alice
[eigenvectors_matrix, ~] = ...
    eig(alice_measurements{2});

% Compute the ket vector for the 1st basis
% for the (sigma_z) measurement choice of Alice
alice_meas_1_ket_basis_0 = [ eigenvectors_matrix(1:2) ].';

% Compute the ket vector for the 2nd basis
% for the (sigma_z) measurement choice of Alice
alice_meas_1_ket_basis_1 = [ eigenvectors_matrix(3:4) ].';

% Compute the bra vector for the 1st basis
% for the (sigma_z) measurement choice of Alice
alice_meas_1_bra_basis_0 = alice_meas_1_ket_basis_0.';

% Compute the bra vector for the 2nd basis
% for the (sigma_z) measurement choice of Alice
alice_meas_1_bra_basis_1 = alice_meas_1_ket_basis_1.';


% Compute the density matrices for the bases of
% the (sigma_z) measurement choice of Alice
alice_meas_1 = { alice_meas_1_ket_basis_0 * alice_meas_1_bra_basis_0, ...
                 alice_meas_1_ket_basis_1 * alice_meas_1_bra_basis_1 };


% Compute the eigenvectors matrix for
% the diagonal measurement choice of Bob
[eigenvectors_matrix, ~] = ...
    eig(bob_measurements{1});

% Compute the ket vector for the 1st basis
% for the diagonal measurement choice of Bob
bob_meas_0_ket_basis_0 = [ eigenvectors_matrix(1:2) ].';

% Compute the ket vector for the 2nd basis
% for the diagonal measurement choice of Bob
bob_meas_0_ket_basis_1 = [ eigenvectors_matrix(3:4) ].';

% Compute the bra vector for the 1st basis
% for the diagonal measurement choice of Bob
bob_meas_0_bra_basis_0 = bob_meas_0_ket_basis_0.';

% Compute the bra vector for the 2nd basis
% for the diagonal measurement choice of Bob
bob_meas_0_bra_basis_1 = bob_meas_0_ket_basis_1.';


% Compute the density matrices for the bases of
% the diagonal measurement choice of Bob
bob_meas_0 = { bob_meas_0_ket_basis_0 * bob_meas_0_bra_basis_0, ...
               bob_meas_0_ket_basis_1 * bob_meas_0_bra_basis_1 };


% Compute the eigenvectors matrix for
% the anti-diagonal measurement choice of Bob
[eigenvectors_matrix, ~] = ...
    eig(bob_measurements{2});

% Compute the ket vector for the 1st basis
% for the anti-diagonal measurement choice of Bob
bob_meas_1_ket_basis_0 = [ eigenvectors_matrix(1:2) ].';

% Compute the ket vector for the 2nd basis
% for the anti-diagonal measurement choice of Bob
bob_meas_1_ket_basis_1 = [ eigenvectors_matrix(3:4) ].';

% Compute the bra vector for the 1st basis
% for the anti-diagonal measurement choice of Bob
bob_meas_1_bra_basis_0 = bob_meas_1_ket_basis_0.';

% Compute the bra vector for the 2nd basis
% for the anti-diagonal measurement choice of Bob
bob_meas_1_bra_basis_1 = bob_meas_1_ket_basis_1.';


% Compute the density matrices for the bases of
% the anti-diagonal measurement choice of Bob
bob_meas_1 = { bob_meas_1_ket_basis_0 * bob_meas_1_bra_basis_0, ...
               bob_meas_1_ket_basis_1 * bob_meas_1_bra_basis_1 };


% Compute the probability of
% Alice measuring the 1st basis,
% when choosing the (sigma_x) measurement
prob_alice_00 = trace( kron( alice_meas_0{1}, eye(2) ) * ...
                       rho_density_matrix );

% Compute the Observable Operator of
% Alice measuring the 1st basis,
% when choosing the (sigma_x) measurement
prob_alice_meas_00 = ( kron( alice_meas_0{1}, eye(2) ) * ...
                       rho_density_matrix * ...
                       kron( alice_meas_0{1}, eye(2) ) ) / ...
                     prob_alice_00;

% Compute the probability of
% Alice measuring the 2nd basis,
% when choosing the (sigma_x) measurement
prob_alice_01 = trace( kron( alice_meas_0{2}, eye(2) ) * ...
                       rho_density_matrix );

% Compute the Observable Operator of
% Alice measuring the 2nd basis,
% when choosing the (sigma_x) measurement
prob_alice_meas_01 = ( kron( alice_meas_0{2}, eye(2) ) * ...
                       rho_density_matrix * ...
                       kron( alice_meas_0{2}, eye(2) ) ) / ...
                     prob_alice_01;


% Compute the probability of
% Alice measuring the 1st basis,
% when choosing the (sigma_x) measurement,
% and Bob measuring the 1st basis,
% when choosing the diagonal measurement
prob_alice_meas_00_bob_meas_00 = ...
    trace( kron( eye(2), bob_meas_0{1} ) * prob_alice_meas_00 );

% Compute the probability of
% Alice measuring the 2nd basis,
% when choosing the (sigma_x) measurement,
% and Bob measuring the 2nd basis,
% when choosing the diagonal measurement
prob_alice_meas_01_bob_meas_01 = ...
    trace( kron( eye(2), bob_meas_0{2} ) * prob_alice_meas_01 );

% Compute the 1st probability for
% winning the game of the CHSH Inequality
% following a quantum strategy
prob_win_0 = prob_alice_00 * prob_alice_meas_00_bob_meas_00 + ...
             prob_alice_01 * prob_alice_meas_01_bob_meas_01;

% Compute the alpha value, considering
% the 1st probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
alpha_0 = ( prob_win_0 - ( 1 - prob_win_0 ) );

% Compute the 1st quantum bound,
% considering the alpha value from
% the 1st probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_0 = 4 * alpha_0;


% Compute the probability of
% Alice measuring the 1st basis,
% when choosing the (sigma_x) measurement,
% and Bob measuring the 1st basis,
% when choosing the anti-diagonal measurement
prob_alice_meas_00_bob_meas_10 = ...
    trace( kron( eye(2), bob_meas_1{1} ) * prob_alice_meas_00 );

% Compute the probability of
% Alice measuring the 2nd basis,
% when choosing the (sigma_x) measurement,
% and Bob measuring the 2nd basis,
% when choosing the anti-diagonal measurement
prob_alice_meas_01_bob_meas_11 = ...
    trace( kron( eye(2), bob_meas_1{2} ) * prob_alice_meas_01 );

% Compute the 2nd probability for
% winning the game of the CHSH Inequality
% following a quantum strategy
prob_win_1 = prob_alice_00 * prob_alice_meas_00_bob_meas_10 + ...
             prob_alice_01 * prob_alice_meas_01_bob_meas_11;

% Compute the alpha value, considering
% the 2nd probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
alpha_1 = ( prob_win_1 - ( 1 - prob_win_1 ) );

% Compute the 2nd quantum bound,
% considering the alpha value from
% the 2nd probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_1 = 4 * alpha_1;


% Compute the probability of
% Alice measuring the 1st basis,
% when choosing the (sigma_z) measurement
prob_alice_10 = trace( kron( alice_meas_1{1}, eye(2) ) * ...
                       rho_density_matrix );

% Compute the Observable Operator of
% Alice measuring the 1st basis,
% when choosing the (sigma_z) measurement
prob_alice_meas_10 = ( kron( alice_meas_1{1}, eye(2) ) * ...
                       rho_density_matrix * ...
                       kron( alice_meas_1{1}, eye(2) ) ) / ...
                     prob_alice_10;

% Compute the probability of
% Alice measuring the 2nd basis,
% when choosing the (sigma_z) measurement
prob_alice_11 = trace( kron( alice_meas_1{2}, eye(2) ) * ...
                       rho_density_matrix );

% Compute the Observable Operator of
% Alice measuring the 2nd basis,
% when choosing the (sigma_z) measurement
prob_alice_meas_11 = ( kron( alice_meas_1{2}, eye(2) ) * ...
                       rho_density_matrix * ...
                       kron( alice_meas_1{2}, eye(2) ) ) / ...
                     prob_alice_11;


% Compute the probability of
% Alice measuring the 1st basis,
% when choosing the (sigma_z) measurement,
% and Bob measuring the 1st basis,
% when choosing the diagonal measurement
prob_alice_meas_10_bob_meas_00 = ...
    trace( kron( eye(2), bob_meas_0{1} ) * prob_alice_meas_10 );

% Compute the probability of
% Alice measuring the 2nd basis,
% when choosing the (sigma_z) measurement,
% and Bob measuring the 2nd basis,
% when choosing the diagonal measurement
prob_alice_meas_11_bob_meas_01 = ...
    trace( kron( eye(2), bob_meas_0{2} ) * prob_alice_meas_11 );

% Compute the 3rd probability for
% winning the game of the CHSH Inequality
% following a quantum strategy
prob_win_2 = prob_alice_10 * prob_alice_meas_10_bob_meas_00 + ...
             prob_alice_11 * prob_alice_meas_11_bob_meas_01;

% Compute the alpha value, considering
% the 3rd probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
alpha_2 = ( prob_win_2 - ( 1 - prob_win_2 ) );

% Compute the 3rd quantum bound,
% considering the alpha value from
% the 3rd probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_2 = 4 * alpha_2;


% Compute the probability of
% Alice measuring the 1st basis,
% when choosing the (sigma_z) measurement,
% and Bob measuring the 1st basis,
% when choosing the anti-diagonal measurement
prob_alice_meas_10_bob_meas_10 = ...
    trace( kron( eye(2), bob_meas_1{1} ) * prob_alice_meas_10 );

% Compute the probability of
% Alice measuring the 2nd basis,
% when choosing the (sigma_z) measurement,
% and Bob measuring the 2nd basis,
% when choosing the anti-diagonal measurement
prob_alice_meas_11_bob_meas_11 = ...
    trace( kron( eye(2), bob_meas_1{2} ) * prob_alice_meas_11 );


% Compute the 4th probability for
% winning the game of the CHSH Inequality
% following a quantum strategy
prob_win_3 = prob_alice_10 * prob_alice_meas_10_bob_meas_10 + ...
             prob_alice_11 * prob_alice_meas_11_bob_meas_11;

% Compute the alpha value, considering
% the 4th probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
alpha_3 = ( prob_win_3 - ( 1 - prob_win_3 ) );

% Compute the 4th quantum bound,
% considering the alpha value from
% the 4th probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_3 = 4 * alpha_3;


% If all the quantum bounds computed,
% considering the alpha values estimated
% from the respective winning probabilities are
% equal to the quantum upper bound for the CHSH Inequality
if ( abs( chsh_inequality_quantum_upper_bound_L - ...
          quantum_bound_0 ) < error_tol ) && ...
   ( abs( chsh_inequality_quantum_upper_bound_L - ...
          quantum_bound_1) < error_tol ) && ...
   ( abs( chsh_inequality_quantum_upper_bound_L - ...
          quantum_bound_2) < error_tol ) && ...
   ( abs( chsh_inequality_quantum_upper_bound_L - ...
          quantum_bound_3) < error_tol )

    % Print of the estimation of the quantum
    % upper bound for the CHSH Inequality
    fprintf(['Quantum Upper Bound L ' ...
             'for CHSH Inequality:\n'])
    fprintf(['  e_00 + e_01 + e_10 - e_11 <= L^(Q) = ' ...
             '2 x sqrt(2) = %d\n'], ...
             chsh_inequality_quantum_upper_bound_L);
    
    % Print a blank line
    fprintf('\n');

end