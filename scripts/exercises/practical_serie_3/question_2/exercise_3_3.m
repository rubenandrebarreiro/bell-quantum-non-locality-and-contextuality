% Clear Command Window
clc;


% Create the ket vector for the two-dimensional
% quantum Einstein-Podolski-Rosen (EPR) pair denoted
% as |{phi}^{+}> = ( 1 / sqrt(2) ) * ([1, 0, 0, 1])^T
ket_epr_pair_phi_plus = ( 1 / sqrt(2) ) * [ 1, 0, 0, 1 ].';


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


% Define the number of constraints
% for the Hardy's Paradox regarding
% the bipartite Bell experiment and
% the CHSH Inequality
num_hardy_paradox_constraints = 5;


% Define the array with the probability
% constraints for the Hardy's Paradox
% regarding the bipartite Bell experiment
% and the CHSH Inequality
probabilities_hardy_paradox_constraints = ...
    zeros(1, num_hardy_paradox_constraints);


% 
alice_bob_meas_00_ket_epr_pair_phi_plus_evol = ...
    kron(alice_measurements{1}, bob_measurements{1}) * ...
        ket_epr_pair_phi_plus;

% 
alice_bob_meas_01_ket_epr_pair_phi_plus_evol = ...
    kron(alice_measurements{1}, bob_measurements{2}) * ...
        ket_epr_pair_phi_plus;

% 
alice_bob_meas_10_ket_epr_pair_phi_plus_evol = ...
    kron(alice_measurements{2}, bob_measurements{1}) * ...
        ket_epr_pair_phi_plus;

% 
alice_bob_meas_11_ket_epr_pair_phi_plus_evol = ...
    kron(alice_measurements{2}, bob_measurements{2}) * ...
        ket_epr_pair_phi_plus;


probabilities_hardy_paradox_constraints(1) = ...
    abs(alice_bob_meas_00_ket_epr_pair_phi_plus_evol(1))^2;

probabilities_hardy_paradox_constraints(2) = ...
    abs(alice_bob_meas_00_ket_epr_pair_phi_plus_evol(4))^2;

probabilities_hardy_paradox_constraints(3) = ...
    abs(alice_bob_meas_01_ket_epr_pair_phi_plus_evol(4))^2;

probabilities_hardy_paradox_constraints(4) = ...
    abs(alice_bob_meas_10_ket_epr_pair_phi_plus_evol(4))^2;

probabilities_hardy_paradox_constraints(5) = ...
    abs(alice_bob_meas_10_ket_epr_pair_phi_plus_evol(4))^2;


% Print the headline for the probability constraints
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf(['* The probalility constraints for \n' ...
         '  the Hardy Paradox are the following ones:\n']);
fprintf('\n');


% Print the 1st probability constraint
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf('   1) P_{win}(a_x, b_y = 00 | x, y = 00) = 0.0\n');

% Print the 2nd probability constraint
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf('   2) P_{win}(a_x, b_y = 11 | x, y = 00) = 0.0\n');

% Print the 3rd probability constraint
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf('   3) P_{win}(a_x, b_y = 11 | x, y = 01) = 0.0\n');

% Print the 4th probability constraint
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf('   4) P_{win}(a_x, b_y = 11 | x, y = 10) = 0.0\n');

% Print the 5th probability constraint
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf('   5) P_{win}(a_x, b_y = 11 | x, y = 11) > 0.0\n');


% Print a blank line
fprintf('\n');

% Print a blank line
fprintf('\n');


% Print the headline for the probability constraints
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf(['* The checking of the probalility constraints\n' ...
         '  for the Hardy Paradox stands as follows:\n']);
fprintf('\n');


% If the 1st probability constraint holds
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(1) == 0.0

    % Print some information about
    % the 1st probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('   1) P_{win}(a_x, b_y = 00 | x, y = 00) = 0.0 holds!\n');

    % Print a blank line
    fprintf('\n');

end


% If the 1st probability constraint does not hold
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(1) ~= 0.0
    
    % Print some information about
    % the 1st probability constraint not holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf(['   1) P_{win}(a_x, b_y = 00 | x, y = 00) = 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 1st probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 00 | x, y = 00) = %.4f\n', ...
            probabilities_hardy_paradox_constraints(1));
    
    % Print a blank line
    fprintf('\n');

end


% If the 2nd probability constraint holds
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(2) == 0.0

    % Print some information about
    % the 2nd probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('   2) P_{win}(a_x, b_y = 11 | x, y = 00) = 0.0 holds!\n');

    % Print a blank line
    fprintf('\n');

end


% If the 2nd probability constraint does not hold
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(2) ~= 0.0
    
    % Print some information about
    % the 2nd probability constraint not holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf(['   2) P_{win}(a_x, b_y = 11 | x, y = 00) = 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 2nd probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 11 | x, y = 00) = %.4f\n', ...
            probabilities_hardy_paradox_constraints(2));
    
    % Print a blank line
    fprintf('\n');

end


% If the 3rd probability constraint holds
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(3) == 0.0

    % Print some information about
    % the 3rd probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('   3) P_{win}(a_x, b_y = 11 | x, y = 01) = 0.0 holds!\n');
    
    % Print a blank line
    fprintf('\n');

end


% If the 3rd probability constraint does not hold
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(3) ~= 0.0
    
    % Print some information about
    % the 3rd probability constraint not holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf(['   3) P_{win}(a_x, b_y = 11 | x, y = 01) = 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 3rd probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 11 | x, y = 01) = %.4f\n', ...
            probabilities_hardy_paradox_constraints(2));

    % Print a blank line
    fprintf('\n');

end


% If the 4th probability constraint holds
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(4) == 0.0

    % Print some information about
    % the 4th probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('   4) P_{win}(a_x, b_y = 11 | x, y = 10) = 0.0 holds!\n');
    
    % Print a blank line
    fprintf('\n');

end


% If the 4th probability constraint does not hold
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(4) ~= 0.0
    
    % Print some information about
    % the 4th probability constraint not holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf(['   4) P_{win}(a_x, b_y = 11 | x, y = 10) = 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 4th probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 11 | x, y = 10) = %.4f\n', ...
            probabilities_hardy_paradox_constraints(4));

    % Print a blank line
    fprintf('\n');

end


% If the 5th probability constraint holds
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(5) > 0.0

    % Print some information about
    % the 5th probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('   5) P_{win}(a_x, b_y = 11 | x, y = 11) > 0.0 holds!\n');
    
    % Print a blank line
    fprintf('\n');

end


% If the 5th probability constraint does not hold
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(5) <= 0.0
    
    % Print some information about
    % the 5th probability constraint not holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf(['   5) P_{win}(a_x, b_y = 11 | x, y = 11) > 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 5th probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 11 | x, y = 11) = %.4f\n', ...
            probabilities_hardy_paradox_constraints(5));

    % Print a blank line
    fprintf('\n');

end
