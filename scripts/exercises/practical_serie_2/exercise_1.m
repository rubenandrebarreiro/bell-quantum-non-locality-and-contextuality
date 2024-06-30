% Clear all the values
% previously used in this scope
clear all; %#ok<CLALL>

% Clear Command Window
clc;


% Disable the warnings for the MATLAB
% graphics rendering features based on
% Open Graphics Language (OpenGL)
warning('off', 'MATLAB:hg:AutoSoftwareOpenGL');


% Define the number of thetas (in radians)
num_thetas_rad = 50;

% Create a linear space for the list of
% thetas (in radians) in [0, 2pi]
thetas_rad_linear_space = ...
    linspace(0, 2*pi, num_thetas_rad);


% Create the list of (optimized) expectation values
% for the Clauser-Horne-Shimony-Holt (CHSH) Inequality,
% corresponding to each theta (in radians), for Alice
optimized_chsh_expectation_values_alice = zeros(1, num_thetas_rad);

% Create the list of (optimized) expectation values
% for the Clauser-Horne-Shimony-Holt (CHSH) Inequality,
% corresponding to each theta (in radians), for Bob
optimized_chsh_expectation_values_bob = zeros(1, num_thetas_rad);


% Define the classical bound for
% the Clauser-Horne-Shimony-Holt (CHSH) Inequality
chsh_classical_bound = 2;

% Define the quantum bound for
% the Clauser-Horne-Shimony-Holt (CHSH) Inequality
chsh_quantum_bound = 2 * sqrt(2);

% Define the algebraic maximum bound for
% the Clauser-Horne-Shimony-Holt (CHSH) Inequality
chsh_algebraic_maximum_bound = 4;


% Define the maximum number of iterations for
% Seesaw (cyclic) Optimization algorithm/procedure
num_max_iterations = 20;


% Define the size of the Hilbert Space
hilbert_space_size = 2;


% Define the number of parties to be considered
% in the Clauser-Horne-Shimony-Holt (CHSH) Inequality
num_parties = 2;


% Define the number of inputs for Alice
num_alice_inputs = 2;

% Define the number of inputs for Bob
num_bob_inputs = 2;


% Define the number of outputs for Alice
num_alice_outputs = 2;

% Define the number of outputs for Bob
num_bob_outputs = 2;


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("****************************" + ...
        "****************************\n");

% Print two blank lines
fprintf("\n\n");



% Print an introductory message for Seesaw (cyclic) Optimization
% algorithm/procedure on the given set of thetas in [0, 2pi]
% for the Clauser-Horne-Shimony-Holt (CHSH) expectation value
fprintf("Starting the <strong>Seesaw (cyclic)</strong> " + ...
        "<strong>Optimization</strong>\n" + ...
        "algorithm/procedure on the given " + ...
        "set of <strong>thetas</strong>\n" + ...
        "in <strong>[0, 2pi]</strong> " + ...
        "for the <strong>Clauser-Horne-Shimony-Holt</strong> " + ...
        "<strong>(CHSH)</strong>\n" + ...
        "expectation value on " + ...
        "the <strong>Bell</strong> state defined as\n" + ...
        "<strong>|psi> = cos(theta)|00> + sin(theta)|11></strong>...\n\n");


% For the index of each theta (in radians)
for curr_theta_rad_idx = 1:num_thetas_rad
    
    % Retrieve the current theta (in radians)
    curr_theta_rad = ...
        thetas_rad_linear_space(curr_theta_rad_idx);
   
    % Print an introductory message about
    % the current theta (in radians) 
    fprintf(" * For <strong>theta_%d = %.8f</strong>:\n", ...
            curr_theta_rad_idx, curr_theta_rad);
    

    % Compute the Seesaw (cyclic) Optimization
    % algorithm/procedure on the current theta (in radians)
    % for the Clauser-Horne-Shimony-Holt (CHSH) expectation value
    [cvx_optval_alice, cvx_optval_bob, ...
     cvx_status_alice, cvx_status_bob, ...
     alice_measurement_projectors_optimized, ...
     bob_measurement_projectors_optimized] = ...
        compute_seesaw_optimization_for_chsh_inequality...
                (hilbert_space_size, num_parties, curr_theta_rad, ...
                 num_alice_inputs, num_alice_outputs, ...
                 num_bob_inputs, num_bob_outputs, ...
                 num_max_iterations, chsh_quantum_bound);
    
    % Print the Clauser-Horne-Shimony-Holt (CHSH)
    % expectation values for Alice and Bob, considering
    % the Seesaw (cyclic) Optimization algorithm/procedure
    % on the current theta (in radians)
    fprintf("   => <strong>S_(Alice) = %.8f</strong> | " + ...
                  "<strong>S_(Bob) = %.8f</strong>\n\n", ...
            cvx_optval_alice, cvx_optval_bob);


    optimized_chsh_expectation_values_alice(curr_theta_rad_idx) = ...
        cvx_optval_alice;


    optimized_chsh_expectation_values_bob(curr_theta_rad_idx) = ...
        cvx_optval_bob;

end


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("****************************" + ...
        "****************************\n");

% Print two blank lines
fprintf("\n\n");


% Plot the (optimized) Clauser-Horne-Shimony-Holt (CHSH)
% expectation values on the Einstein-Podolski-Rosen (EPR) Bell State,
% parameterized by a theta (in radians), for Alice and Bob
plot_optimized_chsh_exp_values_parameterized_bell_state...
         (thetas_rad_linear_space, ...
          optimized_chsh_expectation_values_alice, ...
          optimized_chsh_expectation_values_bob, ...
          chsh_classical_bound, chsh_quantum_bound, ...
          chsh_algebraic_maximum_bound);



% Define a function to compute the ket (column) vector
% for the Einstein-Podolski-Rosen (EPR) Bell State,
% parameterized by a theta (in radians)
function epr_bell_state_ket_vector = ...
         compute_epr_bell_state_ket_vector...
            (theta_rad, hilbert_space_size, num_parties)
    
    % Initialize a ket (column) vector
    % as a zero vector for the considered
    % size of the Hilbert Space and number of parties
    epr_bell_state_ket_vector = ...
        zeros(hilbert_space_size^num_parties, 1);
    
    % Define the parameterization for
    % the state |0...0> based on cosine function
    epr_bell_state_ket_vector(1) = ...
        cos(theta_rad);

    % Define the parameterization for
    % the state |1...1> based on sine function
    epr_bell_state_ket_vector...
        (hilbert_space_size^num_parties) = ...
            sin(theta_rad);

end


% Define a function to compute the bra (row) vector
% for the Einstein-Podolski-Rosen (EPR) Bell State,
% parameterized by a theta (in radians)
function epr_bell_state_bra_vector = ...
         compute_epr_bell_state_bra_vector...
            (theta_rad, hilbert_space_size, num_parties)
    
    % Compute the ket (column) vector
    % for the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by the given theta (in radians)
    epr_bell_state_ket_vector = ...
        compute_epr_bell_state_ket_vector...
            (theta_rad, hilbert_space_size, num_parties);
    
    % Compute the conjugate transpose of the ket (column)
    % vector for the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by the given theta (in radians),
    % to be the corresponding bra (row) vector
    epr_bell_state_bra_vector = ...
        epr_bell_state_ket_vector';

end


% Define a function to compute the density matrix (rho)
% for the Einstein-Podolski-Rosen (EPR) Bell State,
% parameterized by a theta (in radians)
function epr_bell_state_rho_density_matrix = ...
         compute_epr_bell_state_rho_density_matrix...
            (theta_rad, hilbert_space_size, num_parties)
    
    % Compute the ket (column) vector
    % for the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by a theta (in radians)
    epr_bell_state_ket_vector = ...
        compute_epr_bell_state_ket_vector...
            (theta_rad, hilbert_space_size, num_parties);
    
    % Compute the bra (row) vector
    % for the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by a theta (in radians)
    epr_bell_state_bra_vector = ...
        compute_epr_bell_state_bra_vector...
            (theta_rad, hilbert_space_size, num_parties);
    
    % Compute the density matrix (rho)
    % for the Einstein-Podolski-Rosen (EPR)
    % Bell State, parameterized by a theta (in radians),
    % as the multiplication between the respective
    % ket (column) and bra (row) vectors
    epr_bell_state_rho_density_matrix = ...
        epr_bell_state_ket_vector * epr_bell_state_bra_vector;

end


% Define a function to compute the (optimized)
% measurement projectors on the Einstein-Podolski-Rosen (EPR)
% Bell State, parameterized by a theta (in radians), for Alice
function [cvx_optval, cvx_status, ...
          alice_measurement_projectors_optimized] = ...
         compute_alice_measurement_operator_semi_definite_optimization...
            (num_alice_inputs, num_alice_outputs, ...
             num_bob_inputs, num_bob_outputs, ...
             hilbert_space_size, ...
             epr_bell_state_rho_density_matrix, ...
             bob_measurement_projectors) %#ok<STOUT>

    
    % Initialize the (optimized) measurement projectors
    % on the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by a theta (in radians), for Alice,
    % as a zero matrix
    alice_measurement_projectors_optimized = ...
        zeros(hilbert_space_size, hilbert_space_size, ...
              num_alice_inputs, num_alice_outputs);
    
    
    % Start a code block for the ConVeX (CVX)
    % optimization programming framework,
    % in quiet mode (no verbose outputs),
    % using Semi-Definite Programming (SDP)
    cvx_begin sdp quiet
        
        % Define the Self-Dual-Minimization (SeDuMi) solver
        % for the ConVeX (CVX) optimization programming framework 
        cvx_solver sedumi

        % Defin the best precision for the ConVeX (CVX)
        % optimization programming framework 
        cvx_precision best
        

        % Define the Alice's measurement
        % projectors as the variable to be optimized by
        % the ConVeX (CVX) optimization programming framework 
        variable alice_measurement_projectors_optimized(hilbert_space_size, hilbert_space_size, num_alice_inputs, num_alice_outputs) ...
                 hermitian semidefinite
        

        % Initialize the list of
        % Alice's measurement operators
        alice_measurement_operators = ...
            cell(hilbert_space_size, 1);

        % Initialize the list of
        % Bob's measurement operators
        bob_measurement_operators = ...
            cell(hilbert_space_size, 1);
        
    
        % For each output for Alice
        for curr_alice_output = 1:num_alice_outputs
            
            % Define the Alice's measurement operator
            % regarding the current output for Alice
            alice_measurement_operators{curr_alice_output} = ...
                alice_measurement_projectors_optimized...
                    (:, :, 1, curr_alice_output) - ...
                alice_measurement_projectors_optimized...
                    (:, :, 2, curr_alice_output);
    
        end
    

        % For each output for Bob
        for curr_bob_output = 1:num_bob_outputs
            
            % Define the Bob's measurement operator
            % regarding the current output for Bob
            bob_measurement_operators{curr_bob_output} = ...
                bob_measurement_projectors...
                    (:, :, 1, curr_bob_output) - ...
                bob_measurement_projectors...
                    (:, :, 2, curr_bob_output);
            
        end    
        
        
        % Define the Clauser-Horne-Shimony-Holt (CHSH)
        % expetation value as the objective function to
        % be optimized by the ConVeX (CVX) optimization
        % programming framework 
        chsh_expectation_value_objective_function = ...
            trace( kron( alice_measurement_operators{1}, ...
                         bob_measurement_operators{1} ) * ...
                   epr_bell_state_rho_density_matrix ) + ...
            trace( kron( alice_measurement_operators{1}, ...
                         bob_measurement_operators{2} ) * ...
                   epr_bell_state_rho_density_matrix ) + ...
            trace( kron( alice_measurement_operators{2}, ...
                         bob_measurement_operators{1} ) * ...
                   epr_bell_state_rho_density_matrix ) - ...
            trace( kron( alice_measurement_operators{2}, ...
                         bob_measurement_operators{2} ) * ...
                   epr_bell_state_rho_density_matrix );
        

        % Define the optimization problem as
        % a maximization to be optimized by
        % the ConVeX (CVX) optimization programming framework 
        maximise( real(chsh_expectation_value_objective_function) );
        

        % Create a list for the sums of
        % the Alice's measurement projectors
        alice_measurement_projectors_sums = ...
            cell(num_alice_outputs, 1);
        
        % Create a list for the sums of
        % the Bob's measurement projectors
        bob_measurement_projectors_sums = ...
            cell(num_bob_outputs, 1);
 
        
        % For each output for Alice
        for curr_alice_output = 1:num_alice_outputs
                
            % Initialize the sum for the Alice's measurement
            % projectors, regarding the current output
            curr_alice_measurement_projectors_sum = ...
                zeros(hilbert_space_size, hilbert_space_size);
            

            % For each input for Alice
            for curr_alice_input = 1:num_alice_inputs
                
                % Sum the Alice's measurement
                % projector for her current input,
                % to the sum of her measurement projectors,
                % regarding the current output
                curr_alice_measurement_projectors_sum = ...
                        curr_alice_measurement_projectors_sum + ...
                            alice_measurement_projectors_optimized...
                                (:, :, curr_alice_input, ...
                                       curr_alice_output);

            end
            

            % Define the (final) sum for the Alice's
            % measurement projectors, regarding the current output 
            alice_measurement_projectors_sums{curr_alice_output} = ...
                curr_alice_measurement_projectors_sum;
            
        end
        

        % For each output for Bob        
        for curr_bob_output = 1:num_bob_outputs
                
            % Initialize the sum for the Bob's measurement
            % projectors, regarding the current output
            curr_bob_measurement_projectors_sum = ...
                zeros(hilbert_space_size, hilbert_space_size);
            
            % For each input for Bob
            for curr_bob_input = 1:num_bob_inputs
                
                % Sum the Bob's measurement
                % projector for her current input,
                % to the sum of her measurement projectors,
                % regarding the current output
                curr_bob_measurement_projectors_sum = ...
                        curr_bob_measurement_projectors_sum + ...
                            bob_measurement_projectors...
                                (:, :, curr_bob_input, ...
                                       curr_bob_output);

            end


            % Define the (final) sum for the Bob's
            % measurement projectors, regarding the current output 
            bob_measurement_projectors_sums{curr_bob_output} = ...
                curr_bob_measurement_projectors_sum;
            
        end
        
    
        % Define the conditions to which the maximization
        % optimization problem for the ConVeX (CVX) optimization
        % programming framework is subject to
        subject to
            
            % For each output for Alice
            for curr_alice_output = 1:num_alice_outputs
                
                % If the sum for the Alice's measurement
                % projectors, regarding the current output,
                % is equal to the identity matrix
                alice_measurement_projectors_sums...
                    {curr_alice_output} == eye(hilbert_space_size); %#ok<EQEFF>

            end
            

            % For each output for Bob
            for curr_bob_output = 1:num_bob_outputs
                
                % If the sum for the Bob's measurement
                % projectors, regarding the current output,
                % is equal to the identity matrix
                bob_measurement_projectors_sums...
                    {curr_bob_output} == eye(hilbert_space_size); %#ok<EQEFF>

            end

    cvx_end
    
end


% Define a function to compute the (optimized)
% measurement projectors on the Einstein-Podolski-Rosen (EPR)
% Bell State, parameterized by a theta (in radians), for Bob
function [cvx_optval, cvx_status, ...
          bob_measurement_projectors_optimized] = ...
         compute_bob_measurement_operator_semi_definite_optimization...
            (num_bob_inputs, num_bob_outputs, ...
             num_alice_inputs, num_alice_outputs, ...
             hilbert_space_size, ...
             epr_bell_state_rho_density_matrix, ...
             alice_measurement_projectors) %#ok<STOUT>
    
    % Initialize the (optimized) measurement projectors
    % on the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by a theta (in radians), for Alice,
    % as a zero matrix
    bob_measurement_projectors_optimized = ...
        zeros(hilbert_space_size, hilbert_space_size, ...
              num_bob_inputs, num_bob_outputs);
    
    
    % Start a code block for the ConVeX (CVX)
    % optimization programming framework,
    % in quiet mode (no verbose outputs),
    % using Semi-Definite Programming (SDP)
    cvx_begin sdp quiet
    
        % Define the Self-Dual-Minimization (SeDuMi) solver
        % for the ConVeX (CVX) optimization programming framework 
        cvx_solver sedumi
        
        % Defin the best precision for the ConVeX (CVX)
        % optimization programming framework 
        cvx_precision best
        
        
        % Define the Alice's measurement
        % projectors as the variable to be optimized by
        % the ConVeX (CVX) optimization programming framework 
        variable bob_measurement_projectors_optimized(hilbert_space_size, hilbert_space_size, num_bob_inputs, num_bob_outputs) ...
                 hermitian semidefinite
        

        % Initialize the list of
        % Alice's measurement operators
        alice_measurement_operators = ...
            cell(hilbert_space_size, 1);

        % Initialize the list of
        % Bob's measurement operators
        bob_measurement_operators = ...
            cell(hilbert_space_size, 1);
    
    
        % For each output for Alice
        for curr_alice_output = 1:num_alice_outputs
            
            % Define the Alice's measurement operator
            % regarding the current output for Alice
            alice_measurement_operators{curr_alice_output} = ...
                alice_measurement_projectors...
                    (:, :, 1, curr_alice_output) - ...
                alice_measurement_projectors...
                    (:, :, 2, curr_alice_output);
    
        end
        
    
        % For each output for Bob
        for curr_bob_output = 1:num_bob_outputs
            
            % Define the Bob's measurement operator
            % regarding the current output for Bob
            bob_measurement_operators{curr_bob_output} = ...
                bob_measurement_projectors_optimized...
                    (:, :, 1, curr_bob_output) - ...
                bob_measurement_projectors_optimized...
                    (:, :, 2, curr_bob_output);
            
        end        


        % Define the Clauser-Horne-Shimony-Holt (CHSH)
        % expetation value as the objective function to
        % be optimized by the ConVeX (CVX) optimization
        % programming framework 
        chsh_expectation_value_objective_function = ...
            trace( kron( alice_measurement_operators{1}, ...
                         bob_measurement_operators{1} ) * ...
                   epr_bell_state_rho_density_matrix ) + ...
            trace( kron( alice_measurement_operators{1}, ...
                         bob_measurement_operators{2} ) * ...
                   epr_bell_state_rho_density_matrix ) + ...
            trace( kron( alice_measurement_operators{2}, ...
                         bob_measurement_operators{1} ) * ...
                   epr_bell_state_rho_density_matrix ) - ...
            trace( kron( alice_measurement_operators{2}, ...
                         bob_measurement_operators{2} ) * ...
                   epr_bell_state_rho_density_matrix );


        % Define the optimization problem as
        % a maximization to be optimized by
        % the ConVeX (CVX) optimization programming framework 
        maximise( real(chsh_expectation_value_objective_function) )
        

        % Create a list for the sums of
        % the Alice's measurement projectors
        alice_measurement_projectors_sums = ...
            cell(num_alice_outputs, 1);
        
        % Create a list for the sums of
        % the Bob's measurement projectors
        bob_measurement_projectors_sums = ...
            cell(num_bob_outputs, 1);
        

        % For each output for Alice
        for curr_alice_output = 1:num_alice_outputs
                
            % Initialize the sum for the Alice's measurement
            % projectors, regarding the current output
            curr_alice_measurement_projectors_sum = ...
                zeros(hilbert_space_size, hilbert_space_size);
            
            
            % For each input for Alice
            for curr_alice_input = 1:num_alice_inputs
                
                % Sum the Alice's measurement
                % projector for her current input,
                % to the sum of her measurement projectors,
                % regarding the current output
                curr_alice_measurement_projectors_sum = ...
                        curr_alice_measurement_projectors_sum + ...
                            alice_measurement_projectors...
                                (:, :, curr_alice_input, ...
                                       curr_alice_output);

            end

            
            % Define the (final) sum for the Alice's
            % measurement projectors, regarding the current output 
            alice_measurement_projectors_sums{curr_alice_output} = ...
                curr_alice_measurement_projectors_sum;
            
        end
        

        % For each output for Bob
        for curr_bob_output = 1:num_bob_outputs
                
            % Initialize the sum for the Bob's measurement
            % projectors, regarding the current output
            curr_bob_measurement_projectors_sum = ...
                zeros(hilbert_space_size, hilbert_space_size);
            

            % For each input for Bob
            for curr_bob_input = 1:num_bob_inputs
                
                % Sum the Bob's measurement
                % projector for her current input,
                % to the sum of her measurement projectors,
                % regarding the current output
                curr_bob_measurement_projectors_sum = ...
                        curr_bob_measurement_projectors_sum + ...
                            bob_measurement_projectors_optimized...
                                (:, :, curr_bob_input, ...
                                       curr_bob_output);

            end


            % Define the (final) sum for the Bob's
            % measurement projectors, regarding the current output 
            bob_measurement_projectors_sums{curr_bob_output} = ...
                curr_bob_measurement_projectors_sum;
            
        end

        
    
        % Define the conditions to which the maximization
        % optimization problem for the ConVeX (CVX) optimization
        % programming framework is subject to
        subject to
            
            % For each output for Alice
            for curr_alice_output = 1:num_alice_outputs

                % If the sum for the Alice's measurement
                % projectors, regarding the current output,
                % is equal to the identity matrix
                alice_measurement_projectors_sums...
                    {curr_alice_output} == eye(hilbert_space_size); %#ok<EQEFF>
                
            end
            
            
            % For each output for Bob
            for curr_bob_output = 1:num_bob_outputs
    
                % If the sum for the Bob's measurement
                % projectors, regarding the current output,
                % is equal to the identity matrix
                bob_measurement_projectors_sums...
                    {curr_bob_output} == eye(hilbert_space_size); %#ok<EQEFF>
                
            end

    cvx_end

end


% Define a function to compute the (optimized)
% measurement projectors on the Einstein-Podolski-Rosen (EPR)
% Bell State, parameterized by a theta (in radians), for Alice and Bob,
% executing the Seesaw (cyclic) optimization algorithm/procedure
function [cvx_optval_alice, cvx_optval_bob, ...
          cvx_status_alice, cvx_status_bob, ...
          alice_measurement_projectors_optimized, ...
          bob_measurement_projectors_optimized] = ...
            compute_seesaw_optimization_for_chsh_inequality...
                (hilbert_space_size, num_parties, theta_rad, ...
                 num_alice_inputs, num_alice_outputs, ...
                 num_bob_inputs, num_bob_outputs, ...
                 num_max_iterations, chsh_quantum_bound)
    
    % Compute the density matrix (rho)
    % for the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by a theta (in radians)
    epr_bell_state_rho_density_matrix = ...
        compute_epr_bell_state_rho_density_matrix...
            (theta_rad, hilbert_space_size, num_parties);
    
    
    % Create the list for all the (optimized)
    % Alice's measurement projectors, regarding
    % his inputs and outputs
    alice_measurement_projectors_optimized = ...
        zeros(hilbert_space_size, hilbert_space_size, ...
              num_alice_inputs, num_alice_outputs);

    
    % For each output for Alice
    for curr_alice_output = 1:num_alice_outputs
        
        % Generate a random density matrix from
        % which will be built the Alice's measurement projectors
        random_alice_projector = ...
            RandomDensityMatrix(hilbert_space_size);

        % Compute the 1st (optimized) Alice's
        % measurement projectors, for her current output
        alice_measurement_projectors_optimized...
            (:, :, 1, curr_alice_output) = ...
                random_alice_projector;

        % Compute the 2nd (optimized) Alice's
        % measurement projectors, for her current output
        alice_measurement_projectors_optimized...
            (:, :, 2, curr_alice_output) = ...
                eye(hilbert_space_size) - random_alice_projector;
        
    end


    % Create the list for all the (optimized)
    % Bob's measurement projectors, regarding
    % his inputs and outputs
    bob_measurement_projectors_optimized = ...
        zeros(hilbert_space_size, hilbert_space_size, ...
              num_bob_inputs, num_bob_outputs);


    % For each output for Bob
    for curr_bob_output = 1:num_bob_outputs
        
        % Generate a random density matrix from
        % which will be built the Bob's measurement projectors
        random_bob_projector = ...
            RandomDensityMatrix(hilbert_space_size);

        % Compute the 1srt (optimized) Bob's
        % measurement projectors, for his current output
        bob_measurement_projectors_optimized...
            (:, :, 1, curr_bob_output) = ...
                random_bob_projector;

        % Compute the 2nd (optimized) Bob's
        % measurement projectors, for his current output
        bob_measurement_projectors_optimized...
            (:, :, 2, curr_bob_output) = ...
                eye(hilbert_space_size) - random_bob_projector;
        
    end
    

    
    % Define the (initial) current iteration
    curr_iteration = 1;

    
    % Start an 'infinite loop'
    while(true)
        
        % Compute the (optimized) measurement projectors
        % on the Einstein-Podolski-Rosen (EPR) Bell State,
        % parameterized by a theta (in radians), for Alice
        [cvx_optval_alice, cvx_status_alice, ...
         alice_measurement_projectors_optimized] = ...
             compute_alice_measurement_operator_semi_definite_optimization...
                (num_alice_inputs, num_alice_outputs, ...
                 num_bob_inputs, num_bob_outputs, ...
                 hilbert_space_size, ...
                 epr_bell_state_rho_density_matrix, ...
                 bob_measurement_projectors_optimized);
        
        % Compute the (optimized) measurement projectors
        % on the Einstein-Podolski-Rosen (EPR) Bell State,
        % parameterized by a theta (in radians), for Bob
        [cvx_optval_bob, cvx_status_bob, ...
         bob_measurement_projectors_optimized] = ...
             compute_bob_measurement_operator_semi_definite_optimization...
                (num_bob_inputs, num_bob_outputs, ...
                 num_alice_inputs, num_alice_outputs, ...
                 hilbert_space_size, ...
                 epr_bell_state_rho_density_matrix, ...
                 alice_measurement_projectors_optimized);
        
        
        % If the current interation reaches
        % the maximum number of iterations or 
        % Alice's and Bob's expectation values for
        % the Clauser-Horne-Shimony-Holt (CHSH) Inequality reach
        % their optimal value, with some very low tolerance error
        if( curr_iteration == num_max_iterations || ...
            ( ismembertol(cvx_optval_alice, ...
                          chsh_quantum_bound, 1e-6) && ...
              ismembertol(cvx_optval_bob, ...
                          chsh_quantum_bound, 1e-6) ) )

            % Return out of the function,
            % breaking the 'infinite loop'
            return;

        end


        % Increment the current iteration to the next one
        curr_iteration = curr_iteration + 1;

    end

end


% Define a function to plot the (optimized)
% Clauser-Horne-Shimony-Holt (CHSH) expectation values
% on the Einstein-Podolski-Rosen (EPR) Bell State,
% parameterized by a theta (in radians), for Alice and Bob
function plot_optimized_chsh_exp_values_parameterized_bell_state...
         (thetas_rad, optimized_chsh_expectation_values_alice, ...
          optimized_chsh_expectation_values_bob, ...
          chsh_classical_bound, chsh_quantum_bound, ...
          chsh_algebraic_maximum)
    
    % Define the position of the window
    % for the current figure, extending
    % its size to 1000 x 800 
    set(gcf, "position", [40, 40, 1040, 840]);
    

    % Plot the set of all the (optimized)
    % Clauser-Horne-Shimony-Holt (CHSH) expectation values
    % on the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by a theta (in radians), for Alice
    plot(thetas_rad, optimized_chsh_expectation_values_alice, ...
         "-s", "Color", [0.8500 0.3250 0.0980]);
    
    % Define the limits for the x-axis to
    % be within the range [0, 2pi]
    xlim([0, 2*pi]);
    
    % Define the limits for the y-axis to
    % be within the range [1.75, 4.25]
    ylim([1.75, 4.25]);


    % Define the title for the plot graphics
    title("Optimization for CHSH Expectation Value S");


    % Define the label for the x-axis of the plot graphics
    xlabel("theta (in radians)");

    % Define the label for the y-axis of the plot graphics
    ylabel("S (CHSH Expectation Values)");

    
    % Hold on and retain the previous plot,
    % allowing to add new ones
    hold on;
    

    % Plot the set of all the (optimized)
    % Clauser-Horne-Shimony-Holt (CHSH) expectation values
    % on the Einstein-Podolski-Rosen (EPR) Bell State,
    % parameterized by a theta (in radians), for Bob 
    plot(thetas_rad, optimized_chsh_expectation_values_bob, ...
         "-s", "Color", [0.4940 0.1840 0.5560]);


    % Plot the horizontal line for the classical bound
    % for the Clauser-Horne-Shimony-Holt (CHSH) Inequality
    line([0, 2*pi], [chsh_classical_bound, chsh_classical_bound], ...
         "Color", "red", "LineStyle", "-", "LineWidth", 2);

    % Plot the horizontal line for the classical bound
    % for the Clauser-Horne-Shimony-Holt (CHSH) Inequality
    line([0, 2*pi], [chsh_quantum_bound, chsh_quantum_bound], ...
         "Color", "green", "LineStyle", "-", "LineWidth", 2);

    % Plot the horizontal line for the classical bound
    % for the Clauser-Horne-Shimony-Holt (CHSH) Inequality
    line([0, 2*pi], [chsh_algebraic_maximum, chsh_algebraic_maximum], ...
         "Color", "blue", "LineStyle", "-", "LineWidth", 2);

    
    % Add the legends for all the curves
    % and lines plotted previously
    legend("S_{(Alice)}", "S_{(Bob)}", ...
           "Classical Bound ({S}_{C} <= 2)", ...
           "Quantum Bound ({S}_{Q} <= 2 * sqrt(2))", ...
           "Algebraic Maximum Bound ({S}_{A-MAX} <= 4)");

    % Hold off and stop retain all the previous plots
    hold off;
    

    % Add a grid to the figure
    % containing all the plots
    grid;

end
