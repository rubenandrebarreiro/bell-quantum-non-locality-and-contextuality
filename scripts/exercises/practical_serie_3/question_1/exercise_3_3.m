% Clear Command Window
clc;


% Create the Pauli X (sigma_x) Matrix
pauli_sigma_x = full(Pauli('X'));

% Create the Pauli Y (sigma_y) Matrix
pauli_sigma_y = full(Pauli('Y'));


% Define the operator A^(n)_0
a_n_0_operator = pauli_sigma_y;

% Define the operator A^(n)_1
a_n_1_operator = -pauli_sigma_x;


% Print the content of the operator A^(n)_0
fprintf("The operator A^(n)_0 to be used is: \n" + ...
        "* A^(n)_0 = sigma_y = \n");
disp(a_n_0_operator);


% Print the content of the operator A^(n)_1
fprintf("The operator A^(n)_1 to be used is: \n" + ...
        "* A^(n)_1 = -sigma_x = \n");
disp(a_n_1_operator);

% Print a separator line
fprintf("**********************" + ...
        "**********************\n\n");


% Define the Mermin Inequality local classical bound
mermin_inequality_local_classical_bound = vpa(2);


% Define the minimum number n of parties
num_min_parties = 3;

% Define the maximum number n of parties
num_max_parties = 14;


% For each configuration regarding the number n of parties
for num_parties = num_min_parties:num_max_parties
    
    % If the current configuration regarding
    % the number n of parties is odd
    if( mod(num_parties, 2) ~= 0 )
            
        % Print a headline for the current configuration
        % regarding the number n of parties
        fprintf("For a number of parties n = %d:\n\n", num_parties);
    

        % For each value of probability v to
        % configure the noisy GHZ state
        for probability_v = 0:0.025:1.0

            % Compute the expectation value
            % for the Mermin Inequality <M_n>,
            % considering a number n of parties
            % as well as, a probability v
            m_n_werner_noisy_ghz_state_expectation_value = ...
                compute_m_n_werner_noisy_ghz_state_expectation_value...
                    (num_parties, probability_v, ...
                     a_n_0_operator, a_n_1_operator);
            
            % Compute and print the expectation value
            % <M_n> for the current configuration regarding
            % the number n of parties and probability v
            fprintf("* The expectation value <M_%d> for\n " + ...
                    " the operator M_%d and a noisy GHZ state\n " + ...
                    " with n = %d parties and v = %.4f is: <M_%d> = %.4f\n", ...
                    num_parties, num_parties, ...
                    num_parties, probability_v, num_parties, ...
                    m_n_werner_noisy_ghz_state_expectation_value);


            % If the expectation value <M_n>
            % for the current configuration regarding
            % the number n of parties and probability v
            % violates the Mermin Inequality (local classical bound)
            if( m_n_werner_noisy_ghz_state_expectation_value >= ...
                mermin_inequality_local_classical_bound )
                
                % Print the expectation value <M_n>
                % for the current configuration regarding
                % the number n of parties and probability v
                % violating the Mermin Inequality (local classical bound)
                fprintf("* <M_%d> [n = %d ; v = %.4f] >= C_(%d)" + ...
                        " (=) %.4f >= %.4f\n", ...
                        num_parties, num_parties, ...
                        probability_v, num_parties, ...
                        m_n_werner_noisy_ghz_state_expectation_value, ...
                        mermin_inequality_local_classical_bound);

                % Print an informative message about
                % the expectation value <M_n> violating with
                % the Mermin Inequality (local classical bound) C_(n)
                % for the current configuration regarding
                % the number n of parties and probability v
                fprintf("* The expectation value <M_%d> for\n " + ...
                        " the operator M_%d and a noisy GHZ state\n " + ...
                        " with n = %d parties and v = %.4f violates\n " + ...
                        " the Mermin Inequality (local classical bound)!\n", ...
                        num_parties, num_parties, num_parties, ...
                        probability_v);
            
            end
            

            % If the expectation value <M_n>
            % for the current configuration regarding
            % the number n of parties and probability v does not
            % violate the Mermin Inequality (local classical bound)
            if( m_n_werner_noisy_ghz_state_expectation_value < ...
                mermin_inequality_local_classical_bound )
                
                % Print the expectation value <M_n>
                % for the current configuration regarding
                % the number n of parties and probability v not
                % violating the Mermin Inequality (local classical bound)
                fprintf("* <M_%d> [n = %d ; v = %.4f] < C_(%d)" + ...
                        " (=) %.4f < %.4f\n", ...
                        num_parties, num_parties, ...
                        probability_v, num_parties, ...
                        m_n_werner_noisy_ghz_state_expectation_value, ...
                        mermin_inequality_local_classical_bound);
                
            end

            % Print two blank lines
            fprintf("\n\n");

        end

    end

end


% Print a separator
fprintf("**********************" + ...
        "**********************\n\n");

% Print a section about the conclusions of this quantum experiment,
% regarding the configuration for a noisy GHZ state with (white)
% noise for a given number n of parties and a probability v
fprintf("Conclusions:\n")
fprintf("1) Greater the odd number of parties, easier to\n" + ...
        "   achieve the quantum bound Q to violate\n" + ...
        "   the Mermin Inequality (local classical bound) C;\n");
fprintf("2) The (white) noise will affect the GHZ state,\n" + ...
        "   regarding the violation of the Mermin Inequality\n" + ...
        "   (local classical bound) C;\n")
fprintf("3) Greater the odd number of parties,\n" + ...
        "   better tolerance to (white) noise, the GHZ state\n" + ...
        "   will have, regarding the violation of\n" + ...
        "   the Mermin Inequality (local classical bound) C;\n")
fprintf("4) Everytime we increase the number of parties to\n" + ...
        "   the next odd number n, the minimum probability v\n" + ...
        "   that is required to ensure the violation of\n" + ...
        "   the Mermin Inequality (local classical bound) C, halves;\n\n")

% Print a separator
fprintf("**********************" + ...
        "**********************\n\n");


% Define a function to compute the ket (column) vector
% for the GHZ state, given a number n of parties
function ghz_state_ket_vector = ...
         compute_ghz_state_ket_vector(n)
    
    % Create a GHZ state for the given number n of parties
    ghz_state_ket_vector = zeros(2^n, 1);
    
    % Define the first element of
    % the GHZ state as 1 / sqrt(2)
    ghz_state_ket_vector(1,1) = 1/sqrt(2);
    
    % Define the last element of
    % the GHZ state as 1 / sqrt(2)
    ghz_state_ket_vector(2^n,1) = 1/sqrt(2);

end


% Define a function to compute the bra (row) vector
% for the GHZ state, given a number n of parties
function ghz_state_bra_vector = ...
         compute_ghz_state_bra_vector(n)
    
    % Compute the bra (row) vector for
    % the GHZ state, given the number n of parties
    ghz_state_ket_vector = ...
        compute_ghz_state_ket_vector(n);

    % Compute the bra (row) vector for
    % the GHZ state, given the number n of parties 
    ghz_state_bra_vector = ...
        ghz_state_ket_vector';
    
end


% Define a function to compute the density matrix
% for the GHZ state, given a number n of parties
function ghz_state_density_matrix = ...
         compute_ghz_state_density_matrix(n)
    
    % Compute the ket (column) vector for
    % the GHZ state, given the number n of parties 
    ghz_state_ket_vector = ...
        compute_ghz_state_ket_vector(n);

    % Compute the bra (row) vector for
    % the GHZ state, given the number n of parties 
    ghz_state_bra_vector = ...
        compute_ghz_state_bra_vector(n);
    
    % Compute the density matrix for
    % the GHZ state, given a number n of parties
    ghz_state_density_matrix = ...
        ghz_state_ket_vector * ghz_state_bra_vector;

end


% Define a function to compute the density matrix
% for a probabilistic noisy GHZ state, given
% a number n of parties and a probability v
function probabilistic_ghz_state_density_matrix = ...
         compute_probabilistic_ghz_state_density_matrix...
            (n, probability_v)
    
    % Compute the density matrix
    % for a probabilistic noisy GHZ state, given
    % the number n of parties and a probability v
    probabilistic_ghz_state_density_matrix = ...
        probability_v * compute_ghz_state_density_matrix(n);

end


% Define a function to compute the white noise density matrix
% given a number n of parties and a probability v
function probabilistic_white_noise_density_matrix = ...
    compute_probabilistic_white_noise_density_matrix(n, probability_v)
    
    % Compute the white noise density matrix
    % given the number n of parties and a probability v
    probabilistic_white_noise_density_matrix = ...
        (1 - probability_v) * (1 / 2^n) * eye(2^n);

end


% Define a function to compute the (final)
% density matrix for a probabilistic noisy GHZ state,
% given a number n of parties and a probability v
function probabilistic_werner_noisy_ghz_state_density_matrix = ...
    compute_probabilistic_werner_noisy_ghz_state_density_matrix...
        (n, probabilistic_v)
    
    % Compute the (final) density matrix for
    % a probabilistic noisy GHZ state, given
    % the number n of parties and a probability v
    probabilistic_werner_noisy_ghz_state_density_matrix = ...
        compute_probabilistic_ghz_state_density_matrix...
            (n, probabilistic_v) + ...
                compute_probabilistic_white_noise_density_matrix...
                    (n, probabilistic_v);

end


% Define a function to compute the quantum bound Q_n
% for the expectation value <M_n>, given a number n of parties
function quantum_bound_expectation_value = ...
    compute_quantum_bound_expectation_value(n)
    
    % Compute the quantum bound Q_n for
    % the expectation value <M_n>, given a number n of parties
    quantum_bound_expectation_value = ...
        2^( (n - 1) / 2 );

end


% Define a function to compute the A_(2) global operator
function a_2_global_operator = ...
    compute_a_2_global_operator(a_operator)
    
    % Compute the A_(2) global operator
    a_2_global_operator = kron(a_operator, a_operator);
    
end


% Define a function to compute the M_(n-1) Bell operator,
% given a number n of parties, as well as, A_0 and A_1 operators
function m_n_minus_1_1_bell_operator = ...
    compute_m_n_minus_1_1_bell_operator(n, a_operator_1, a_operator_2)
    
    % Compute the (initial) M_(n-1) Bell operator,
    % from the A_(2) global operator (basis case)
    m_n_minus_1_1_bell_operator = ...
        compute_a_2_global_operator(a_operator_1);

    
    % If the number n of parties is greater than 3
    if( n > 3 )
        
        % Initialize the counter for consecutive Bell operators
        num_consecutive_bell_operators = 0;
        
        % Initialize the signal for the Bell operator pair
        bell_operator_pair_signal = -1;
    
        
        % For the remaining number n of parties,
        % within the range [3, N-1]
        for i = 3:n-1
            
           % If were applied two consecutive Bell operators
           if( num_consecutive_bell_operators == 2 )
               
               % Flip the signal for the (next) Bell operator pair
               bell_operator_pair_signal = ...
                   -1 * bell_operator_pair_signal;
                
               % Reset the counter for consecutive Bell operators
               num_consecutive_bell_operators = 0;
    
           end
    
            
           % If were not applied two consecutive Bell operators
           if( num_consecutive_bell_operators ~= 2 )
                
               % Compute the M_(n-1) Bell operator
               m_n_minus_1_1_bell_operator = ...
                   kron( m_n_minus_1_1_bell_operator, ...
                         ( a_operator_1 + ...
                           bell_operator_pair_signal * a_operator_2 ) );
               
               % Increase the counter for consecutive Bell operators
               num_consecutive_bell_operators = ...
                   num_consecutive_bell_operators + 1;
    
           end
    
        end

    end

    
    % Normalize the M_(n-1) Bell operator with
    % the corresponding norm factor
    m_n_minus_1_1_bell_operator = ...
        ( 1 / sqrt( 2^(n - 2) ) ) * m_n_minus_1_1_bell_operator;

    % Compute the quantum bound Q_(n-1) for
    % the expectation value <M_(n-1)>
    quantum_bound_expectation_value = ...
        compute_quantum_bound_expectation_value(n - 1);
    
    % Multiply the quantum bound Q_(n-1)
    % for the expectation value <M_(n-1)>,
    % regarding the corresponding M_(n-1) Bell operator
    m_n_minus_1_1_bell_operator = ...
        quantum_bound_expectation_value * m_n_minus_1_1_bell_operator;

end


% Define a function to compute the M'_(n-1) Bell operator,
% given a number n of parties, as well as, A_0 and A_1 operators
function m_n_minus_1_2_bell_operator = ...
    compute_m_n_minus_1_2_bell_operator(n, a_operator_1, a_operator_2)
    
    % Compute the (initial) M'_(n-1) Bell operator,
    % from the A_(2) global operator (basis case)
    m_n_minus_1_2_bell_operator = ...
        compute_a_2_global_operator(a_operator_2);


    % If the number n of parties is greater than 3
    if( n > 3 )
        

        % Initialize the counter for consecutive Bell operators
        num_consecutive_bell_operators = 0;
        
        % Initialize the signal for the Bell operator pair
        bell_operator_pair_signal = 1;

        

        % For the remaining number n of parties,
        % within the range [3, N-1]
        for i = 3:n-1
           
           % If were applied two consecutive Bell operators
           if( num_consecutive_bell_operators == 2 )
               
               % Flip the signal for the (next) Bell operator pair
               bell_operator_pair_signal = ...
                   -1 * bell_operator_pair_signal;
                
               % Reset the counter for consecutive Bell operators
               num_consecutive_bell_operators = 0;
    
           end
    

           % If were not applied two consecutive Bell operators
           if( num_consecutive_bell_operators ~= 2 )
               
               % Compute the M'_(n-1) Bell operator
               m_n_minus_1_2_bell_operator = ...
                   kron( m_n_minus_1_2_bell_operator, ...
                         ( a_operator_1 + ...
                           bell_operator_pair_signal * a_operator_2 ) );
               
               % Increase the counter for consecutive Bell operators
               num_consecutive_bell_operators = ...
                   num_consecutive_bell_operators + 1;
    
           end
    
        end

    end


    % Normalize the M_(n-1) Bell operator with
    % the corresponding norm factor
    m_n_minus_1_2_bell_operator = ...
        ( 1 / sqrt( 2^(n - 2) ) ) * m_n_minus_1_2_bell_operator;

    % Compute the quantum bound Q'_(n-1) for
    % the expectation value <M'_(n-1)>
    quantum_bound_expectation_value = ...
        compute_quantum_bound_expectation_value(n - 1);
    
    % Multiply the quantum bound Q'_(n-1)
    % for the expectation value <M'_(n-1)>,
    % regarding the corresponding M'_(n-1) Bell operator
    m_n_minus_1_2_bell_operator = ...
        quantum_bound_expectation_value * m_n_minus_1_2_bell_operator;
    
end


% Define the function to compute the M_(n) Bell operator,
% given a number n of parties, as well as, A_0 and A_1 operators
function m_n_bell_operator = ...
    compute_m_n_bell_operator(n, a_operator_1, a_operator_2)
    
    % If the number n of parties is greater or equal to 3
    if( n >= 3 )
        
        % Compute the M_(n-1) Bell operator,
        % given the number n of parties, as well as,
        % A_0 and A_1 operators
        m_n_minus_1_1_bell_operator = ...
            compute_m_n_minus_1_1_bell_operator...
                (n, a_operator_1, a_operator_2);
        
        % Compute the M'_(n-1) Bell operator,
        % given the number n of parties, as well as,
        % A_0 and A_1 operators
        m_n_minus_1_2_bell_operator = ...
            compute_m_n_minus_1_2_bell_operator...
                (n, a_operator_1, a_operator_2);
        
        % Compute the first term of
        % the sum for the M_(n) Bell operator,
        % given the number n of parties, as well as,
        % A_0 and A_1 operators
        m_n_bell_operator_term_1 = ...
            kron(m_n_minus_1_1_bell_operator, ...
                 a_operator_1 + a_operator_2);

        % Compute the second term of
        % the sum for the M_(n) Bell operator,
        % given the number n of parties, as well as,
        % A_0 and A_1 operators
        m_n_bell_operator_term_2 = ...
            kron(m_n_minus_1_2_bell_operator, ...
                 a_operator_1 - a_operator_2);
        
        % Compute the M_(n) Bell operator,
        % from the sum of the respective terms
        m_n_bell_operator = ...
            m_n_bell_operator_term_1 + m_n_bell_operator_term_2;

    end

end


% Define the function to compute the expectation value <M_(n)>
% for the GHZ state for n parties and the M_(n) Bell operator,
% given a number n of parties, as well as, A_0 and A_1 operators
function m_n_werner_noisy_ghz_state_expectation_value = ...
    compute_m_n_werner_noisy_ghz_state_expectation_value...
        (n, probability_v, a_operator_1, a_operator_2)
    
    % Compute the density matrix for
    % the GHZ state, given the number n of parties
    probabilistic_werner_noisy_ghz_state_density_matrix = ...
        compute_probabilistic_werner_noisy_ghz_state_density_matrix...
            (n, probability_v);
    
    % Compute the M_(n) Bell operator,
    % given the number n of parties, as well as, A_0 and A_1 operators
    m_n_bell_operator = ...
        compute_m_n_bell_operator(n, a_operator_1, a_operator_2);
    
    % Conpute the expectation value <M_(n)> from
    % the trace for the multiplication of the density matrix for
    % the GHZ state for n parties and the  M_(n) Bell operator
    m_n_werner_noisy_ghz_state_expectation_value = ...
        trace(probabilistic_werner_noisy_ghz_state_density_matrix * ...
              m_n_bell_operator);
    
end