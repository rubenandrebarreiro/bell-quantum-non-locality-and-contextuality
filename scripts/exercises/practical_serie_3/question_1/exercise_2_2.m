% Clear Command Window
clc;


% Define the initial total number of
% observers for the Mermin Inequality
num_initial_total_observers = 2;

% Define the final total number of
% observers for the Mermin Inequality
num_final_total_observers = 10;


% Print a blank line
fprintf('\n');


parties_measurements = ...
    generate_parties_measurements(num_initial_total_observers);

rho_density_matrix_for_epr_pair_or_ghz_state = ...
     generate_rho_density_matrix_for_epr_pair_or_ghz_state...
     (num_initial_total_observers);




function ket_epr_pair_phi_plus = ...
         generate_ket_epr_pair_phi_plus()
    
    num_total_observers = 2;

    ket_epr_pair_phi_plus = zeros(1, 2^num_total_observers);
    
    ket_epr_pair_phi_plus(1) = ( 1 / sqrt(2) );

    ket_epr_pair_phi_plus(2^num_total_observers) = ( 1 / sqrt(2) );
    
    ket_epr_pair_phi_plus = ket_epr_pair_phi_plus.';

end

function ket_ghz_state = ...
         generate_ket_ghz_state(num_total_observers)
    
    ket_ghz_state = zeros(1, 2^num_total_observers);
    
    ket_ghz_state(1) = ( 1 / sqrt(2) );

    ket_ghz_state(2^num_total_observers) = ( 1 / sqrt(2) );
    
    ket_ghz_state = ket_ghz_state.';

end

function ket_epr_pair_or_ghz_state = ...
         generate_ket_epr_pair_or_ghz_state(num_total_observers)

    if num_total_observers == 2
    
        ket_epr_pair_or_ghz_state = ...
            generate_ket_epr_pair_phi_plus();

    end


    if num_total_observers > 2
    
        ket_epr_pair_or_ghz_state = ...
            generate_ket_ghz_state(num_total_observers);
        
    end

end

function pauli_sigmas_vec = ...
         generate_pauli_sigmas_vec()
    
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
    
end

function rho_density_matrix_for_epr_pair_or_ghz_state = ...
         generate_rho_density_matrix_for_epr_pair_or_ghz_state...
         (num_total_observers)

    ket_epr_pair_or_ghz_state = ...
        generate_ket_epr_pair_or_ghz_state(num_total_observers);
    
    bra_epr_pair_or_ghz_state = ...
        ket_epr_pair_or_ghz_state.';

    % Compute the density
    % matrix rho = |psi><psi| for
    % the two-dimensional quantum state psi
    rho_density_matrix_for_epr_pair_or_ghz_state = ...
        ket_epr_pair_or_ghz_state * bra_epr_pair_or_ghz_state;
    
end

function parties_measurements = ...
         generate_parties_measurements(num_total_observers)    
    
    pauli_sigmas_vec = ...
        generate_pauli_sigmas_vec();

    pauli_sigma_x = ...
        pauli_sigmas_vec{1};

    pauli_sigma_z = ...
        pauli_sigmas_vec{3};

    
    odd_party_measurements = { pauli_sigma_x, ...
                               pauli_sigma_z };
        
    even_party_measurements = { ( 1 / sqrt(2) ) * ( pauli_sigma_x + ...
                                                    pauli_sigma_z ), ...
                                ( 1 / sqrt(2) ) * ( pauli_sigma_x - ...
                                                    pauli_sigma_z ) };
    

    parties_measurements = cell(num_total_observers);
    

    for current_num_parties = 1:num_total_observers

        if mod(current_num_parties, 2) ~= 0
            
            parties_measurements{current_num_parties} = ...
                odd_party_measurements;
        
        end


        if mod(current_num_parties, 2) == 0
            
            parties_measurements{current_num_parties} = ...
                even_party_measurements;
            
        end
        
    end

end