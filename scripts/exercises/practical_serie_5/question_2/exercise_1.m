% Clear Command Window
clc;


% Define the size of the Hilbert Space
hilbert_space_size = 2;


% Define the 1st list of projectors to test
% Note(s):
% * This list of projectors do not
%   sum the Identity Matrix I_{2x2} on
%   purpose to fail the test conditions,
%   resulting on no coloring map v;
projectors_list_1 = { [  0.5  0.5 ; 
                         0.5  0.5 ];
                      [  0.5  0.5 ;
                        -0.5  0.5 ] };

% Compute the number of projectors in
% the 1st list of projectors to test
num_projectors_1 = ...
    length(projectors_list_1);

% Compute the list of unitary vectors from
% the 1st list of projectors to test
unitary_vectors_1 = ...
    projectors_to_unitary_vectors(projectors_list_1);

% Compute the number of unitary vectors
% from the 1st list of projectors to test
num_unitary_vectors_1 = ...
    length(unitary_vectors_1);


% Define the 2nd list of projectors to test
% Note(s):
% * This list of projectors corresponds
%   to the bases |+> and |->
projectors_list_2 = { [  0.5  0.5 ; 
                         0.5  0.5 ];
                      [  0.5 -0.5 ;
                        -0.5  0.5 ] };

% Compute the number of projectors in
% the 2nd list of projectors to test
num_projectors_2 = ...
    length(projectors_list_2);

% Compute the list of unitary vectors
% from the 2nd list of projectors to test
unitary_vectors_2 = ...
    projectors_to_unitary_vectors(projectors_list_2);

% Compute the number of unitary vectors
% from the 2nd list of projectors to test
num_unitary_vectors_2 = ...
    length(unitary_vectors_2);


% Define the 3rd list of projectors to test
% Note(s):
% * This list of projectors corresponds
%   to the bases |0> and |1>
projectors_list_3 = { [ 1 0 ; 
                        0 0 ];
                      [ 0 0 ;
                        0 1 ] };

% Compute the number of projectors in
% the 3rd list of projectors to test
num_projectors_3 = ...
    length(projectors_list_3);

% Compute the list of unitary vectors
% from the 3rd list of projectors to test
unitary_vectors_3 = ...
    projectors_to_unitary_vectors(projectors_list_3);

% Compute the number of unitary vectors
% from the 3rd list of projectors to test
num_unitary_vectors_3 = ...
    length(unitary_vectors_3);


% Print a blank line
fprintf("\n");


% Analyse if the 1st list of projectors and their respective
% unitary vectors represents a possible Kochen-Specker (KS) set
analyse_possible_kochen_specker_set...
    (projectors_list_1, num_projectors_1, ...
     unitary_vectors_1, num_unitary_vectors_1, ...
     hilbert_space_size);


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("****************************" + ...
        "****************************\n");

% Print two blank lines
fprintf("\n\n");


% Analyse if the 2nd list of projectors and their respective
% unitary vectors represents a possible Kochen-Specker (KS) set
analyse_possible_kochen_specker_set...
    (projectors_list_2, num_projectors_2, ...
     unitary_vectors_2, num_unitary_vectors_2, ...
     hilbert_space_size)


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("****************************" + ...
        "****************************\n");

% Print two blank lines
fprintf("\n\n");


% Analyse if the 3rd list of projectors and their respective
% unitary vectors represents a possible Kochen-Specker (KS) set
analyse_possible_kochen_specker_set...
    (projectors_list_3, num_projectors_3, ...
     unitary_vectors_3, num_unitary_vectors_3, ...
     hilbert_space_size)


% Print a blank line
fprintf("\n");


% Define a function to check if
% a given matrix M is a projector
function is_matrix_projector = ...
         check_matrix_projector(matrix)
    
    % Check if the given matrix M is a projector,
    % by verifying if M = M^2 = M x M and M = M^(dagger)
    is_matrix_projector = ...
        isequal(matrix, matrix * matrix) && ...
            isequal(matrix, matrix') && ...
                ismembertol(trace(matrix), 1, 1e-6);

end


% Define a function to check if
% a given matrix M has rank equal to 1
function is_matrix_rank_1 = ...
         check_matrix_rank_1(matrix)
    
    % Check if the rank of the given
    % matrix M is equal to one
    is_matrix_rank_1 = ...
        rank(matrix) == 1;

end


% Define a function to compute
% the Kronecker Delta value,
% for two given values
function kronecker_delta = ...
         compute_kronecker_delta(i, j)
    
    % Compute the Kronecker Delta value,
    % for the two given values
    kronecker_delta = (i == j);

end


% Define a function to check if
% two given projectors are orthogonal
function are_projectors_orthogonal = ...
         check_projectors_orthogonality(projector_1, projector_2)
    
    % Check if the two given projectors are orthogonal
    are_projectors_orthogonal = ...
        ismembertol(trace( ( projector_1 * projector_2 ) ), 0, 1e-6);

end


% Define a function to compute the list of
% the respective associated unitary vectors
% from a given list of projectors
function unitary_vectors = ...
    projectors_to_unitary_vectors(projectors_list)

    % Compute the number of the projectors,
    % from which will be computed the list of
    % the respective associated unitary vectors
    num_projectors = length(projectors_list);

    % Create a list for the unitary vectors
    % to be computed from the projectors
    unitary_vectors = cell(num_projectors, 1);

    
    % For the index of each one of the projectors,
    % from which will be computed the list of
    % the respective associated unitary vectors
    for curr_projector_idx = 1:num_projectors
        
        % Retrieve the current projector
        % from the respective associated list
        projector = projectors_list{curr_projector_idx};

        % Compute the eigen decomposition of
        % the current projector retrieved,
        % extracting its respective eigenvectos
        % and the diagonal matrix with eigenvalues
        [ eigenvectors, ...
          eigenvalues_diagonal_matrix ] = eig(projector);
       
        % Compute the index of the maximum eigenvalue
        % from the diagonal vector which is obtained on
        % the diagonal matrix with eigenvalues
        [~, max_eigenvalue_index] = ...
            max(diag(eigenvalues_diagonal_matrix));
        
        % Compute the unitary vector from
        % the eigenvector corresponding to
        % the index of the maximum eigenvalue
        unitary_vectors{curr_projector_idx} = ...
            eigenvectors(:, max_eigenvalue_index);
        
    end

end


% Define a function to find the coloring set
% for a list of projectors, according to
% a certain size of the Hilbert Space
function [coloring_set, conditions_met_flag] = ...
         find_coloring_set(projectors_list, hilbert_space_size)
    
    % Define a boolean flag for
    % the necessary conditions to be met
    % in order to consider a set of projectors
    % to be a possible candidate to be colored
    conditions_met_flag = true;


    % Compute the number of the projectors,
    % from which will be computed the list of
    % the respective associated unitary vectors
    num_projectors = length(projectors_list);
    
    % Create a matrix for the sum of the projectors
    % retrieved from the respective associated list
    projectors_sum = ...
        zeros(hilbert_space_size, hilbert_space_size);
    
    
    % For the index of each one of the projectors,
    % from which will be computed the list of
    % the respective associated unitary vectors
    for curr_projector_idx = 1:num_projectors

        % Retrieve the current projector
        % from the respective associated list
        projector = projectors_list{curr_projector_idx};
        
        
        % If the matrix corresponding to
        % the current projector is not really
        % a projector, neither has a rank equal to one
        if( ( ~check_matrix_projector(projector) ) || ...
            ( ~check_matrix_rank_1(projector) ) )
            
            % Set the coloring map as an empty set
            coloring_set = [];

            % In this case, the necessary conditions to
            % be met in order to consider a set of projectors
            % to be a possible candidate to be colored,
            % were not all fulfilled
            conditions_met_flag = false;
            
            % Return out of the function
            return;
        
        % If the matrix corresponding to
        % the current projector is really
        % a projector and has a rank equal to one
        else
            
            % Sum the current projector to
            % its corresponding sum of projectors
            projectors_sum = ...
                projectors_sum + projector;

        end
   
    end


    % For the index of each one of the projectors,
    % in the respective associated 1st list,
    % from which will be computed the list of
    % the respective associated unitary vectors
    for curr_projector_1_idx = 1:num_projectors

        % Retrieve the current 1st projector,
        % from the corresponding list of projectors
        projector_1 = projectors_list{curr_projector_1_idx};
    

        % For the index of each one of the projectors,
        % in the respective associated 2nd list,
        % from which will be computed the list of
        % the respective associated unitary vectors
        for curr_projector_2_idx = 1:num_projectors
            
            % Retrieve the current 2nd projector,
            % from the corresponding list of projectors
            projector_2 = projectors_list{curr_projector_2_idx};
            
            % Compute the matrix multiplication between
            % the current 1st and 2nd retrieved projectors
            projectors_1_2_mult = projector_1 * projector_2;
            
            % Compute the Kronecker Delta value
            % for the indices of the current
            % 1st and 2nd retrieved projectors
            kronecker_delta_1_2 = ...
                compute_kronecker_delta(curr_projector_1_idx, ...
                                        curr_projector_2_idx);
            
            % Compute the multiplication
            % between the Kronecker Delta value
            % for the indices of the current
            % 1st and 2nd retrieved projectors,
            % and the 1st retrieved projector
            kronecker_delta_1_2_projector_1_mult = ...
                kronecker_delta_1_2 * projector_1;
              

            % If the matrix multiplication between
            % the current 1st and 2nd retrieved projectors
            % is not equal to the multiplication
            % between the Kronecker Delta value
            % for the indices of the current
            % 1st and 2nd retrieved projectors,
            % and the 1st retrieved projector
            if( ~ismembertol(projectors_1_2_mult, ...
                             kronecker_delta_1_2_projector_1_mult, ...
                             1e-6, 'ByRows', true) )
                
                % Set the coloring map as an empty set
                coloring_set = [];
                
                % In this case, the necessary conditions to
                % be met in order to consider a set of projectors
                % to be a possible candidate to be colored,
                % were not all fulfilled
                conditions_met_flag = false;

                % Return out of the function
                return;
        
            end   
    
        end

    end

    % If the sum of projectors is not
    % equal to an identity matrix I_{2x2}
    if( ~ismembertol(projectors_sum, ...
                     eye(hilbert_space_size), ...
                     1e-6, 'ByRows', true) )
        
        % Set the coloring map as an empty set
        coloring_set = [];

        % In this case, the necessary conditions to
        % be met in order to consider a set of projectors
        % to be a possible candidate to be colored,
        % were not all fulfilled
        conditions_met_flag = false;
        
        % Return out of the function
        return;

    end    

    
    % Initialize the coloring set with
    % an uncolored value, according to
    % the number of projectors
    coloring_set = -1 * ones(1, num_projectors);

    
    % For the index of each one of the projectors,
    % in the respective associated 1st list,
    % from which will be computed the list of
    % the respective associated unitary vectors
    for curr_projector_1_idx = 1:num_projectors

        % If the color currently associated to
        % the current 1st projector is still uncolored
        if( coloring_set(curr_projector_1_idx) == -1 )
            
            % Retrieve the current 1st projector,
            % from the corresponding list of projectors
            projector_1 = projectors_list{curr_projector_1_idx};
            
            % Define the color currently associated to
            % the current 1st projector is the 1st color
            coloring_set(curr_projector_1_idx) = 0;
    
            
            % For the index of each one of the projectors,
            % in the respective associated 2nd list,
            % from which will be computed the list of
            % the respective associated unitary vectors
            for curr_projector_2_idx = 1:num_projectors
                
                % Retrieve the current 2nd projector,
                % from the corresponding list of projectors
                projector_2 = projectors_list{curr_projector_2_idx};
            
                
                % If the indexes of the two current projectors
                % are different and the corresponding two projectors
                % are orthogonal between themselves
                if( ( curr_projector_1_idx ~= curr_projector_2_idx ) && ...
                    ( check_projectors_orthogonality(projector_1, projector_2) ) )
                    
                    % If the color currently associated to
                    % the current 2nd projector is still uncolored
                    if( coloring_set(curr_projector_2_idx) == -1 )
                        
                        % Define the color currently associated to
                        % the current 2nd projector is the 2nd color
                        coloring_set(curr_projector_2_idx) = 1;
                    
                    % If the color currently associated to
                    % the current 2nd projector is the 1st color
                    % Note(s):
                    % * In this case, for any unitary vector
                    %   corresponding to a projector mapped with
                    %   the 1st color, all the other unitary vectors
                    %   being orthogonal to the former, needs to be
                    %   mapped to the 2nd color (and never mapped to
                    %   the 1st color), i.e., two orthogonal unitary
                    %   vectors can never both be mapped to the 1st color
                    elseif( coloring_set(curr_projector_2_idx) == 0 )
                        
                        % Set the coloring map as an empty set
                        coloring_set = [];
                        
                        % Return out of the function
                        return;

                    end

                end

            end

        end

    end
    
    
    % If the coloring map is not an empty set
    if( ~isempty(coloring_set) )

        % If the sum of all the mapped colors in
        % the coloring map does not sum to 1
        if( sum(coloring_set) ~= 1 )
            
            % Set the coloring map as an empty set
            coloring_set = [];
        
        end

    end
    
end


% Define a function to check if the unitary vectors
% corresponding to a list of projectors represent
% a Kochen Specker (KS) set, verifying the coloring map
% for those same unitary vectors obtained from the projectors
function [is_kochen_specker_set, conditions_met_flag, coloring_set] = ...
         check_unitary_vectors_from_projectors_are_kochen_specker_set...
         (projectors_list, hilbert_space_size)
    
    % Find the coloring set for the given
    % list of projectors, according to the given
    % corresponding size of the Hilbert Space
    [coloring_set, conditions_met_flag] = ...
        find_coloring_set(projectors_list, hilbert_space_size);
    
    
    % If the coloring map is an empty set,
    % meeting all the required conditions
    % and previously defined
    if( isempty(coloring_set) && conditions_met_flag )
        
        % In this case, when the coloring map
        % is an empty set, the unitary vectors
        % represent a Kochen-Specker (KS) set
        is_kochen_specker_set = true;
    
    % If the coloring map is not an empty set,
    % meeting all the required conditions
    % and previously defined
    elseif( ~isempty(coloring_set) && conditions_met_flag )
        
        % In this case, when the coloring map
        % is not an empty set, the unitary vectors
        % does not represent a Kochen-Specker (KS) set
        is_kochen_specker_set = false;
    
    % If the coloring map is an empty set,
    % not even meeting all the required conditions
    % and previously defined
    elseif( isempty(coloring_set) && ~conditions_met_flag )
    
        % In this case, when the coloring map
        % is an empty set but the respective
        % projectors do not meet all the required
        % conditions, the unitary vectors does not
        % represent a Kochen-Specker (KS) set
        is_kochen_specker_set = false;

    end

end


% Define a function to analyse if
% a given list of projectors, as well as
% their respective unitary vectors and
% size of the Hilbert Space
function analyse_possible_kochen_specker_set...
    (projectors_list, num_projectors, ...
     unitary_vectors, num_unitary_vectors, ...
     hilbert_space_size)
    
    % Check if the unitary vectors
    % corresponding to the given list of projectors
    % represent a Kochen Specker (KS) set, verifying
    % the coloring map for those same unitary vectors
    % obtained from the respective projectors and
    % size of the Hilbert Space
    [is_kochen_specker_set, conditions_met_flag, coloring_set] = ...
        check_unitary_vectors_from_projectors_are_kochen_specker_set...
            (projectors_list, hilbert_space_size);
    
    % If the necessary conditions to
    % be met in order to consider a set of projectors
    % to be a possible candidate to be colored,
    % were not all fulfilled
    if(~conditions_met_flag)
        
        % Print the headline information about
        % the unitary vectors corresponding to
        % the given list of projectors representing
        % a Kochen Specker (KS) set (for n = 2)
        fprintf("The following set <strong>U</strong> of " + ...
                "unitary vectors\n" + ...
                "<strong>does not</strong> even meet the required\n" + ...
                "<strong>conditions</strong> to " + ...
                "be a <strong>candidate</strong> set to \n" + ...
                "a <strong>Kochen-Specker (KS)</strong> set " + ...
                "for n = %d:\n\n", ...
                hilbert_space_size);
    
        
        % For the index of each one of the unitary vectors,
        % corresponding to the respective list of projectors
        for curr_unitary_vector_idx = 1:num_unitary_vectors
            
            % Retrieve the current unitary vector,
            % from the corresponding list of unitary vectors
            unitary_vector = ...
                unitary_vectors(curr_unitary_vector_idx);
            
            % Print the introductory notation
            % for the current unitary vector
            fprintf(" * <strong>u_%d</strong> = [ ", ...
                    ( curr_unitary_vector_idx - 1 ));
            
            
            % For the index of each one of the complex
            % coefficients of the current unitary vector
            for complex_coeff_idx = 1:hilbert_space_size
                
                % Retrieve the current complex
                % coefficient of the current unitary vector 
                complex_coeff = unitary_vector{1}(complex_coeff_idx);
                            
                % Print the current complex
                % coefficient of the current unitary vector
                fprintf("%.4f", complex_coeff);
    
                
                % If the current complex coefficient of
                % the current unitary vector is not the last one
                if( complex_coeff_idx < hilbert_space_size )
                    
                    % Print a comma separator for
                    % the next complex coefficient
                    fprintf(", ");
                
                % If the current complex coefficient of
                % the current unitary vector is the last one
                else
                    
                    % Print a closing bracket, once all
                    % the complex coefficients are already printed
                    fprintf(" ];\n");
    
                end
        
            end
        
        end
        
        
        % Print two blank lines
        fprintf("\n\n");
    
        
        % Print the headline information
        % about the projectors corresponding to
        % the given list of projectors representing
        % a Kochen Specker (KS) set (for n = 2)
        fprintf("For the following set " + ...
                "<strong>P</strong> of projectors:\n\n");
        
        
        % For the index of each one of the projectors
        for curr_projector_idx = 1:num_projectors
            
            % Retrieve the current projector
            projector = ...
                projectors_list{curr_projector_idx};
            
            % Print the introductory notation
            % for the current projector
            fprintf(" * <strong>P_%d</strong> = \n", ...
                    ( curr_projector_idx - 1 ));
            
            % Display the content of the current projector
            disp(projector);
            
        end
        

        % Print a blank line
        fprintf("\n");
    
        % Print the information about no possible coloring map
        % candidates for the given list of projectors
        % and their corresponding unitary vectors
        fprintf("<strong>No</strong> possible coloring map " + ...
                "<strong>v</strong> candidates!\n");
        
        % Return out of the function
        return;

    end

    
    % If the unitary vectors corresponding to
    % the given list of projectors represent
    % a Kochen Specker (KS) set
    if(is_kochen_specker_set)
        
        % Print the headline information about
        % the unitary vectors corresponding to
        % the given list of projectors representing
        % a Kochen Specker (KS) set (for n = 2)
        fprintf("The following set <strong>U</strong> of " + ...
                "unitary vectors is \n" + ...
                "a <strong>Kochen-Specker (KS)</strong> set " + ...
                "for n = %d:\n\n", ...
                hilbert_space_size);
    
        
        % For the index of each one of the unitary vectors,
        % corresponding to the respective list of projectors
        for curr_unitary_vector_idx = 1:num_unitary_vectors
            
            % Retrieve the current unitary vector,
            % from the corresponding list of unitary vectors
            unitary_vector = ...
                unitary_vectors(curr_unitary_vector_idx);
            
            % Print the introductory notation
            % for the current unitary vector
            fprintf(" * <strong>u_%d</strong> = [ ", ...
                    ( curr_unitary_vector_idx - 1 ));
            
            
            % For the index of each one of the complex
            % coefficients of the current unitary vector
            for complex_coeff_idx = 1:hilbert_space_size
                
                % Retrieve the current complex
                % coefficient of the current unitary vector 
                complex_coeff = unitary_vector{1}(complex_coeff_idx);
                            
                % Print the current complex
                % coefficient of the current unitary vector
                fprintf("%.4f", complex_coeff);
    
                
                % If the current complex coefficient of
                % the current unitary vector is not the last one
                if( complex_coeff_idx < hilbert_space_size )
                    
                    % Print a comma separator for
                    % the next complex coefficient
                    fprintf(", ");
                
                % If the current complex coefficient of
                % the current unitary vector is the last one
                else
                    
                    % Print a closing bracket, once all
                    % the complex coefficients are already printed
                    fprintf(" ];\n");
    
                end
        
            end
        
        end
        
        
        % Print two blank lines
        fprintf("\n\n");
    
        
        % Print the headline information
        % about the projectors corresponding to
        % the given list of projectors representing
        % a Kochen Specker (KS) set (for n = 2)
        fprintf("For the following set " + ...
                "<strong>P</strong> of projectors:\n\n");
        
        
        % For the index of each one of the projectors
        for curr_projector_idx = 1:num_projectors
            
            % Retrieve the current projector
            projector = ...
                projectors_list{curr_projector_idx};
            
            % Print the introductory notation
            % for the current projector
            fprintf(" * <strong>P_%d</strong> = \n", ...
                    ( curr_projector_idx - 1 ));
            
            % Display the content of the current projector
            disp(projector);
            
        end
        

        % Print a blank line
        fprintf("\n");
    
        % Print the information about no coloring map
        % being found for the given list of projectors
        % and their corresponding unitary vectors
        fprintf("<strong>No</strong> coloring map " + ...
                "<strong>v</strong> was found!\n");
    
    
    % If the unitary vectors corresponding to
    % the given list of projectors do not
    % represent a Kochen Specker (KS) set
    else
        
        % Print the headline information about
        % the unitary vectors corresponding to
        % the given list of projectors not representing
        % a Kochen Specker (KS) set (for n = 2)
        fprintf("The following set <strong>U</strong> of " + ...
                "unitary vectors is\n" + ...
                "<strong>not</strong> a <strong>Kochen-Specker</strong> " + ...
                "<strong>(KS)</strong> " + ...
                "set for <strong>n = %d</strong>:\n\n", ...
                hilbert_space_size);
    
        
        % For the index of each one of the unitary vectors,
        % corresponding to the respective list of projectors
        for curr_unitary_vector_idx = 1:num_unitary_vectors
            
            % Retrieve the current unitary vector,
            % from the corresponding list of unitary vectors
            unitary_vector = ...
                unitary_vectors(curr_unitary_vector_idx);
    
            % Print the introductory notation
            % for the current unitary vector
            fprintf(" * <strong>u_%d</strong> = [ ", ...
                    ( curr_unitary_vector_idx - 1 ));
            
    
            % For the index of each one of the complex
            % coefficients of the current unitary vector
            for complex_coeff_idx = 1:hilbert_space_size
                
                % Retrieve the current complex
                % coefficient of the current unitary vector
                complex_coeff = unitary_vector{1}(complex_coeff_idx);
                            
                % Print the current complex
                % coefficient of the current unitary vector
                fprintf("%.4f", complex_coeff);
    
    
                % If the current complex coefficient of
                % the current unitary vector is not the last one
                if( complex_coeff_idx < hilbert_space_size )
                    
                    % Print a comma separator for
                    % the next complex coefficient
                    fprintf(", ");
    
                % If the current complex coefficient of
                % the current unitary vector is the last one
                else
                    
                    % Print a closing bracket, once all
                    % the complex coefficients are already printed
                    fprintf(" ];\n");
    
                end
        
            end
        
        end
        
        
        % Print two blank lines
        fprintf("\n\n");
    
    
        % Print the headline information
        % about the projectors corresponding to
        % the given list of projectors not representing
        % a Kochen Specker (KS) set (for n = 2)
        fprintf("For the following set " + ...
                "<strong>P</strong> of projectors:\n\n");
        

        % For the index of each one of the projectors
        for curr_projector_idx = 1:num_projectors
            
            % Retrieve the current projector
            projector = ...
                projectors_list{curr_projector_idx};
    
            % Print the introductory notation
            % for the current projector
            fprintf(" * <strong>P_%d</strong> = \n", ...
                    ( curr_projector_idx - 1 ));
            
            % Display the content of the current projector
            disp(projector);
            
        end
        
        
        % Print a blank line
        fprintf("\n");
    
        
        % Print the introductory information about
        % the coloring map found for the given list of
        % projectors and their corresponding unitary vectors
        fprintf("The coloring map <strong>v</strong> found is:\n");
        fprintf(" * <strong>v</strong> = [ ");
            
        
        % For the index of each color whihc was assigned to
        % each projector/unitary vector in the coloring map found
        for color_idx = 1:num_projectors
            
            % Retrieve the current color assignment
            color = coloring_set(color_idx);
                       
            % Print the information about
            % the current color assignment
            fprintf("%d", color);
    
            
            % If the current color assigment for
            % the current projector/unitary vector
            % is not the last one
            if( color_idx < num_projectors )
                
                % Print a comma separator for
                % the next color assigment
                fprintf(", ");
    
            % If the current color assigment for
            % the current projector/unitary vector
            % is the last one
            else
                
                % Print a closing bracket, once all
                % the color assignments are already printed
                fprintf(" ];\n");
    
            end
    
        end
    
    end

end