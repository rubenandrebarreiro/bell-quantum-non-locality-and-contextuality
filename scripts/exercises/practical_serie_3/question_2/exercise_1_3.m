% Clear Command Window
clc;


% Define the number of rows for the (binary)
% predicate table for the CHSH XOR Game,
% which will be equal to the total
% number of strategies
num_rows = 2^4;

% Define the number of columns for the (binary)
% predicate table for the CHSH XOR Game
num_cols = 7;

% Create a zero matrix for the (binary)
% predicate table for the CHSH XOR Game
chsh_xor_game_predicate_table = ...
    zeros(num_rows, num_cols);


% Initialize the index of the current strategy
% for the (binary) predicate table for the CHSH XOR Game
curr_strategy = 1;


% For all the possible (binary)
% values of the input x for Alice
for x = 0:1
    
    % For all the possible (binary)
    % values of the input y for Bob
    for y = 0:1
        
        % For all the possible (binary)
        % values of the output a_x for Alice
        for ax = 0:1
            
            % For all the possible (binary)
            % values of the output b_y for Bob
            for by = 0:1
                
                % Construct the cell values for the (binary)
                % predicate table for the CHSH XOR Game
                chsh_xor_game_predicate_table(curr_strategy, 1) = x;
                chsh_xor_game_predicate_table(curr_strategy, 2) = y;
                chsh_xor_game_predicate_table(curr_strategy, 3) = x && y;
                chsh_xor_game_predicate_table(curr_strategy, 4) = ax;
                chsh_xor_game_predicate_table(curr_strategy, 5) = by;
                chsh_xor_game_predicate_table(curr_strategy, 6) = xor(ax, by);
                chsh_xor_game_predicate_table(curr_strategy, 7) = ...
                    chsh_xor_game_predicate_table(curr_strategy, 3) == ...
                        chsh_xor_game_predicate_table(curr_strategy, 6);
                
                % Increment the current strategy
                % index for Alice and Bob
                curr_strategy = curr_strategy + 1;
                
            end

        end
    
    end

end


% Print the information about the (binary)
% predicate table for the CHSH XOR Game
fprintf(['* The (binary) predicate table for\n' ...
         '  the CHSH XOR Game is the following:\n']);
fprintf('\n');
fprintf('     x     y   x ∧ y   ax    by  ax⊕by  win?\n');
disp(chsh_xor_game_predicate_table);
fprintf('\n');


% Define the number of constraints
% for the Hardy's Paradox regarding
% the bipartite Bell experiment and
% the CHSH Inequality as a XOR Game
num_hardy_paradox_constraints = 4;

% Define the number of combinations for the (binary)
% values for the outputs a_x and b_y for Alice and Bob
num_combinations_ax_by = 4;


% Define the array with the probability
% constraints for the Hardy's Paradox
% regarding the bipartite Bell experiment
% and the CHSH Inequality as a XOR Game
probabilities_hardy_paradox_constraints = ...
    zeros(1, num_hardy_paradox_constraints);


% For all the current strategies for
% Alice and Bob in the CHSH XOR Game
for curr_strategy = 1:num_rows
    
    % Get the (binary) value of the input x for Alice,
    % regarding the current strategy being considered
    x = chsh_xor_game_predicate_table(curr_strategy, 1);
    
    % Get the (binary) value for the input y for Bob,
    % regarding the current strategy being considered
    y = chsh_xor_game_predicate_table(curr_strategy, 2);
    

    % Get the (binary) value of the output a_x for Alice,
    % regarding the current strategy being considered
    ax = chsh_xor_game_predicate_table(curr_strategy, 4);

    % Get the (binary) value for the output b_y for Bob,
    % regarding the current strategy being considered
    by = chsh_xor_game_predicate_table(curr_strategy, 5);
    

    % If the (binary) value of the outputs x, y, a_x and b_y
    % for Alice and Bob are x = 0, y = 0, a_x = 0, and b_y = 0
    if x == 0 && y == 0 && ax == 0 && by == 0
        
        % If the current strategy is a winning
        % one considering the (binary) value of
        % the inputs x and y and the outputs a_x and b_y
        % for Alice and Bob as x = 0, y = 0, a_x = 0, and b_y = 0
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Compute the probability of winning
            % the CHSH XOR Game when the (binary) value of
            % the inputs x and y and the outputs a_x and b_y
            % for Alice and Bob as x = 0, y = 0, a_x = 0, and b_y = 0
            probabilities_hardy_paradox_constraints(1) = ...
                1 / num_combinations_ax_by;

        end

    end
    

    % If the (binary) value of the outputs x, y, a_x and b_y
    % for Alice and Bob are x = 0, y = 1, a_x = 1, and b_y = 1
    if x == 0 && y == 1 && ax == 1 && by == 1
        
        % If the current strategy is a winning
        % one considering the (binary) value of
        % the inputs x and y and the outputs a_x and b_y
        % for Alice and Bob as x = 0, y = 1, a_x = 1, and b_y = 1
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Compute the probability of winning
            % the CHSH XOR Game when the (binary) value of
            % the inputs x and y and the outputs a_x and b_y
            % for Alice and Bob as x = 0, y = 0, a_x = 0, and b_y = 0
            probabilities_hardy_paradox_constraints(2) = ...
                1 / num_combinations_ax_by;

        end

    end
    

    % If the (binary) value of the outputs x, y, a_x and b_y
    % for Alice and Bob are x = 1, y = 0, a_x = 1, and b_y = 1
    if x == 1 && y == 0 && ax == 1 && by == 1
        
        % If the current strategy is a winning
        % one considering the (binary) value of
        % the inputs x and y and the outputs a_x and b_y
        % for Alice and Bob as x = 1, y = 0, a_x = 1, and b_y = 1
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Compute the probability of winning
            % the CHSH XOR Game when the (binary) value of
            % the inputs x and y and the outputs a_x and b_y
            % for Alice and Bob as x = 1, y = 0, a_x = 1, and b_y = 1
            probabilities_hardy_paradox_constraints(3) = ...
                1 / num_combinations_ax_by;

        end

    end
    

    % If the (binary) value of the outputs x, y, a_x and b_y
    % for Alice and Bob are x = 1, y = 1, a_x = 1, and b_y = 1
    if x == 1 && y == 1 && ax == 1 && by == 1
        
        % If the current strategy is a winning
        % one considering the (binary) value of
        % the inputs x and y and the outputs a_x and b_y
        % for Alice and Bob as x = 1, y = 1, a_x = 1, and b_y = 1
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Compute the probability of winning
            % the CHSH XOR Game when the (binary) value of
            % the inputs x and y and the outputs a_x and b_y
            % for Alice and Bob as x = 1, y = 1, a_x = 1, and b_y = 1
            probabilities_hardy_paradox_constraints(4) = ...
                1 / num_combinations_ax_by;

        end

    end

end


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
fprintf('   2) P_{win}(a_x, b_y = 11 | x, y = 01) = 0.0\n');

% Print the 3rd probability constraint
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf('   3) P_{win}(a_x, b_y = 11 | x, y = 10) = 0.0\n');

% Print the 4th probability constraint
% of the Hardy Paradox considering the (binary) values
% for the inputs x and y, and outputs a_x and b_y for Alice and Bob
fprintf('   4) P_{win}(a_x, b_y = 11 | x, y = 11) > 0.0\n');


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
    fprintf('   2) P_{win}(a_x, b_y = 11 | x, y = 01) = 0.0 holds!\n');
    
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
    fprintf(['   2) P_{win}(a_x, b_y = 11 | x, y = 01) = 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 2nd probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 11 | x, y = 01) = %.4f\n', ...
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
    fprintf('   3) P_{win}(a_x, b_y = 11 | x, y = 10) = 0.0 holds!\n');
    
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
    fprintf(['   3) P_{win}(a_x, b_y = 11 | x, y = 10) = 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 3rd probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 11 | x, y = 10) = %.4f\n', ...
            probabilities_hardy_paradox_constraints(3));

    % Print a blank line
    fprintf('\n');

end


% If the 4th probability constraint holds
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(4) > 0.0

    % Print some information about
    % the 4th probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('   4) P_{win}(a_x, b_y = 11 | x, y = 11) > 0.0 holds!\n');
    
    % Print a blank line
    fprintf('\n');

end


% If the 4th probability constraint does not hold
% for the Hardy Paradox considering the (binary) values
% with the inputs x and y, and outputs a_x and b_y for Alice and Bob
if probabilities_hardy_paradox_constraints(4) <= 0.0
    
    % Print some information about
    % the 4th probability constraint not holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf(['   4) P_{win}(a_x, b_y = 11 | x, y = 11) > 0.0 ' ...
             'does not hold since:\n']);
    
    % Print some information about
    % the 4th probability constraint holding
    % for the Hardy Paradox considering the (binary) values
    % with the inputs x and y, and outputs a_x and b_y for Alice and Bob
    fprintf('      => P_{win}(a_x, b_y = 11 | x, y = 11) = %.4f\n', ...
            probabilities_hardy_paradox_constraints(4));

    % Print a blank line
    fprintf('\n');

end