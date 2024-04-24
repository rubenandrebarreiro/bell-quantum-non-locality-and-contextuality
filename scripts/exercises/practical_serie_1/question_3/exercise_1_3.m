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


% Define the number of combinations for the (binary)
% values for the inputs x and y for Alice and Bob
num_combinations_x_y = 4;
num_combinations_ax_by = 4;


% Define the array with the winning cases
% for the combinations of the (binary) values of
% the outputs for Alice and Bob in the CHSH XOR Game
win_cases_chsh_xor_game = zeros(1, num_combinations_ax_by);

% Define the array with the winning probabilities
% for the combinations of the (binary) values of
% the outputs for Alice and Bob in the CHSH XOR Game
win_probabilities_chsh_xor_game = zeros(1, num_combinations_ax_by);


% For all the current strategies for
% Alice and Bob in the CHSH XOR Game
for curr_strategy = 1:num_rows
    
    % Get the (binary) value of the output a_x for Alice,
    % regarding the current strategy being considered
    ax = chsh_xor_game_predicate_table(curr_strategy, 4);

    % Get the (binary) value for the output b_y for Bob,
    % regarding the current strategy being considered
    by = chsh_xor_game_predicate_table(curr_strategy, 5);
    
    
    % If the (binary) value of the outputs a_x and b_y
    % for Alice and Bob are a_x = 0 and b_y = 0
    if ax == 0 && by == 0
        
        % If the current strategy considering
        % the (binary) value of the outputs a_x and b_y
        % for Alice and Bob as a_x = 0 and b_y = 0
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Increment the winning cases count for
            % the combination of the (binary) values of
            % the outputs a_x and b_y for Alice and Bob
            % are a_x = 0 and b_y = 0
            win_cases_chsh_xor_game(1) = ...
                win_cases_chsh_xor_game(1) + 1;

        end

    end
    

    % If the (binary) value of the outputs a_x and b_y
    % for Alice and Bob are a_x = 0 and b_y = 1
    if ax == 0 && by == 1
        
        % If the current strategy considering
        % the (binary) value of the outputs a_x and b_y
        % for Alice and Bob as a_x = 0 and b_y = 1
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Increment the winning cases count for
            % the combination of the (binary) values of
            % the outputs a_x and b_y for Alice and Bob
            % are a_x = 0 and b_y = 1
            win_cases_chsh_xor_game(2) = ...
                win_cases_chsh_xor_game(2) + 1;
        
        end

    end
    
    
    % If the (binary) value of the outputs a_x and b_y
    % for Alice and Bob are a_x = 1 and b_y = 0
    if ax == 1 && by == 0
        
        % If the current strategy considering
        % the (binary) value of the outputs a_x and b_y
        % for Alice and Bob as a_x = 1 and b_y = 0
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Increment the winning cases count for
            % the combination of the (binary) values of
            % the outputs a_x and b_y for Alice and Bob
            % are a_x = 1 and b_y = 0
            win_cases_chsh_xor_game(3) = ...
                win_cases_chsh_xor_game(3) + 1;
            
        end

    end


    % If the (binary) value of the outputs a_x and b_y
    % for Alice and Bob are a_x = 1 and b_y = 1
    if ax == 1 && by == 1
        
        % If the current strategy considering
        % the (binary) value of the outputs a_x and b_y
        % for Alice and Bob as a_x = 1 and b_y = 1
        if chsh_xor_game_predicate_table(curr_strategy, 7)
            
            % Increment the winning cases count for
            % the combination of the (binary) values of
            % the outputs a_x and b_y for Alice and Bob
            % are a_x = 1 and b_y = 1
            win_cases_chsh_xor_game(4) = ...
                win_cases_chsh_xor_game(4) + 1;

        end

    end

end


% Print the headline for the winning probabilities for
% all the combinations of the (binary) values
% for the outputs a_x and b_y for Alice and Bob
fprintf(['* The winning probabilities for all\n' ...
         '  the combinations of a_x and b_y are the following ones:\n']);
fprintf('\n');


% For all the combinations of the (binary) values
% for the outputs a_x and b_y for Alice and Bob
for curr_combination_ax_by = 1:num_combinations_ax_by
    
    % Compute the winning probability for
    % the current combination of the (binary) values
    % for the outputs a_x and b_y for Alice and Bob
    win_probabilities_chsh_xor_game(curr_combination_ax_by) = ...
        win_cases_chsh_xor_game(curr_combination_ax_by) / num_combinations_x_y;
    

    % If the current combination of the (binary) values
    % for the outputs a_x and b_y for Alice and Bob is the 1st one
    if curr_combination_ax_by == 1
        
        % Print the winning probability for
        % the 1st combination of the (binary) values
        % for the outputs a_x and b_y for Alice and Bob
        fprintf('   => P_{win}(a_x, b_y = 00|xy) = %.4f\n', ...
                win_probabilities_chsh_xor_game(curr_combination_ax_by));

        % Print a blank line
        fprintf('\n');

    end
    

    % If the current combination of the (binary) values
    % for the outputs a_x and b_y for Alice and Bob is the 2nd one
    if curr_combination_ax_by == 2
        
        % Print the winning probability for
        % the 2nd combination of the (binary) values
        % for the outputs a_x and b_y for Alice and Bob
        fprintf('   => P_{win}(a_x, b_y = 01|xy) = %.4f\n', ...
                win_probabilities_chsh_xor_game(curr_combination_ax_by));

        % Print a blank line
        fprintf('\n');

    end
    

    % If the current combination of the (binary) values
    % for the outputs a_x and b_y for Alice and Bob is the 3rd one
    if curr_combination_ax_by == 3
        
        % Print the winning probability for
        % the 3rd combination of the (binary) values
        % for the outputs a_x and b_y for Alice and Bob
        fprintf('   => P_{win}(a_x, b_y = 10|xy) = %.4f\n', ...
                win_probabilities_chsh_xor_game(curr_combination_ax_by));

        % Print a blank line
        fprintf('\n');

    end
    
    
    % If the current combination of the (binary) values
    % for the outputs a_x and b_y for Alice and Bob is the 4th one
    if curr_combination_ax_by == 4
        
        % Print the winning probability for
        % the 4th combination of the (binary) values
        % for the outputs a_x and b_y for Alice and Bob
        fprintf('   => P_{win}(a_x, b_y = 11|xy) = %.4f\n', ...
                win_probabilities_chsh_xor_game(curr_combination_ax_by));

        % Print a blank line
        fprintf('\n');

    end

end


% Print a blank line
fprintf('\n');


% Compute the alpha value, considering
% the 1st probability for winning
% the game of the CHSH Inequality
% following a local (classical) strategy
alpha_0 = ( win_probabilities_chsh_xor_game(1) - ...
            ( 1 - win_probabilities_chsh_xor_game(1) ) );

% Compute the 1st local (classical) bound,
% considering the alpha value from
% the 1st probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_0 = 4 * alpha_0;


% Compute the alpha value, considering
% the 2nd probability for winning
% the game of the CHSH Inequality
% following a local (classical) strategy
alpha_1 = ( win_probabilities_chsh_xor_game(2) - ...
            ( 1 - win_probabilities_chsh_xor_game(2) ) );

% Compute the 2nd local (classical) bound,
% considering the alpha value from
% the 2nd probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_1 = 4 * alpha_1;


% Compute the alpha value, considering
% the 3rd probability for winning
% the game of the CHSH Inequality
% following a local (classical) strategy
alpha_2 = ( win_probabilities_chsh_xor_game(3) - ...
            ( 1 - win_probabilities_chsh_xor_game(3) ) );

% Compute the 3rd local (classical) bound,
% considering the alpha value from
% the 3rd probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_2 = 4 * alpha_2;


% Compute the alpha value, considering
% the 4th probability for winning
% the game of the CHSH Inequality
% following a local (classical) strategy
alpha_3 = ( win_probabilities_chsh_xor_game(4) - ...
            ( 1 - win_probabilities_chsh_xor_game(4) ) );

% Compute the 4th local (classical) bound,
% considering the alpha value from
% the 4th probability for winning
% the game of the CHSH Inequality
% following a quantum strategy
quantum_bound_3 = 4 * alpha_3;


% Compute the (classical) local
% upper bound L for the CHSH Inequality
chsh_inequality_local_upper_bound_L = ...
    max( [quantum_bound_0, quantum_bound_1, ...
          quantum_bound_2, quantum_bound_3]);

% Compute the maximum winning probability
% for the CHSH XOR Game, considering all strategies 
max_win_probabilities_chsh_xor_game = ...
    max(win_probabilities_chsh_xor_game);


% Print the headline for the winning strategies for
% all the combinations of the (binary) values of
% the outputs a_x and b_y for Alice and Bob 
fprintf(['* The winning strategies for all\n' ...
         '  the combinations of a_x and b_y are the following ones:\n']);
fprintf('\n');


% If the 1st strategy holds the maximum
% winning probability for the CHSH XOR Game
if win_probabilities_chsh_xor_game(1) == ...
        max_win_probabilities_chsh_xor_game
    
    % Print the 1st strategy if it holds the maximum
    % winning probability for the CHSH XOR Game
    fprintf(['  => A winning strategy is to Alice and Bob\n' ...
             '     agree on a_x = 0 and b_y = 0 a priori!\n']);
    
    % Print a blank line
    fprintf('\n');

end


% If the 2nd strategy holds the maximum
% winning probability for the CHSH XOR Game
if win_probabilities_chsh_xor_game(2) == ...
        max_win_probabilities_chsh_xor_game
    
    % Print the 2nd strategy if it holds the maximum
    % winning probability for the CHSH XOR Game
    fprintf(['  => A winning strategy is to Alice and Bob\n' ...
             '     agree on a_x = 0 and b_y = 1 a priori!\n']);

    % Print a blank line
    fprintf('\n');

end


% If the 3rd strategy holds the maximum
% winning probability for the CHSH XOR Game
if win_probabilities_chsh_xor_game(3) == ...
        max_win_probabilities_chsh_xor_game
    
    % Print the 3rd strategy if it holds the maximum
    % winning probability for the CHSH XOR Game
    fprintf(['  => A winning strategy is to Alice and Bob\n' ...
             '     agree on a_x = 1 and b_y = 0 a priori!\n']);
    
    % Print a blank line
    fprintf('\n');

end


% If the 4th strategy holds the maximum
% winning probability for the CHSH XOR Game
if win_probabilities_chsh_xor_game(4) == ...
        max_win_probabilities_chsh_xor_game
    
    % Print the 4th strategy if it holds the maximum
    % winning probability for the CHSH XOR Game
    fprintf(['  => A winning strategy is to Alice and Bob\n' ...
             '     agree on a_x = 1 and b_y = 1 a priori!\n']);
    
    % Print a blank line
    fprintf('\n');

end


% Print a blank line
fprintf('\n');

% Print of the complete solution of
% the maximization optimization problem
% for the estimation of the upper bound
% for the CHSH Inequality
fprintf(['(Classical) Local Upper Bound L ' ...
         'for CHSH Inequality:\n'])
fprintf(['  a_0 x b_0 + a_0 x b_0 + a_1 x b_0 - a_1 x b_1 ' ...
         '<= L^(C) = %d\n'], ...
         chsh_inequality_local_upper_bound_L);

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
