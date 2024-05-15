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


% Print information about all the configurations
% which will be used to solve the Mermin Inequality
% using its exponent/power formulation,
% with the respective ranging values
fprintf(['Solving the Mermin Inequality ' ...
         'using its exponent/power formulation\n' ...
         'for a number of observers a_n, ' ...
         'ranging from n = %d to n = %d:\n'], ...
         num_initial_total_observers, ...
         num_final_total_observers);


% Print a blank line
fprintf('\n');


% For each defined configuration
% for the total number of observers
for num_total_observers = ...
    num_initial_total_observers:num_final_total_observers
    
    % Compute the (Classical) Local Upper Bound L
    % for the Mermin Inequality configured
    % with the total number of observers
    % defined previously
    f_n_mermin_inequality_local_upper_bound_L = ...
        compute_f_n_mermin_inequality_local_upper_bound_L...
            (num_total_observers);
    
    % Print some information about solving
    % the (Classical) Local Upper Bound L of
    % the Mermin Inequality with the current
    % total number of observers configured
    fprintf(['* (Classical) Local Upper Bound L ' ...
             'for Mermin Inequality\n' ...
             '  with a number of observers n = %d:\n'], ...
             num_total_observers);
    
    % Print a blank line
    fprintf('\n');


    % If the current total number of observers
    % configured for the Mermin Inequality is even
    if mod(num_total_observers, 2) == 0
        
        % Print more information about solving
        % the (Classical) Local Upper Bound L of
        % the Mermin Inequality with the current
        % (even) total number of observers configured
        fprintf(['  => F_n = 2^(n / 2), ' ...
                 'with n even (n = %d):\n'], ...
                num_total_observers);

        % Print more information about solving
        % the (Classical) Local Upper Bound L of
        % the Mermin Inequality with the current
        % (even) total number of observers configured,
        % with all the respective calculations
        fprintf(['     * F_n = 2^(n / 2) (=) ' ...
                 'F_%d = 2^(%d / 2) (=)\n' ...
                 '       F_%d = 2^%d (=) F_%d = %d\n'], ...
                num_total_observers, num_total_observers, ...
                num_total_observers, ...
                ( num_total_observers / 2 ), ...
                num_total_observers, ...
                f_n_mermin_inequality_local_upper_bound_L);
        
        % Print a blank line
        fprintf('\n');

    end
    
    
    % If the current total number of observers
    % configured for the Mermin Inequality is odd
    if mod(num_total_observers, 2) ~= 0
        
        % Print more information about solving
        % the (Classical) Local Upper Bound L of
        % the Mermin Inequality with the current
        % (odd) total number of observers configured
        fprintf(['  => F_n = 2^( ( n - 1 ) / 2 ), ' ...
                 'with n odd (n = %d):\n'], ...
                num_total_observers);

        % Print more information about solving
        % the (Classical) Local Upper Bound L of
        % the Mermin Inequality with the current
        % (odd) total number of observers configured,
        % with all the respective calculations
        fprintf(['     * F_n = 2^( ( n - 1 ) / 2 ) (=) ' ...
                 'F_%d = 2^( ( %d - 1 ) / 2 ) (=)\n' ...
                 '       F_%d = 2^%d (=) F_%d = %d\n'], ...
                num_total_observers, num_total_observers, ...
                num_total_observers, ...
                ( ( num_total_observers - 1 ) / 2 ), ...
                num_total_observers, ...
                f_n_mermin_inequality_local_upper_bound_L);
        
        % Print a blank line
        fprintf('\n');
        
    end

end


% Define the function to compute 
% the (Classical) Local Upper Bound L of
% the Mermin Inequality with a given
% total number of observers configured,
% using its exponent/power formulation
function f_n_mermin_inequality_local_upper_bound_L = ...
         compute_f_n_mermin_inequality_local_upper_bound_L(n)
    
    % If the given total number of observers is even
    if mod(n, 2) == 0
    
        % Compute the (Classical) Local Upper Bound L of
        % the Mermin Inequality with a given (even)
        % total number of observers configured,
        % using its exponent/power formulation 
        f_n_mermin_inequality_local_upper_bound_L = ...
            2^(n / 2);

    end

    
    % If the given total number of observers is odd
    if mod(n, 2) ~= 0
        
        % Compute the (Classical) Local Upper Bound L of
        % the Mermin Inequality with a given (odd)
        % total number of observers configured,
        % using its exponent/power formulation 
        f_n_mermin_inequality_local_upper_bound_L = ...
            2^( ( n - 1 ) / 2 );

    end

end