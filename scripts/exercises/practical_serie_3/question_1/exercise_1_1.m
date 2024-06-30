% Clear Command Window
clc;


% Define the final number of observers
num_final_total_observers = 12;


% Print a blank line
fprintf("\n");

% Print an introductory message about
% solving the (local) classical bounds for
% the Mermin Inequality for multiple parties
fprintf("Solving the (local) classical bounds for\n" + ...
        "the <strong>Mermin Inequality</strong> " + ...
        "for multiple parties <strong>(n = %d)</strong>...\n", ...
        num_final_total_observers);


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("************************************" + ...
        "************************************\n");

% Print two blank lines
fprintf("\n\n");


% Print an introductory message about
% solving the basis case (n = 2) for
% the (local) classical bounds for
% the Mermin Inequality for multiple parties
fprintf("Let's start with the basis case " + ...
        "<strong>(n = 2)</strong>,\n" + ...
        "which is equivalent to the well-known\n" + ...
        "<strong>Clauser-Horne-Shimony-Holt</strong> " + ...
        "<strong>(CHSH) Inequality</strong>:\n");

% Print a blank line
fprintf("\n");

% Print the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F_2, for n = 2
fprintf(" * <strong>F_2 = a_1 * a_2 + a_1 * a'_2 + </strong>" + ...
                 "<strong>a'_1 * a_2 - a'_1 * a'_2</strong>\n");

% Print the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F'_2, for n = 2
fprintf(" * <strong>F'_2 = a'_1 * a'_2 + a'_1 * a_2 + </strong>" + ...
                 "<strong>a_1 * a'_2 - a_1 * a_2</strong>\n");

% Print a blank line
fprintf("\n");

% Print an introductory message about
% solving the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F_2, for n = 2
fprintf("   => Solving the optimization problem on\n" + ...
        "      the (local) classical bound <strong>F_2</strong>...\n");

% Print a blank line
fprintf("\n");

% Compute the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F_2, for n = 2
f_2_1_chsh_inequality_local_bound = ...
    compute_f_2_1_chsh_inequality_local_upper_bound_L();

% Print an introductory message about
% solving the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F'_2, for n = 2
fprintf("   => Solving the optimization problem on\n" + ...
        "      the (local) classical bound <strong>F'_2</strong>...\n");

% Print a blank line
fprintf("\n");

% Compute the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F'_2, for n = 2
f_2_2_chsh_inequality_local_bound = ...
    compute_f_2_2_chsh_inequality_local_upper_bound_L();

% Print a message about the result of
% the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F'_2, for n = 2
fprintf("   => The (local) classical bound <strong>F_2</strong> of\n" + ...
        "      the <strong>Mermin Inequality</strong> " + ...
              "is equal to: %.4f\n", ...
        f_2_1_chsh_inequality_local_bound);


% Print two blank lines
fprintf("\n\n");

% Print a separator
fprintf("************************************" + ...
        "************************************\n");

% Print two blank lines
fprintf("\n\n");


% If the final number of
% observers is equal to two
if(num_final_total_observers == 2)
    
    % Print an informative message about
    % not existing more remaining cases to solve
    fprintf("There are no more remaining cases to solve...\n");

% If the final number of
% observers is greater than two
elseif(num_final_total_observers > 2)

    % Print an informative message about
    % continuing to solve the remaining cases
    fprintf("Let's continue with the remaining cases " + ...
            "<strong>(n > 2)</strong>...\n");
    
    % Print two blank lines
    fprintf("\n\n");
    
    % Print a separator
    fprintf("************************************" + ...
            "************************************\n");
    
    % Print two blank lines
    fprintf("\n\n");
    
    
    % Initialize a vector to keep
    % the (local) classical bounds of
    % the Mermin Inequality F_n,
    % for the each on of the several
    % configurations on the number of parties
    f_n_1_chsh_inequality_local_bounds = ...
        zeros(num_final_total_observers);
    
    % Store the (local) classical bounds of
    % the Mermin Inequality F_2, for n = 2
    f_n_1_chsh_inequality_local_bounds(2) = ...
        f_2_1_chsh_inequality_local_bound;


    % Initialize a vector to keep
    % the (local) classical bounds of
    % the Mermin Inequality F'_n,
    % for the each on of the several
    % configurations on the number of parties
    f_n_2_chsh_inequality_local_bounds = ...
        zeros(num_final_total_observers);
    
    % Store the (local) classical bounds of
    % the Mermin Inequality F'_2, for n = 2
    f_n_2_chsh_inequality_local_bounds(2) = ...
        f_2_2_chsh_inequality_local_bound;

    
    % For each configuration for the number of
    % observers, greater than two (basis casis)
    for curr_num_total_observers = 3:num_final_total_observers
        
        % Print an introductory message about the current
        % configuration for the number of observers
        fprintf("For <strong>(n = %d)</strong>:\n", ...
                curr_num_total_observers);

        % Print a blank line
        fprintf("\n");
        
        % Print the mathematical expression
        % for the (local) classical bound of
        % the Mermin Inequality F_n, for n > 2
        fprintf(" * <strong>F_%d = 1/2 * (a_%d + a'_%d) * F_%d + </strong>" + ...
                          "<strong>1/2 * (a_%d - a'_%d) * F'_%d</strong>\n", ...
                curr_num_total_observers, curr_num_total_observers, ...
                curr_num_total_observers, (curr_num_total_observers - 1), ...
                curr_num_total_observers, curr_num_total_observers, ...
                (curr_num_total_observers - 1));
        
        % Print the mathematical expression
        % for the (local) classical bound of
        % the Mermin Inequality F'_n, for n > 2
        fprintf(" * <strong>F'_%d = 1/2 * (a'_%d + a_%d) * F'_%d + </strong>" + ...
                           "<strong>1/2 * (a'_%d - a_%d) * F_%d</strong>\n", ...
                curr_num_total_observers, curr_num_total_observers, ...
                curr_num_total_observers, (curr_num_total_observers - 1), ...
                curr_num_total_observers, curr_num_total_observers, ...
                (curr_num_total_observers - 1));


        % Print a blank line
        fprintf("\n");

        
        % Retrieve the (local) classical bound of
        % the Mermin Inequality F_(n - 1), for n > 2,
        % for the previous configuration
        % for the number of observers
        f_n_minus_1_1_chsh_inequality_local_bound = ...
            f_n_1_chsh_inequality_local_bounds...
                ( (curr_num_total_observers - 1 ) );
            
        % Retrieve the (local) classical bound of
        % the Mermin Inequality F'_(n - 1), for n > 2,
        % for the previous configuration
        % for the number of observers
        f_n_minus_1_2_chsh_inequality_local_bound = ...
            f_n_2_chsh_inequality_local_bounds...
                ( (curr_num_total_observers - 1 ) );
        
        
        % Print an introductory message about
        % solving the mathematical expression
        % for the (local) classical bound of
        % the Mermin Inequality F_n, for n > 2
        fprintf("   => Solving the optimization problem on\n" + ...
                "      the (local) classical bound " + ...
                      "<strong>F_%d</strong>...\n", ...
                curr_num_total_observers);
        
        % Print a blank line
        fprintf("\n");

        % Compute the (local) classical bound of
        % the Mermin Inequality F_n by induction, for n > 2
        f_n_1_chsh_inequality_local_bound = ...
            compute_f_n_1_chsh_inequality_local_upper_bound_L...
                (curr_num_total_observers, ...
                 f_n_minus_1_1_chsh_inequality_local_bound, ...
                 f_n_minus_1_2_chsh_inequality_local_bound);
        
        
        % Print an introductory message about
        % solving the mathematical expression
        % for the (local) classical bound of
        % the Mermin Inequality F'_n, for n > 2
        fprintf("   => Solving the optimization problem on\n" + ...
                "      the (local) classical bound " + ...
                      "<strong>F'_%d</strong>...\n", ...
                curr_num_total_observers);
        
        % Print a blank line
        fprintf("\n");
        
        % Compute the (local) classical bound of
        % the Mermin Inequality F_n by induction, for n > 2
        f_n_2_chsh_inequality_local_bound = ...
            compute_f_n_2_chsh_inequality_local_upper_bound_L...
                (curr_num_total_observers, ...
                 f_n_minus_1_1_chsh_inequality_local_bound, ...
                 f_n_minus_1_2_chsh_inequality_local_bound);
        
        
        % Print a message about the result of
        % the basis mathematical expression
        % for the (local) classical bound of
        % the Mermin Inequality F_n, for n > 2
        fprintf("   => The (local) classical bound " + ...
                      "<strong>F_%d</strong> of\n" + ...
                "      the <strong>Mermin Inequality</strong> " + ...
                      "is equal to: %.4f\n", ...
                curr_num_total_observers, ...
                f_n_1_chsh_inequality_local_bound);
        
        % Print two blank lines
        fprintf("\n\n");
        
        
        % Store the (local) classical bounds of
        % the Mermin Inequality F_n, for n > 2,
        % considering the current configuration
        % for the number of observers
        f_n_1_chsh_inequality_local_bounds...
            (curr_num_total_observers) = ...
                f_n_1_chsh_inequality_local_bound;
        
        % Store the (local) classical bounds of
        % the Mermin Inequality F'_n, for n > 2,
        % considering the current configuration
        % for the number of observers
        f_n_2_chsh_inequality_local_bounds...
            (curr_num_total_observers) = ...
                f_n_2_chsh_inequality_local_bound;
        
    end
    
end


% Print a separator
fprintf("************************************" + ...
        "************************************\n");

% Print a blank line
fprintf("\n");


% Define a function to compute
% the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F_2, for n = 2
function f_2_1_chsh_inequality_local_upper_bound_L = ...
         compute_f_2_1_chsh_inequality_local_upper_bound_L()
    
    % Definition of the maximization optimization
    % problem for the estimation value of the (local)
    % classical upper bound for the Mermin Inequality F_2
    f_2_1_chsh_inequality_local_bound_optimization_prob = ...
        optimproblem('ObjectiveSense', 'max');
    
    % Definition of the optimization variables/unknowns
    % to be optimized in the maximization optimization
    % problem for the estimation value of the (local)
    % classical upper bound for the Mermin Inequality F_2
    a_1_1 = optimvar('a_1_1', 1, 1, 'Type', 'integer');
    a_1_2 = optimvar('a_1_2', 1, 1, 'Type', 'integer');
    a_2_1 = optimvar('a_2_1', 1, 1, 'Type', 'integer');
    a_2_2 = optimvar('a_2_2', 1, 1, 'Type', 'integer');
    
    % Definition of the objective function of
    % the maximization optimization problem
    % for the estimation value of the (local)
    % classical upper bound for the Mermin Inequality F_2
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Objective = a_1_1(1) * a_2_1(1) + a_1_1(1) * a_2_2(1) + ... 
                     a_1_2(1) * a_2_1(1) - a_1_2(1) * a_2_2(1);
    
    
    % Definition of the constraints
    % for the a_1 value E {-1,+1}
    constraint_1 = a_1_1(1) >= -1;
    constraint_2 = a_1_1(1) <= 1;
    constraint_3 = a_1_1(1)^2 <= 1;
    constraint_4 = a_1_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a'_1 value E {-1,+1}
    constraint_5 = a_1_2(1) >= -1;
    constraint_6 = a_1_2(1) <= 1;
    constraint_7 = a_1_2(1)^2 <= 1;
    constraint_8 = a_1_2(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_2 value E {-1,+1}
    constraint_9 = a_2_1(1) >= -1;
    constraint_10 = a_2_1(1) <= 1;
    constraint_11 = a_2_1(1)^2 <= 1;
    constraint_12 = a_2_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a'_2 value E {-1,+1}
    constraint_13 = a_2_2(1) >= -1;
    constraint_14 = a_2_2(1) <= 1;
    constraint_15 = a_2_2(1)^2 <= 1;
    constraint_16 = a_2_2(1)^2 >= 1;
    
    
    % Setup of all the constraints
    % for the correlator estimation
    % values a_1, a'_1, a_2, and a'_2,
    % for the maximization optimization
    % problem for the estimation of the (local)
    % classical upper bound for the Mermin Inequality F_2
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons1 = constraint_1;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons2 = constraint_2;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons3 = constraint_3;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons4 = constraint_4;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons5 = constraint_5;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons6 = constraint_6;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons7 = constraint_7;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons8 = constraint_8;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons9 = constraint_9;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons10 = constraint_10;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons11 = constraint_11;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons12 = constraint_12;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons13 = constraint_13;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons14 = constraint_14;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons15 = constraint_15;
    f_2_1_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons16 = constraint_16;
    

    % Initialize the vectors for the initial point,
    % considering the variables regarding the correlator
    % estimation values a_1, a'_1, a_2, and a'_2,
    % for the maximization optimization problem
    % for the estimation of the (local) classical
    % upper bound for the Mermin Inequality F_2, for n = 2
    x0.a_1_1 = ones(size(a_1_1)) * (randi(13) - 13);
    x0.a_1_2 = ones(size(a_1_1)) * (randi(13) - 13);
    x0.a_2_1 = ones(size(a_1_1)) * (randi(13) - 13);
    x0.a_2_2 = ones(size(a_1_1)) * (randi(13) - 13);
    

    % Define some optimization options for the optimization
    % algorithm/procedure, in order to disable the verbose outputs
    optimization_options = optimoptions("ga", "Display", "none");

    % Solve the maximization optimization problem
    % for the estimation of the (local) classical
    % upper bound for the Mermin Inequality F_2, for n = 2
    f_2_1_chsh_inequality_local_bound_optimization_prob_sol = ...
        solve(f_2_1_chsh_inequality_local_bound_optimization_prob, ...
              x0, "Options", optimization_options);
    
    % Compute the a_1, a'_1, a_2, and a'_2 values
    % for the Mermin Inequality F_2 that allow to
    % estimate the (local) classical upper bound L, for n = 2
    f_2_1_chsh_inequality_local_bound_a_1_1_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_1_1;
    f_2_1_chsh_inequality_local_bound_a_1_2_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_1_2;
    f_2_1_chsh_inequality_local_bound_a_2_1_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_2_1;
    f_2_1_chsh_inequality_local_bound_a_2_2_value = ...
        f_2_1_chsh_inequality_local_bound_optimization_prob_sol.a_2_2;
    
    % Compute the individual expectation values
    % for the Mermin Inequality F_2 that allow to
    % estimate the (local) classical upper bound L, for n = 2
    a_1_1 = f_2_1_chsh_inequality_local_bound_a_1_1_value(1);
    a_1_2 = f_2_1_chsh_inequality_local_bound_a_1_2_value(1);
    a_2_1 = f_2_1_chsh_inequality_local_bound_a_2_1_value(1);
    a_2_2 = f_2_1_chsh_inequality_local_bound_a_2_2_value(1);
    
    % Compute the (local) classical upper bound L
    % for the Mermin Inequality F_2, for n = 2
    f_2_1_chsh_inequality_local_upper_bound_L = ...
        a_1_1 * a_2_1 + a_1_1 * a_2_2 + ...
        a_1_2 * a_2_1 - a_1_2 * a_2_2;

end


% Define a function to compute
% the basis mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F'_2, for n = 2
function f_2_2_chsh_inequality_local_upper_bound_L = ...
         compute_f_2_2_chsh_inequality_local_upper_bound_L()
        
    % Definition of the maximization optimization
    % problem for the estimation value of the (local)
    % classical upper bound for the Mermin Inequality F'_2
    f_2_2_chsh_inequality_local_bound_optimization_prob = ...
        optimproblem('ObjectiveSense', 'max');
    
    % Definition of the optimization variables/unknowns
    % to be optimized in the maximization optimization
    % problem for the estimation value of the (local)
    % classical upper bound for the Mermin Inequality F'_2
    a_1_1 = optimvar('a_1_1', 1, 1, 'Type', 'integer');
    a_1_2 = optimvar('a_1_2', 1, 1, 'Type', 'integer');
    a_2_1 = optimvar('a_2_1', 1, 1, 'Type', 'integer');
    a_2_2 = optimvar('a_2_2', 1, 1, 'Type', 'integer');
    
    % Definition of the objective function of
    % the maximization optimization problem
    % for the estimation value of the (local)
    % classical upper bound for the Mermin Inequality F'_2
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Objective = a_1_2(1) * a_2_2(1) + a_1_2(1) * a_2_1(1) + ... 
                     a_1_1(1) * a_2_2(1) - a_1_1(1) * a_2_1(1);
    
    % Definition of the constraints
    % for the a_1_1 value E {-1,+1}
    constraint_1 = a_1_1(1) >= -1;
    constraint_2 = a_1_1(1) <= 1;
    constraint_3 = a_1_1(1)^2 <= 1;
    constraint_4 = a_1_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_1_2 value E {-1,+1}
    constraint_5 = a_1_2(1) >= -1;
    constraint_6 = a_1_2(1) <= 1;
    constraint_7 = a_1_2(1)^2 <= 1;
    constraint_8 = a_1_2(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_2_1 value E {-1,+1}
    constraint_9 = a_2_1(1) >= -1;
    constraint_10 = a_2_1(1) <= 1;
    constraint_11 = a_2_1(1)^2 <= 1;
    constraint_12 = a_2_1(1)^2 >= 1;
    
    % Definition of the constraints
    % for the a_2_2 value E {-1,+1}
    constraint_13 = a_2_2(1) >= -1;
    constraint_14 = a_2_2(1) <= 1;
    constraint_15 = a_2_2(1)^2 <= 1;
    constraint_16 = a_2_2(1)^2 >= 1;
    
    
    % Setup of all the constraints
    % for the correlator estimation
    % values a_1, a'_1, a_2, and a'_2,
    % for the maximization optimization
    % problem for the estimation of the (local)
    % classical upper bound for the Mermin Inequality F'_2
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons1 = constraint_1;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons2 = constraint_2;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons3 = constraint_3;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons4 = constraint_4;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons5 = constraint_5;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons6 = constraint_6;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons7 = constraint_7;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons8 = constraint_8;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons9 = constraint_9;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons10 = constraint_10;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons11 = constraint_11;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons12 = constraint_12;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons13 = constraint_13;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons14 = constraint_14;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons15 = constraint_15;
    f_2_2_chsh_inequality_local_bound_optimization_prob...
        .Constraints.cons16 = constraint_16;
    

    % Initialize the vectors for the initial point,
    % considering the variables regarding the correlator
    % estimation values a_1, a'_1, a_2, and a'_2,
    % for the maximization optimization problem
    % for the estimation of the (local) classical
    % upper bound for the Mermin Inequality F_2, for n = 2
    x0.a_1_1 = ones(size(a_1_1)) * (randi(13) - 13);
    x0.a_1_2 = ones(size(a_1_1)) * (randi(13) - 13);
    x0.a_2_1 = ones(size(a_1_1)) * (randi(13) - 13);
    x0.a_2_2 = ones(size(a_1_1)) * (randi(13) - 13);
    
    
    % Define some optimization options for the optimization
    % algorithm/procedure, in order to disable the verbose outputs
    optimization_options = optimoptions("ga", "Display", "none");
    

    % Solve the maximization optimization problem
    % for the estimation of the (local) classical
    % upper bound for the Mermin Inequality F'_2, for n = 2
    f_2_2_chsh_inequality_local_bound_optimization_prob_sol = ...
        solve(f_2_2_chsh_inequality_local_bound_optimization_prob, ...
              x0, "Options", optimization_options);
    
    % Compute the a_1, a'_1, a_2, and a'_2 values
    % for the Mermin Inequality F'_2 that allow to
    % estimate the (local) classical upper bound L, for n = 2
    f_2_2_chsh_inequality_local_bound_a_1_1_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_1_1;
    f_2_2_chsh_inequality_local_bound_a_1_2_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_1_2;
    f_2_2_chsh_inequality_local_bound_a_2_1_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_2_1;
    f_2_2_chsh_inequality_local_bound_a_2_2_value = ...
        f_2_2_chsh_inequality_local_bound_optimization_prob_sol.a_2_2;
    

    % Compute the individual expectation values
    % for the Mermin Inequality F'_2 that allow to
    % estimate the (local) classical upper bound L, for n = 2
    a_1_1 = f_2_2_chsh_inequality_local_bound_a_1_1_value(1);
    a_1_2 = f_2_2_chsh_inequality_local_bound_a_1_2_value(1);
    a_2_1 = f_2_2_chsh_inequality_local_bound_a_2_1_value(1);
    a_2_2 = f_2_2_chsh_inequality_local_bound_a_2_2_value(1);

    % Compute the (local) classical upper bound L
    % for the Mermin Inequality F'_2, for n = 2
    f_2_2_chsh_inequality_local_upper_bound_L = ...
        a_1_2 * a_2_2 + a_1_2 * a_2_1 + ...
        a_1_1 * a_2_2 - a_1_1 * a_2_1;
    
    
end


% Define a function to compute
% the general mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F_n, for n > 2
function f_n_1_chsh_inequality_local_upper_bound_L = ...
         compute_f_n_1_chsh_inequality_local_upper_bound_L...
         (num_total_observers, ...
          f_n_minus_1_1_chsh_inequality_local_bound, ...
          f_n_minus_1_2_chsh_inequality_local_bound)
    
    % If the given total number of
    % observers is greater than two
    if num_total_observers > 2

        % Definition of the maximization optimization
        % problem for the estimation value of the (local)
        % classical upper bound for the Mermin Inequality F_n
        f_n_1_chsh_inequality_local_bound_optimization_prob = ...
            optimproblem('ObjectiveSense', 'max');
       

        % Definition of the optimization variables/unknowns
        % to be optimized in the maximization optimization
        % problem for the estimation value of the (local)
        % classical upper bound for the Mermin Inequality F_n
        a_n_1 = optimvar('a_n_1', 1, 1, 'Type', 'integer');
        a_n_2 = optimvar('a_n_2', 1, 1, 'Type', 'integer');
 

        % Definition of the objective function of
        % the maximization optimization problem
        % for the estimation value of the (local)
        % classical upper bound for the Mermin Inequality F_n
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Objective = (1 / 2) * (a_n_1 + a_n_2) * ...
                            f_n_minus_1_1_chsh_inequality_local_bound + ...
                         (1 / 2) * (a_n_1 - a_n_2) * ...
                            f_n_minus_1_2_chsh_inequality_local_bound;
        
        % Definition of the constraints
        % for the a_n value E {-1,+1}
        constraint_1 = a_n_1(1) >= -1;
        constraint_2 = a_n_1(1) <= 1;
        constraint_3 = a_n_1(1)^2 <= 1;
        constraint_4 = a_n_1(1)^2 >= 1;
        
        % Definition of the constraints
        % for the a'_n value E {-1,+1}
        constraint_5 = a_n_2(1) >= -1;
        constraint_6 = a_n_2(1) <= 1;
        constraint_7 = a_n_2(1)^2 <= 1;
        constraint_8 = a_n_2(1)^2 >= 1;
        

        % Setup of all the constraints for the correlator
        % estimation values a_n, and a'_n, for the maximization
        % optimization problem for the estimation of the (local)
        % classical upper bound for the Mermin Inequality F_n
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons1 = constraint_1;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons2 = constraint_2;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons3 = constraint_3;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons4 = constraint_4;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons5 = constraint_5;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons6 = constraint_6;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons7 = constraint_7;
        f_n_1_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons8 = constraint_8;
        

        % Initialize the vectors for the initial point,
        % considering the variables regarding the correlator
        % estimation values a_n and a'_n, for the maximization
        % optimization problem for the estimation of the (local)
        % classical upper bound for the Mermin Inequality F_n, for n > 2
        x0.a_n_1 = ones(size(a_n_1)) * (randi(13) - 13);
        x0.a_n_2 = ones(size(a_n_1)) * (randi(13) - 13);
        

        % Define some optimization options for the optimization
        % algorithm/procedure, in order to disable the verbose outputs 
        optimization_options = optimoptions("ga", "Display", "none");
        
        
        % Solve the maximization optimization problem
        % for the estimation of the (local) classical
        % upper bound for the Mermin Inequality F_n, for n > 2
        f_n_1_chsh_inequality_local_bound_optimization_prob_sol = ...
            solve(f_n_1_chsh_inequality_local_bound_optimization_prob, ...
                  x0, 'Options', optimization_options);
        
        % Compute the a_n and a'_n values for the
        % Mermin Inequality F_n that allow to estimate
        % the (local) classical upper bound L, for n > 2
        f_n_1_chsh_inequality_local_bound_a_n_1_value = ...
            f_n_1_chsh_inequality_local_bound_optimization_prob_sol.a_n_1;
        f_n_1_chsh_inequality_local_bound_a_n_2_value = ...
            f_n_1_chsh_inequality_local_bound_optimization_prob_sol.a_n_2;
        
        % Compute the individual expectation values
        % for the Mermin Inequality F_n that allow to
        % estimate the (local) classical upper bound L, for n > 2
        a_n_1 = f_n_1_chsh_inequality_local_bound_a_n_1_value(1);
        a_n_2 = f_n_1_chsh_inequality_local_bound_a_n_2_value(1);
        
        % Compute the (local) classical upper bound L
        % for the Mermin Inequality F_n, for n > 2
        f_n_1_chsh_inequality_local_upper_bound_L = ...
            (1 / 2) * (a_n_1 + a_n_2) * ...
                f_n_minus_1_1_chsh_inequality_local_bound + ...
            (1 / 2) * (a_n_1 - a_n_2) * ...
                f_n_minus_1_2_chsh_inequality_local_bound;
        
    end

end


% Define a function to compute
% the general mathematical expression
% for the (local) classical bound of
% the Mermin Inequality F'_n, for n > 2
function f_n_2_chsh_inequality_local_upper_bound_L = ...
        compute_f_n_2_chsh_inequality_local_upper_bound_L...
        (num_total_observers, ...
         f_n_minus_1_1_chsh_inequality_local_bound, ...
         f_n_minus_1_2_chsh_inequality_local_bound)
    
    % If the given total number of
    % observers is greater than two
    if num_total_observers > 2

        % Definition of the maximization optimization
        % problem for the estimation value of the (local)
        % classical upper bound for the Mermin Inequality F'_n
        f_n_2_chsh_inequality_local_bound_optimization_prob = ...
            optimproblem('ObjectiveSense', 'max');
        
        % Definition of the optimization variables/unknowns
        % to be optimized in the maximization optimization
        % problem for the estimation value of the (local)
        % classical upper bound for the Mermin Inequality F'_n
        a_n_1 = optimvar('a_n_1', 1, 1, 'Type', 'integer');
        a_n_2 = optimvar('a_n_2', 1, 1, 'Type', 'integer');
        
        % Definition of the objective function of
        % the maximization optimization problem
        % for the estimation value of the (local)
        % classical upper bound for the Mermin Inequality F'_n
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Objective = (1 / 2) * (a_n_2 + a_n_1) * ...
                            f_n_minus_1_2_chsh_inequality_local_bound + ...
                         (1 / 2) * (a_n_2 - a_n_1) * ...
                            f_n_minus_1_1_chsh_inequality_local_bound;
        
        % Definition of the constraints
        % for the a_n_1 value E {-1,+1}
        constraint_1 = a_n_1(1) >= -1;
        constraint_2 = a_n_1(1) <= 1;
        constraint_3 = a_n_1(1)^2 <= 1;
        constraint_4 = a_n_1(1)^2 >= 1;
        
        % Definition of the constraints
        % for the a_n_2 value E {-1,+1}
        constraint_5 = a_n_2(1) >= -1;
        constraint_6 = a_n_2(1) <= 1;
        constraint_7 = a_n_2(1)^2 <= 1;
        constraint_8 = a_n_2(1)^2 >= 1;
        
        % Setup of all the constraints for the correlator
        % estimation values a_n, and a'_n, for the maximization
        % optimization problem for the estimation of the (local)
        % classical upper bound for the Mermin Inequality F'_n
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons1 = constraint_1;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons2 = constraint_2;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons3 = constraint_3;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons4 = constraint_4;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons5 = constraint_5;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons6 = constraint_6;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons7 = constraint_7;
        f_n_2_chsh_inequality_local_bound_optimization_prob...
            .Constraints.cons8 = constraint_8;


        % Initialize the vectors for the initial point,
        % considering the variables regarding the correlator
        % estimation values a_n and a'_n, for the maximization
        % optimization problem for the estimation of the (local)
        % classical upper bound for the Mermin Inequality F'_n, for n > 2
        x0.a_n_1 = ones(size(a_n_1)) * (randi(13) - 13);
        x0.a_n_2 = ones(size(a_n_1)) * (randi(13) - 13);
    

        % Define some optimization options for the optimization
        % algorithm/procedure, in order to disable the verbose outputs 
        optimization_options = optimoptions("ga", "Display", "none");
        
        
        % Solve the maximization optimization problem
        % for the estimation of the (local) classical
        % upper bound for the Mermin Inequality F'_n, for n > 2
        f_n_2_chsh_inequality_local_bound_optimization_prob_sol = ...
            solve(f_n_2_chsh_inequality_local_bound_optimization_prob, ...
                  x0, 'Options', optimization_options);

        % Compute the a_n and a'_n values for the
        % Mermin Inequality F'_n that allow to estimate
        % the (local) classical upper bound L, for n > 2
        f_n_2_chsh_inequality_local_bound_a_n_1_value = ...
            f_n_2_chsh_inequality_local_bound_optimization_prob_sol.a_n_1;
        f_n_2_chsh_inequality_local_bound_a_n_2_value = ...
            f_n_2_chsh_inequality_local_bound_optimization_prob_sol.a_n_2;
        
        % Compute the individual expectation values
        % for the Mermin Inequality F'_n that allow to
        % estimate the (local) classical upper bound L, for n > 2
        a_n_1 = f_n_2_chsh_inequality_local_bound_a_n_1_value(1);
        a_n_2 = f_n_2_chsh_inequality_local_bound_a_n_2_value(1);
        
        % Compute the (local) classical upper bound L
        % for the Mermin Inequality F'_n, for n > 2
        f_n_2_chsh_inequality_local_upper_bound_L = ...
            (1 / 2) * (a_n_2 + a_n_1) * ...
                f_n_minus_1_2_chsh_inequality_local_bound + ...
            (1 / 2) * (a_n_2 - a_n_1) * ...
                f_n_minus_1_1_chsh_inequality_local_bound;
        
    end

end