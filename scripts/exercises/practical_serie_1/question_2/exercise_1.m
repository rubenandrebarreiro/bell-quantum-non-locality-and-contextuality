% Clear Command Window
clc;


% Create the symbolic numbers
% for the complex coefficients of
% a two-dimensional quantum state
% (i.e., alpha and beta)
syms alpha_0 alpha_1 beta_0 beta_1;

% Create the symbolic trigonometric angles
% theta and phi for the two-dimensional
% quantum state in the trigonometric form
syms theta phi;


% Create the assumptions regarding
% the real and imaginary parts of
% the complex coefficient alpha
% being real numbers
assume(alpha_0, 'real');
assume(alpha_1, 'real');

% Create the assumptions regarding
% the real and imaginary parts of
% the complex coefficient beta
% being real numbers
assume(beta_0, 'real');
assume(beta_1, 'real');


% Create the complex number alpha
% as the 1st complex coefficient
% for the two-dimensional quantum state
alpha = alpha_0 * 1 + alpha_1 * 1j;

% Create the complex number beta
% as the 2nd complex coefficient for
% the two-dimensional quantum state
beta = beta_0 * 1 + beta_1 * 1j;

% Create the assumption for
% the linear and probabilistic combination of
% magnitudes for the complex coefficients of
% a two-dimensional quantum state 
assume( abs(alpha)^2 + abs(beta)^2 == 1.0 );

% Create the assumptions for the interval of
% the theta angle of the two-dimensional
% quantum state in the trigonometric form
assume(theta >= 0);
assume(theta < pi);

% Create the assumptions for the interval of
% the phi angle of the two-dimensional
% quantum state in the trigonometric form
assume(phi >= 0);
assume(phi < pi);

% Compute the complex coefficients of
% a two-dimensional quantum state as
% trigonometric identities in function
% of the theta and phi angles
alpha_trig = cos(theta / 2);
beta_trig = sin(theta / 2) * exp(phi * 1j);


% Create the ket vector for
% the two-dimensional quantum state
% denoted as |psi> = ([alpha, beta])^T
ket_psi = [alpha beta].';

% Compute the two-dimensional quantum state
% in its trigonometric form in function of
% the respective theta and phi angles
ket_psi_trig = [alpha_trig beta_trig].';

% Compute the two-dimensional quantum states
% |0> and |1> corresponding to the classical states
% (i.e., bits) 0 and 1, respectively
ket_0 = [1 0].';
ket_1 = [0 1].';

% Print a blank line
fprintf('\n');

% If the linear complex combination of
% alpha x |0> + beta x |1> is equal to
% the general form of a two-dimensional
% quantum state |psi>
if isequal(alpha * ket_0 + ...
           beta * ket_1, ...
           ket_psi)

    % Print information about
    % the linear complex combination of
    % alpha x |0> + beta x |1> being equal to
    % the general form of a two-dimensional
    % quantum state |psi>
    fprintf([' alpha * |0> + beta * |1> =\n' ...
             ' = alpha * [1 0]^T + beta * [0 1]^T =' ...
             ' [alpha beta]^T\n']);

    % Print a blank line
    fprintf('\n');
    
    % Print more information about
    % the linear complex combination of
    % alpha x |0> + beta x |1> being equal to
    % the general form of a two-dimensional
    % quantum state |psi>
    fprintf(['The condition of the linear complex\n' ...
             'combination of alpha x |0> + beta x |1>\n' ...
             'being equal to the general form of\n' ...
             'a two-dimensional quantum state |psi> holds']);

    % Print a blank line
    fprintf('\n');

end


% Print a blank line
fprintf('\n');


% If the linear and probabilistic combination of
% magnitudes for the complex coefficients of
% a two-dimensional quantum state
% in the respective trigonometric form
% as as sum is equal to one
if simplify( abs(alpha_trig)^2 + ...
             abs(beta_trig)^2 ) == 1.0
    
    % Print information about
    % the linear and probabilistic combination of
    % magnitudes for the complex coefficients of
    % a two-dimensional quantum state
    % in the respective trigonometric form
    % as sum being equal to one
    fprintf([' (|alpha_trig|)^2 + (|beta_trig|)^2 =\n' ...
             ' = (|cos(theta/2)|)^2 + ' ...
             '(|sin(theta / 2) * exp(phi * 1j)|)^2 = 1.0\n']);
    
    % Print a blank line
    fprintf('\n');

    % Print more information about
    % the linear and probabilistic combination of
    % magnitudes for the complex coefficients of
    % a two-dimensional quantum state
    % in the respective trigonometric form
    % as sum being equal to one
    fprintf(['The condition of the linear and probabilistic\n' ...
             'combination of magnitudes for the complex\n' ...
             'coefficients of a two-dimensional quantum state |psi>\n' ...
             'in the respective trigonometric form as sum\n' ...
             'being equal to one holds']);
    
    % Print a blank line
    fprintf('\n');

end


% Print a blank line
fprintf('\n');


% Compute the radius for
% the cartesian coordinates
% from the spherical coordinates of
% the two-dimensional quantum state |psi>
% in its trigonometric form
radius = simplify( norm(ket_psi_trig) );

% Compute the cartesian coordinates
% from the spherical coordinates of
% the two-dimensional quantum state |psi>
% in its trigonometric form
[ cart_coord_x,...
  cart_coord_y,...
  cart_coord_z] = sph2cart(theta, ...
                           phi, ...
                           radius);

% Print the information of
% the cartesian coordinates computed
% from the spherical coordinates of
% the two-dimensional quantum state |psi>
% in its trigonometric form
fprintf(['The Cartesian Coordinates (x, y, z) of ' ...
         'the quantum state |psi> are:\n']);
fprintf(' x = %s\n', cart_coord_x);
fprintf(' y = %s\n', cart_coord_y);
fprintf(' z = %s\n', cart_coord_z);

% Print a blank line
fprintf('\n');