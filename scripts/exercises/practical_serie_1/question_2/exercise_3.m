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
assume(theta < 2* pi);

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
beta_trig = sin(theta / 2) * exp(1j * phi);


% Create the ket vector for
% the two-dimensional quantum state
% denoted as |psi> = ([alpha, beta])^T
ket_psi = [alpha beta].';

% Compute the two-dimensional quantum state
% in its trigonometric form in function of
% the respective theta and phi angles
ket_psi_trig = [alpha_trig beta_trig].';

% Create the bra vector for
% the two-dimensional quantum state
% denoted as <psi| = [alpha*, beta*]
bra_psi = conj(transpose(ket_psi));

% Create the bra vector for
% the two-dimensional quantum state
% in its trigonometric form in function of
% the respective theta and phi angles
% denoted as <psi| = [alpha*, beta*]
bra_psi_trig = conj(transpose(ket_psi_trig));


% Create the Pauli I (sigma_i) Matrix
pauli_sigma_i = full(Pauli('I'));

% Create the Pauli X (sigma_x) Matrix
pauli_sigma_x = full(Pauli('X'));

% Create the Pauli Y (sigma_y) Matrix
pauli_sigma_y = full(Pauli('Y'));

% Create the Pauli Z (sigma_z) Matrix
pauli_sigma_z = full(Pauli('Z'));


% Create the multiplication between
% the bra vector for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli X (sigma_x) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
m_x_bra_psi_trig_mult_pauli_sigma_x_mult_ket_psi_trig = ...
    simplify(bra_psi_trig * pauli_sigma_x * ket_psi_trig);

% Create the multiplication between
% the bra vector for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli Y (sigma_y) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
m_y_bra_psi_trig_mult_pauli_sigma_y_mult_ket_psi_trig = ...
    simplify(bra_psi_trig * pauli_sigma_y * ket_psi_trig);

% Create the multiplication between
% the bra vector for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli Z (sigma_z) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
m_z_bra_psi_trig_mult_pauli_sigma_z_mult_ket_psi_trig = ...
    simplify(bra_psi_trig * pauli_sigma_z * ket_psi_trig);

% Create the vector with the multiple
% multiplication operations between
% the bra vectors for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli X, Y, Z (sigma_x, sigma_y, sigma_z) Matrices
% and the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
m_vec = { m_x_bra_psi_trig_mult_pauli_sigma_x_mult_ket_psi_trig, ...
          m_y_bra_psi_trig_mult_pauli_sigma_y_mult_ket_psi_trig, ...
          m_z_bra_psi_trig_mult_pauli_sigma_z_mult_ket_psi_trig };


% Print a blank line 
fprintf('\n');

% Print the content of the ket vector
% for the two-dimensional quantum state
% denoted as |psi> = ([alpha, beta])^T
fprintf('|psi> =\n  ');
disp(ket_psi);

% Print a blank line 
fprintf('\n');

% Print the content of the ket vector
% for the two-dimensional quantum state
% denoted as <psi| = [alpha*, beta*]
fprintf('<psi| =\n  ');
disp(bra_psi);

% Print a blank line 
fprintf('\n');


% Print the content of
% the Pauli I (sigma_i) Matrix
fprintf('Pauli I Matrix (sigma_i) =\n');
disp(pauli_sigma_i);

% Print the content of
% the Pauli X (sigma_x) Matrix
fprintf('Pauli X Matrix (sigma_x) =\n');
disp(pauli_sigma_x);

% Print the content of
% the Pauli Y (sigma_y) Matrix
fprintf('Pauli Y Matrix (sigma_y) =\n');
disp(pauli_sigma_y);

% Print the content of
% the Pauli Z (sigma_z) Matrix
fprintf('Pauli Z Matrix (sigma_z) =\n');
disp(pauli_sigma_z);


% Print the content of the multiplication between
% the bra vector for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli X (sigma_x) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
fprintf('m_x = <psi| x (sigma_x) x |psi> =\n');
disp(m_x_bra_psi_trig_mult_pauli_sigma_x_mult_ket_psi_trig);

% Print the content of the multiplication between
% the bra vector for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli Y (sigma_y) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
fprintf('m_y = <psi| x (sigma_y) x |psi> =\n');
disp(m_y_bra_psi_trig_mult_pauli_sigma_y_mult_ket_psi_trig);

% Print the content of the multiplication between
% the bra vector for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli Z (sigma_z) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
fprintf('m_z = <psi| x (sigma_z) x |psi> =\n');
disp(m_z_bra_psi_trig_mult_pauli_sigma_z_mult_ket_psi_trig);


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
cart_coord_x = radius * sin(theta) * cos(phi);
cart_coord_y = radius * sin(theta) * sin(phi);
cart_coord_z = radius * cos(theta);


% If the cartesian coordinates computed
% from the spherical coordinates of
% the two-dimensional quantum state |psi>
% in its trigonometric form are
% equal to the vector with the multiple
% multiplication operations between
% the bra vectors for the two-dimensional
% quantum state denoted as |psi> = [alpha*, beta*],
% the Pauli X, Y, Z (sigma_x, sigma_y, sigma_z) Matrices
% and the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
if isequal(cart_coord_x, ...
           m_x_bra_psi_trig_mult_pauli_sigma_x_mult_ket_psi_trig) ...
   && isequal(cart_coord_y, ...
              m_y_bra_psi_trig_mult_pauli_sigma_y_mult_ket_psi_trig) ...
   && isequal(cart_coord_z, ...
              m_z_bra_psi_trig_mult_pauli_sigma_z_mult_ket_psi_trig)

    % Print the information of
    % the cartesian coordinates computed
    % from the spherical coordinates of
    % the two-dimensional quantum state |psi>
    % in its trigonometric form
    fprintf(['The Cartesian Coordinates (x, y, z) of ' ...
             'the quantum state |psi> are:\n']);
    fprintf(' x = %s\n', simplify(cart_coord_x));
    fprintf(' y = %s\n', simplify(cart_coord_y));
    fprintf(' z = %s\n', simplify(cart_coord_z));
    
    % Print a blank line
    fprintf('\n');
   
end