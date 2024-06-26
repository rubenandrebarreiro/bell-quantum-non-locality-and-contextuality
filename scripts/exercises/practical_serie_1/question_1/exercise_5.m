% Clear Command Window
clc;


% Create the symbolic numbers
% for the complex coefficients of
% a two-dimensional quantum state
% (i.e., alpha and beta)
syms alpha_0 alpha_1 beta_0 beta_1;


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


% Create the ket vector for
% the two-dimensional quantum state
% denoted as |psi> = ([alpha, beta])^T
ket_psi = [alpha beta].';


% Create the Pauli I (sigma_i) Matrix
pauli_sigma_i = full(Pauli('I'));

% Create the Pauli X (sigma_x) Matrix
pauli_sigma_x = full(Pauli('X'));

% Create the Pauli Y (sigma_y) Matrix
pauli_sigma_y = full(Pauli('Y'));

% Create the Pauli Z (sigma_z) Matrix
pauli_sigma_z = full(Pauli('Z'));


% Create the multiplication between
% the Pauli I (sigma_i) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
pauli_sigma_i_mult_ket_psi = ...
    pauli_sigma_i * ket_psi;

% Create the multiplication between
% the Pauli X (sigma_x) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
pauli_sigma_x_mult_ket_psi = ...
    pauli_sigma_x * ket_psi;

% Create the multiplication between
% the Pauli Y (sigma_y) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
pauli_sigma_y_mult_ket_psi = ...
    pauli_sigma_y * ket_psi;

% Create the multiplication between
% the Pauli Z (sigma_z) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
pauli_sigma_z_mult_ket_psi = ...
    pauli_sigma_z * ket_psi;


% Print a blank line 
fprintf('\n');

% Print the content of the ket vector
% for the two-dimensional quantum state
% denoted as |psi> = ([alpha, beta])^T
fprintf('|psi> =\n  ');
disp(ket_psi);

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


% Print the content of the multiplication
% between the Pauli I (sigma_i) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
fprintf('(sigma_i) x |psi> =\n  ');
disp(pauli_sigma_i_mult_ket_psi);

% Print the content of the multiplication
% between the Pauli X (sigma_x) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
fprintf('(sigma_x) x |psi> =\n');
disp(pauli_sigma_x_mult_ket_psi);

% Print the content of the multiplication
% between the Pauli Y (sigma_y) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
fprintf('(sigma_y) x |psi> =\n');
disp(pauli_sigma_y_mult_ket_psi);

% Print the content of the multiplication
% between the Pauli Z (sigma_z) Matrix and
% the ket vector for the two-dimensional
% quantum state denoted as
% |psi> = ([alpha, beta])^T
fprintf('(sigma_z) x |psi> =\n');
disp(pauli_sigma_z_mult_ket_psi);
