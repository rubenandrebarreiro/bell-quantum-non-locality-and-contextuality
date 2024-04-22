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

% Create the assumptions regarding
% the trigonometric angles theta and phi
% being real numbers
assume(theta, 'real');
assume(phi, 'real');


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

% Create vector of the Pauli X, Y, and Z
% (sigma_x, sigma_y, sigma_z) Matrices
pauli_sigmas_vec = { pauli_sigma_x, ...
                     pauli_sigma_y, ...
                     pauli_sigma_z };


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
% the Pauli X, Y, and Z (sigma_x, sigma_y, and sigma_z) Matrices
% and the ket vector for the two-dimensional
% quantum state denoted as |psi> = ([alpha, beta])^T
m_vec = { m_x_bra_psi_trig_mult_pauli_sigma_x_mult_ket_psi_trig, ...
          m_y_bra_psi_trig_mult_pauli_sigma_y_mult_ket_psi_trig, ...
          m_z_bra_psi_trig_mult_pauli_sigma_z_mult_ket_psi_trig };


% Create a 2x2 Identity Matrix
identity_matrix = eye(2);

% Create the P^{vec{m}}_{+} Matrix as being equal to
% ( ( I + vec(m) . vec(sigma) ) / 2 )
p_vec_m_plus = 1 / 2 * ( identity_matrix + ...
                         ( ( m_x_bra_psi_trig_mult_pauli_sigma_x_mult_ket_psi_trig * ... 
                             pauli_sigma_x ) + ...
                           ( m_y_bra_psi_trig_mult_pauli_sigma_y_mult_ket_psi_trig * ...
                             pauli_sigma_y ) + ...
                           ( m_z_bra_psi_trig_mult_pauli_sigma_z_mult_ket_psi_trig * ...
                             pauli_sigma_z ) ) ...
                        );

% Simplify the P^{vec{m}}_{+} Matrix
% created previously
p_vec_m_plus = simplify(p_vec_m_plus);

% Perform the substitution of ( ( 1 + cos(2 * theta) ) / 2 )
% by cos(theta)^2 on the P^{vec{m}}_{+} Matrix
% computed previously
p_vec_m_plus = subs( p_vec_m_plus, ...
                     ( ( 1 + cos(theta) ) / 2 ), ...
                     ( ( cos(theta / 2) )^2 ) );

% Perform the substitution of ( ( 1 + cos(2 * theta) ) / 2 )
% by sin(theta)^2 on the P^{vec{m}}_{+} Matrix
% computed previously
p_vec_m_plus = subs( p_vec_m_plus, ...
                     ( ( 1 - cos(theta) ) / 2 ), ...
                     ( ( sin(theta / 2) )^2 ) );

% Perform the substitution of ( cos(phi) + sin(phi) * 1j )
% by e^(phi * 1j) on the P^{vec{m}}_{+} Matrix
% computed previously
p_vec_m_plus = subs( p_vec_m_plus, ...
                     ( cos(phi) + sin(phi) * 1j ), ...
                     ( exp(phi * 1j) ) );

% Perform the substitution of ( cos(phi) - sin(phi) * 1j )
% by ( e^(-phi * 1j) ) on the P^{vec{m}}_{+} Matrix
% computed previously
p_vec_m_plus = subs( p_vec_m_plus, ...
                     ( cos(phi) - sin(phi) * 1j ), ...
                     ( exp(-phi * 1j) ) );

% Perform the substitution of ( sin(2 * theta) )
% by ( 2 * sin(theta) * cos(theta) ) on
% the P^{vec{m}}_{+} Matrix
% computed previously
p_vec_m_plus = subs( p_vec_m_plus, ...
                     ( sin(theta) / 2 ), ...
                     ( sin(theta / 2) * cos(theta / 2) ) );

% Perform the substitution of ( e^(phi * 1j) )
% by ( e^(-phi * ij) ) on the P^{vec{m}}_{+} Matrix
% computed previously
p_vec_m_plus = subs( p_vec_m_plus, ...
                     exp(phi * 1j), exp(-phi * 1i) );


% Compute the density
% matrix rho = |psi><psi| for
% the two-dimensional quantum state psi
rho_density_matrix = ket_psi_trig * bra_psi_trig;

% Perform the substitution of ( e^(phi * 1j) )
% by ( e^(-phi * ij) ) on the density matrix
% rho = |psi><psi| for the two-dimensional
% quantum state psi computed previously
rho_density_matrix = subs( rho_density_matrix, ...
                           exp(phi * 1j), exp(-phi * 1i) );


% Print a blank line 
fprintf('\n');

% Print the content of the ket vector
% for the two-dimensional quantum state
% denoted as |psi> = ([alpha, beta])^T
fprintf('|psi> =\n  ');
disp(ket_psi_trig);

% Print a blank line 
fprintf('\n');

% Print the content of the ket vector
% for the two-dimensional quantum state
% denoted as <psi| = [alpha*, beta*]
fprintf('<psi| =\n  ');
disp(bra_psi_trig);

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


% If the P^{vec{m}}_{+} Matrix
% created previously is equal to
% the density matrix rho = |psi><psi|
% for the two-dimensional quantum state psi
if isequal( simplify(p_vec_m_plus), ...
            simplify(rho_density_matrix) )
    
    % Print the content of the the P^{vec{m}}_{+} Matrix
    fprintf('P^{vec{m}}_{+} =\n');
    fprintf('  [%s   %s]\n', ...
            p_vec_m_plus(1), ...
            p_vec_m_plus(2));
    fprintf('  [ %s     %s ]\n', ...
            p_vec_m_plus(3), ...
            p_vec_m_plus(4));
    
    % Print a blank line 
    fprintf('\n');
    
    
    % Print the content of the density matrix
    % rho = |psi><psi| for the two-dimensional
    % quantum state psi
    fprintf('rho = |psi><psi| =\n');
    fprintf('  [%s   %s]\n', ...
            rho_density_matrix(1), ...
            rho_density_matrix(2));
    fprintf('  [ %s     %s ]\n', ...
            rho_density_matrix(3), ...
            rho_density_matrix(4));
    
    % Print a blank line 
    fprintf('\n');
    
    % Print some information about
    % the P^{vec{m}}_{+} Matrix being equal to
    % the density matrix rho = |psi><psi|
    % for the two-dimensional quantum state psi
    fprintf(['The following equality holds:\n' ...
             '=> P^{vec{m}}_{+} = ' ...
             '( ( I + vec(m) . vec(sigma) ) / 2 ) = rho\n' ...
             '   with vec(m) = <psi|vec(sigma)|psi>\n']);
    
    % Print a blank line 
    fprintf('\n');

end