% Create the Pauli I (sigma_i) Matrix
pauli_sigma_i = full(Pauli('I'));

% Create the Pauli X (sigma_x) Matrix
pauli_sigma_x = full(Pauli('X'));

% Create the Pauli Y (sigma_y) Matrix
pauli_sigma_y = full(Pauli('Y'));

% Create the Pauli Z (sigma_z) Matrix
pauli_sigma_z = full(Pauli('Z'));


% Create the matrix multiplication
% between the Pauli X (sigma_x) and
% the Pauli Y (sigma_y) Matrices
pauli_sigma_x_mult_pauli_sigma_y = ...
    pauli_sigma_x * pauli_sigma_y;

% Create the matrix multiplication
% between the Pauli X (sigma_x) and
% the Pauli Y (sigma_y) Matrices
pauli_sigma_y_mult_pauli_sigma_x = ...
    pauli_sigma_y * pauli_sigma_x;


% Create the scalar multiplication between
% (0 + 1j) and the Pauli Z (sigma_z) Matrix
imag_num_mult_pauli_sigma_z = ...
    1j * pauli_sigma_z;

% Create the scalar multiplication between
% (0 - 1j) and the Pauli Z (sigma_z) Matrix
minus_imag_num_mult_pauli_sigma_z = ...
    -1j * pauli_sigma_z;


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


% Print the content of
% the matrix multiplication
% between the Pauli X (sigma_x) and
% the Pauli Y (sigma_y) Matrices
fprintf('(sigma_x) x (sigma_y) =\n');
disp(pauli_sigma_x_mult_pauli_sigma_y);

% Print the content of
% the matrix multiplication
% between the Pauli Y (sigma_y) and
% the Pauli X (sigma_x) Matrices
fprintf('(sigma_y) x (sigma_x) =\n');
disp(pauli_sigma_y_mult_pauli_sigma_x);


% Print the content of
% the scalar multiplication between
% (0 + 1j) and the Pauli Z (sigma_z) Matrix
fprintf('(0 + 1j) x (sigma_z) =\n');
disp(imag_num_mult_pauli_sigma_z);

% Print the content of
% the scalar multiplication between
% (0 - 1j) and the Pauli Z (sigma_z) Matrix
fprintf('(0 - 1j) x (sigma_z) =\n');
disp(minus_imag_num_mult_pauli_sigma_z);


% Print a blank line 
fprintf('\n');


% If the matrix multiplication
% between the Pauli X (sigma_x) and
% the Pauli Y (sigma_y) Matrices
% is equal to the scalar multiplication between
% (0 + 1j) and the Pauli Z (sigma_z) Matrix
if isequal(pauli_sigma_x_mult_pauli_sigma_y, ...
           imag_num_mult_pauli_sigma_z)
    
    % Print information about
    % the matrix multiplication
    % between the Pauli X (sigma_x)
    % and the Pauli Y (sigma_y)
    % Matrices is equal to the scalar
    % multiplication between (0 + 1j) and
    % the Pauli Z (sigma_z) Matrix
    fprintf(['(0 + 1j) x (sigma_z) = ' ...
             '(sigma_x) x (sigma_y)']);
    
    % Print a blank line 
    fprintf('\n');

end


% Print a blank line 
fprintf('\n');


% If the matrix multiplication
% between the Pauli Y (sigma_y) and
% the Pauli X (sigma_x) Matrices
% is equal to the scalar multiplication between
% (0 - 1j) and the Pauli Z (sigma_z) Matrix
if isequal(pauli_sigma_y_mult_pauli_sigma_x, ...
           minus_imag_num_mult_pauli_sigma_z)
    
    % Print information about
    % the matrix multiplication
    % between the Pauli Y (sigma_y)
    % and the Pauli X (sigma_x)
    % Matrices is equal to the scalar
    % multiplication between (0 - 1j) and
    % the Pauli Z (sigma_z) Matrix
    fprintf(['(0 - 1j) x (sigma_z) = ' ...
             '(sigma_y) x (sigma_x)']);
    
    % Print a blank line 
    fprintf('\n');

end


% Print a blank line 
fprintf('\n');
