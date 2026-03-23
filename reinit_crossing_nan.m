function phi_sd = reinit_crossing_nan(phi0,dx,dy,lbrt)
% Reinitialise implicit function using crossing times.

% Default dy.
if nargin < 3
    dy = dx;
end

% Size variables.
[m,n] = size(phi0);

% Intialise phi.
phi_sd = zeros(m,n);

% Iterations in small band.
K1 = 10;

% Iterations before checking.
K2 = 10;

% Outside.
t = 0;
phi = phi0;
% Small band.
dt = min(dx,dy)/ 4;
% Initialise storage.
PHI = zeros(m,n,K1+1);
PHI(:,:,1) = phi;
T = zeros(1,K1+1);
T(1) = t;
% Iterate.
for k1 = 1:K1
    % Pad phi.
    k = 2;
    phi_pad = pad(phi,k,lbrt);
    % Compute spatial derivatives.
    i_ind = (1:m)+k;
    j_ind = (1:n)+k;
    phi_x_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind+1,j_ind  ) - phi_pad(i_ind+2,j_ind  )) / (2*dx);
    phi_x_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind-1,j_ind  ) + phi_pad(i_ind-2,j_ind  )) / (2*dx);
    phi_y_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind  ,j_ind+1) - phi_pad(i_ind  ,j_ind+2)) / (2*dy);
    phi_y_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind  ,j_ind-1) + phi_pad(i_ind  ,j_ind-2)) / (2*dy);
    % Godunov's method to square derivatives.
    phi_x_2 = max(max(phi_x_m,0).^2,min(phi_x_p,0).^2);
    phi_y_2 = max(max(phi_y_m,0).^2,min(phi_y_p,0).^2);
    % Compute grad.
    grad_phi = sqrt(phi_x_2 + phi_y_2);
    grad_phi(grad_phi==0) = 1;
    % Update phi.
    phi = phi - dt * (phi_x_2 + phi_y_2) ./ grad_phi;
    t = t + dt;
    % Store solutions.
    PHI(:,:,k1+1) = phi;
    T(k1+1) = t;
end
% Determine whether interface has crossed nodes.
[I,J] = find((PHI(:,:,1)>0) .* (PHI(:,:,end)<0));
% Add time to signed distance function.
for ii = 1:length(I)
    i = I(ii); j = J(ii);
    kk = sum(PHI(i,j,:)>0);
    phi_sd(i,j) = T(kk) + dt * PHI(i,j,kk) ./ (PHI(i,j,kk) - PHI(i,j,kk+1));
end
% Remaining nodes.
dt = min(dx,dy)/ 2;
while max(max(phi)) >= 0
    % Initialise storage.
    PHI = zeros(m,n,K2+1);
    PHI(:,:,1) = phi;
    T = zeros(1,K2+1);
    T(1) = t;
    % Iterate.
    for k2 = 1:K2
        % Pad phi.
        k = 1;
        phi_pad = pad(phi,k,lbrt);
        % Compute spatial derivatives.
        i_ind = (1:m)+k;
        j_ind = (1:n)+k;
        phi_x_p = (phi_pad(i_ind+1,j_ind  ) - phi_pad(i_ind  ,j_ind  )) / dx;
        phi_x_m = (phi_pad(i_ind  ,j_ind  ) - phi_pad(i_ind-1,j_ind  )) / dx;
        phi_y_p = (phi_pad(i_ind  ,j_ind+1) - phi_pad(i_ind  ,j_ind  )) / dy;
        phi_y_m = (phi_pad(i_ind  ,j_ind  ) - phi_pad(i_ind  ,j_ind-1)) / dy;
        % Godunov's method to square derivatives.
        phi_x_2 = max(max(phi_x_m,0).^2,min(phi_x_p,0).^2);
        phi_y_2 = max(max(phi_y_m,0).^2,min(phi_y_p,0).^2);
        % Compute grad.
        grad_phi = sqrt(phi_x_2 + phi_y_2);
        grad_phi(grad_phi==0) = 1;
        % Update phi.
        phi = phi - dt * (phi_x_2 + phi_y_2) ./ grad_phi;
        t = t + dt;
        % Store solutions.
        PHI(:,:,k2+1) = phi;
        T(k2+1) = t;
    end
    % Determine whether interface has crossed nodes.
    [I,J] = find((PHI(:,:,1)>0) .* (PHI(:,:,end)<0));
    % Add time to signed distance function.
    for ii = 1:length(I)
        i = I(ii); j = J(ii);
        kk = sum(PHI(i,j,:)>0);
        phi_sd(i,j) = T(kk) + dt * PHI(i,j,kk) ./ (PHI(i,j,kk) - PHI(i,j,kk+1));
    end
end

% Inside.
t = 0;
phi = phi0;
% Small band.
dt = min(dx,dy)/ 4;
% Initialise storage.
PHI = zeros(m,n,K1+1);
PHI(:,:,1) = phi;
T = zeros(1,K1+1);
T(1) = t;
% Iterate.
for k1 = 1:K1
    % Pad phi.
    k = 2;
    phi_pad = pad(phi,k,lbrt);
    i_ind = (1:m)+k;
    j_ind = (1:n)+k;
    phi_x_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind+1,j_ind  ) - phi_pad(i_ind+2,j_ind  )) / (2*dx);
    phi_x_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind-1,j_ind  ) + phi_pad(i_ind-2,j_ind  )) / (2*dx);
    phi_y_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind  ,j_ind+1) - phi_pad(i_ind  ,j_ind+2)) / (2*dy);
    phi_y_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind  ,j_ind-1) + phi_pad(i_ind  ,j_ind-2)) / (2*dy);
    % Godunov's method to square derivatives.
    phi_x_2 = max(min(phi_x_m,0).^2,max(phi_x_p,0).^2);
    phi_y_2 = max(min(phi_y_m,0).^2,max(phi_y_p,0).^2);
    % Compute grad.
    grad_phi = sqrt(phi_x_2 + phi_y_2);
    grad_phi(grad_phi==0) = 1;
    % Update phi.
    phi = phi + dt * (phi_x_2 + phi_y_2) ./ grad_phi;
    t = t + dt;
    % Store solutions.
    PHI(:,:,k1+1) = phi;
    T(k1+1) = t;
end
% Determine whether interface has crossed nodes.
[I,J] = find((PHI(:,:,1)<0) .* (PHI(:,:,end)>0));
% Add time to signed distance function.
for ii = 1:length(I)
    i = I(ii); j = J(ii);
    kk = sum(PHI(i,j,:)<0);
    phi_sd(i,j) = - (T(kk) + dt * PHI(i,j,kk) ./ (PHI(i,j,kk) - PHI(i,j,kk+1)));
end
% Remaining nodes.
dt = min(dx,dy)/ 2;
while min(min(phi)) <= 0
    % Initialise storage.
    PHI = zeros(m,n,K2+1);
    PHI(:,:,1) = phi;
    T = zeros(1,K2+1);
    T(1) = t;
    % Iterate.
    for k2 = 1:K2
        % Pad phi.
        k = 2;
        phi_pad = pad(phi,k,lbrt);
        i_ind = (1:m)+k;
        j_ind = (1:n)+k;
        phi_x_p = (phi_pad(i_ind+1,j_ind  ) - phi_pad(i_ind  ,j_ind  )) / dx;
        phi_x_m = (phi_pad(i_ind  ,j_ind  ) - phi_pad(i_ind-1,j_ind  )) / dx;
        phi_y_p = (phi_pad(i_ind  ,j_ind+1) - phi_pad(i_ind  ,j_ind  )) / dy;
        phi_y_m = (phi_pad(i_ind  ,j_ind  ) - phi_pad(i_ind  ,j_ind-1)) / dy;
        % Godunov's method to square derivatives.
        phi_x_2 = max(min(phi_x_m,0).^2,max(phi_x_p,0).^2);
        phi_y_2 = max(min(phi_y_m,0).^2,max(phi_y_p,0).^2);
        % Compute grad.
        grad_phi = sqrt(phi_x_2 + phi_y_2);
        grad_phi(grad_phi==0) = 1;
        % Update phi.
        phi = phi + dt * (phi_x_2 + phi_y_2) ./ grad_phi;
        t = t + dt;
        % Store solutions.
        PHI(:,:,k2+1) = phi;
        T(k2+1) = t;
    end
    % Determine whether interface has crossed nodes.
    [I,J] = find((PHI(:,:,1)<0) .* (PHI(:,:,end)>0));
    % Add time to signed distance function.
    for ii = 1:length(I)
        i = I(ii); j = J(ii);
        kk = sum(PHI(i,j,:)<0);
        phi_sd(i,j) = - (T(kk) + dt * PHI(i,j,kk) ./ (PHI(i,j,kk) - PHI(i,j,kk+1)));
    end
end





























