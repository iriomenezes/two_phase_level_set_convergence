function [dP] = upwindb(phi,F1,F2,dx,dy,lbrt,c)
% Computes upwind spatial derivatives.

% Size variables.
[m,n] = size(phi);

k = 2;

% Indexing vectors.
i_ind = (1:m)+k;
j_ind = (1:n)+k;

% Pad phi.
phi_pad = pad(phi,k,lbrt);

% Diagonal spatial increment.
dez = sqrt(dx^2 + dy^2);

% xy gradient.
% Derivatives.
phi_x_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind+1,j_ind  ) - phi_pad(i_ind+2,j_ind  )) / (2*dx);
phi_x_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind-1,j_ind  ) + phi_pad(i_ind-2,j_ind  )) / (2*dx);
phi_y_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind  ,j_ind+1) - phi_pad(i_ind  ,j_ind+2)) / (2*dy);
phi_y_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind  ,j_ind-1) + phi_pad(i_ind  ,j_ind-2)) / (2*dy);

% Godunov's method to square derivatives for V>0.
phi_x_2_p = max(max(phi_x_m,0).^2,min(phi_x_p,0).^2);
phi_y_2_p = max(max(phi_y_m,0).^2,min(phi_y_p,0).^2);
% Godunov's method to square derivatives for V<0.
phi_x_2_m = max(min(phi_x_m,0).^2,max(phi_x_p,0).^2);
phi_y_2_m = max(min(phi_y_m,0).^2,max(phi_y_p,0).^2);

grad_xy_p = sqrt(phi_x_2_p + phi_y_2_p);
grad_xy_m = sqrt(phi_x_2_m + phi_y_2_m);

% Form upwind approximation.

F1grad_xy = max(F1,0).*grad_xy_p + min(F1,0).*grad_xy_m;


% ez gradient.
% Derivatives.
phi_e_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind+1,j_ind+1  ) - phi_pad(i_ind+2,j_ind+2  )) / (2*dez);
phi_e_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind-1,j_ind-1  ) + phi_pad(i_ind-2,j_ind-2  )) / (2*dez);
phi_z_p = (- 3*phi_pad(i_ind,j_ind) + 4*phi_pad(i_ind+1  ,j_ind-1) - phi_pad(i_ind+2  ,j_ind-2)) / (2*dez);
phi_z_m = (  3*phi_pad(i_ind,j_ind) - 4*phi_pad(i_ind-1  ,j_ind+1) + phi_pad(i_ind-2  ,j_ind+2)) / (2*dez);

% Godunov's method to square derivatives for V>0.
phi_e_2_p = max(max(phi_e_m,0).^2,min(phi_e_p,0).^2);
phi_z_2_p = max(max(phi_z_m,0).^2,min(phi_z_p,0).^2);
% Godunov's method to square derivatives for V<0.
phi_e_2_m = max(min(phi_e_m,0).^2,max(phi_e_p,0).^2);
phi_z_2_m = max(min(phi_z_m,0).^2,max(phi_z_p,0).^2);

grad_ez_p = sqrt(phi_e_2_p + phi_z_2_p);
grad_ez_m = sqrt(phi_e_2_m + phi_z_2_m);

% Form upwind approximation.

F2grad_ez = max(F2,0).*grad_ez_p + min(F2,0).*grad_ez_m;
% dP = 0.5*(F1grad_xy + F2grad_ez);
dP = 0.5.*(max(F2,0).*(grad_ez_p-c)+min(F2,0).*(grad_ez_m-c) + ...
    max(F1,0).*(grad_xy_p-c)+min(F1,0).*(grad_xy_m-c));






























