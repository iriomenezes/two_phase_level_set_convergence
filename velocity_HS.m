function [dP,F] = velocity_HSb(phi,T1,T2,dx,dy,lbrt,eta1,eta2)
% Computes the velocity field for the Stefan problem.

% Size variables.
[m,n] = size(phi);

% Pad phi, T.
k = 3;

phi_pad = pad(phi,k,'eeee');
T1_pad = pad(T1,k,'eeee');
T2_pad = pad(T2,k,'eeee');

% Indexing vectors.
ivec = (1:m)+k;
jvec = (1:n)+k;

% Diagonal spatial increment.
dez = sqrt(dx^2 + dy^2);

% Initial derivative stencils.
% u1_x = (T1_pad(ivec+1,jvec  ) - T1_pad(ivec-1,jvec  )) / (2*dx);
% u1_y = (T1_pad(ivec  ,jvec+1) - T1_pad(ivec  ,jvec-1)) / (2*dy);
% u1_e = (T1_pad(ivec+1,jvec+1) - T1_pad(ivec-1,jvec-1)) / (2*dez);
% u1_z = (T1_pad(ivec+1,jvec-1) - T1_pad(ivec-1,jvec+1)) / (2*dez);

u2_x = (T2_pad(ivec+1,jvec  ) - T2_pad(ivec-1,jvec  )) / (2*dx);
u2_y = (T2_pad(ivec  ,jvec+1) - T2_pad(ivec  ,jvec-1)) / (2*dy);
u2_e = (T2_pad(ivec+1,jvec+1) - T2_pad(ivec-1,jvec-1)) / (2*dez);
u2_z = (T2_pad(ivec+1,jvec-1) - T2_pad(ivec-1,jvec+1)) / (2*dez);

% Central derivatives of phi

phi_xc = (phi_pad(ivec+1,jvec  ) - phi_pad(ivec-1,jvec  )) / (2*dx);
phi_yc = (phi_pad(ivec  ,jvec+1) - phi_pad(ivec  ,jvec-1)) / (2*dy);
phi_ec = (phi_pad(ivec+1,jvec+1) - phi_pad(ivec-1,jvec-1)) / (2*dez);
phi_zc = (phi_pad(ivec+1,jvec-1) - phi_pad(ivec-1,jvec+1)) / (2*dez);

grad_xyc = sqrt(phi_xc.^2 + phi_yc.^2 + dx*dy);
grad_ezc = sqrt(phi_ec.^2 + phi_zc.^2 + dx*dy);

% F1_xy = (-phi_xc.*u1_x - phi_yc.*u1_y)./grad_xyc;
% F1_ez = (-phi_ec.*u1_e - phi_zc.*u1_z)./grad_ezc;

F2_xy = (-phi_xc.*u2_x - phi_yc.*u2_y)./grad_xyc;
F2_ez = (-phi_ec.*u2_e - phi_zc.*u2_z)./grad_ezc;

% Biharmonic extension.
% [F1_xy_e] = extend_biharmonicb(phi<0,F1_xy,dx,dy,lbrt);
% [F1_ez_e] = extend_biharmonicb(phi<0,F1_ez,dx,dy,lbrt);
[F2_xy_e] = extend_biharmonicb(phi>0,F2_xy,dx,dy,lbrt);
[F2_ez_e] = extend_biharmonicb(phi>0,F2_ez,dx,dy,lbrt);

% ENO spatial derivatives.
% [dP1] = upwindb(phi,F1_xy_e,F1_ez_e,dx,dy,lbrt,0);
[dP2] = upwindb(phi,F2_xy_e,F2_ez_e,dx,dy,lbrt,0);

% F = max(max(max(abs(F2_xy(:))));
F = max(abs(F2_xy(:)));

% dP = (eta1*dP1 + eta2*dP2)/(eta1+eta2);

% F = F2_xy_e;
dP = dP2;

end





























