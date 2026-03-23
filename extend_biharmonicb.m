function [F_e] = extend_biharmonicb(ppp,F,dx,dy,lbrt)
% Extend values using biharmonic function.

% Default dy.
if nargin < 7
    dy = dx;
end

% Size variables.
[m,n] = size(ppp);

p = ones(m,n);

p(isnan(F)) = 0;

p = logical(p);

% I index.
I = reshape(1:m*n,m,n);

% Pad I.
k = 2;
I_pad = pad(I,k,lbrt);

% Indexing vectors.
ivec = (1:m)+k;
jvec = (1:n)+k;

% J indices.
%   Current node.
Jp = I_pad(ivec  ,jvec  );
%   First neighbours.
Je = I_pad(ivec+1,jvec  );
Jw = I_pad(ivec-1,jvec  );
Jn = I_pad(ivec  ,jvec+1);
Js = I_pad(ivec  ,jvec-1);
%   Second neighbours.
Jee = I_pad(ivec+2,jvec  );
Jww = I_pad(ivec-2,jvec  );
Jnn = I_pad(ivec  ,jvec+2);
Jss = I_pad(ivec  ,jvec-2);
%   Corner neighbours.
Jne = I_pad(ivec+1,jvec+1);
Jnw = I_pad(ivec-1,jvec+1);
Jse = I_pad(ivec+1,jvec-1);
Jsw = I_pad(ivec-1,jvec-1);

% Matrix values.
dxdy = dx / dy;
%   Current node.
Vp = (6 + 6*dxdy^4 + 8*dxdy^2) * ones(m,n);
%   First neighbours.
Ve = (- 4 - 4*dxdy^2) * ones(m,n);
Vw = (- 4 - 4*dxdy^2) * ones(m,n);
Vn = (- 4*dxdy^4 - 4*dxdy^2) * ones(m,n);
Vs = (- 4*dxdy^4 - 4*dxdy^2) * ones(m,n);
%   Second neighbours.
Vee = (1) * ones(m,n);
Vww = (1) * ones(m,n);
Vnn = (1*dxdy^4) * ones(m,n);
Vss = (1*dxdy^4) * ones(m,n);
%   Corner neighbours.
Vne = (2*dxdy^2) * ones(m,n);
Vnw = (2*dxdy^2) * ones(m,n);
Vse = (2*dxdy^2) * ones(m,n);
Vsw = (2*dxdy^2) * ones(m,n);

% Replace stencil for predetermined values.
%   Current node.
Vp(p) = 1;
%   First neighbours.
Ve(p) = 0;
Vw(p) = 0;
Vn(p) = 0;
Vs(p) = 0;
%   Second neighbours.
Vee(p) = 0;
Vww(p) = 0;
Vnn(p) = 0;
Vss(p) = 0;
%   Corner neighbours.
Vne(p) = 0;
Vnw(p) = 0;
Vse(p) = 0;
Vsw(p) = 0;

% Form complete indexing and value matrices.
I_mat = [I I I I I I I I I I I I I];
J_mat = [Jp Je Jw Jn Js Jee Jww Jnn Jss Jne Jnw Jse Jsw];
V_mat = [Vp Ve Vw Vn Vs Vee Vww Vnn Vss Vne Vnw Vse Vsw];

% Form sparse matrix.
A = sparse(I_mat(:),J_mat(:),V_mat(:));

% Factorise matrix.
[L,U,P,Q] = lu(A);

b = zeros(m,n);
b(p) = F(p);
b = b(:);
x = Q * (U \ (L \ (P * b)));
F_e = reshape(x,m,n);
























