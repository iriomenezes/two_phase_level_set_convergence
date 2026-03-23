function [T1,T2] = pressure_pol(phi,sigma,X,Y,domain,eta1,eta2,Q)
% Recreating the cartesian grid
warning('off', 'MATLAB:griddata:DuplicateDataPoints')
x = X(1,:);
y = Y(:,1);
[X, Y] = meshgrid(x, y);

% Convert to polars
Nodes = length(x);
t1    = -pi+2*pi*rand();
theta = linspace(t1, t1 + 2.*pi, Nodes+1)';
theta(end) = [];
dth = theta(2) - theta(1);
%r = (dth : (theta(2) - theta(1)) : domain)';%linspace(0,domain,Nodes);
r = linspace(0,domain,Nodes);
dr = r(2) - r(1);
nth = length(theta);
[R, Theta] = meshgrid(r, theta);
[RR, TT]   = pol2cart(Theta, R); % Unstructed coordiantes

phi_pol = interp2(X, Y, phi, RR, TT, 'cubic');


weight = zeros(nth,nth);
tt = theta(:,ones(1,nth));
tttt = (theta');
tt2 = tttt(ones(nth,1),:);
for kk = 1:10
    weight = weight + kk.*cos(kk*(tt2 - tt));
end

% Injection Source term
rroo = 0.05;
source = -Q*eta1*(1 + cos(pi*R/rroo))/rroo^2;
source(R>rroo) = 0;
b = source;


%Solve laplace`s equation

[m,n] = size(phi_pol);

% I index.
I = reshape(1:m*n,m,n);
% % 

I_farf = I(:,end)*(ones(1,m));
J_farf = (ones(m,1))*I(:,end)';
V_farf = -2*(dth/dr/pi/R(end,end))*weight;

% Pad I, phi.
k = 2;
I_pad = pad_pol(I,k,'pnpn');
phi_pad = pad_pol(phi_pol,k,'pepe');
R_pad = pad_pol(R,k,'pzpe');


% Indexing vectors.
ivec = (1:m)+k;
jvec = (1:n)+k;

% J indices.
%   Current node.
Jp = I_pad(ivec  ,jvec  );
%   Neighbours.
Je = I_pad(ivec+1,jvec  );
Jw = I_pad(ivec-1,jvec  );
Jn = I_pad(ivec  ,jvec+1);
Js = I_pad(ivec  ,jvec-1);

% RHS values.
% b = source;
% b = zeros(m,n);

% Matrix values.
%   Current node.
Vp = -(2/dr^2) * ones(m,n) - 2 ./ (R.^2) / dth^2;
%   Neighbours.
Ve = 1 ./ R.^2 / dth^2;
Vw = 1 ./ R.^2 / dth^2;
Vn = (ones(m,n) + dr ./ (2.*R)) / dr^2;
Vs = (ones(m,n) - dr ./ (2.*R)) / dr^2;


% Check for nodes near the interface.
for ip = ivec
    for jp = jvec
        i = ip - k; j = jp - k;
        % Extract node.
        phip = phi_pad(ip,jp);
        rp = R_pad(ip,jp);
        % Check if on the interface.
%         if phip == 0
%             Tp = front(phi_pad(ip-1:ip+1,jp-1:jp+1),rp,'p');
%             Vp(i,j) = 1;
%             Vw(i,j) = 0;
%             Ve(i,j) = 0;
%             Vn(i,j) = 0;
%             Vs(i,j) = 0;
%             b(i,j) = Tp;
%         else
            % theta direction.

            % Extract neighbours.
            phie = phi_pad(ip+1,jp);
            phiw = phi_pad(ip-1,jp);
            he=dth;hw=dth;hs=dr;hn=dr;
            % Check if near the interface.
            if phie*phip <= 0 || phiw*phip <= 0
                % Check east side.
                if phie * phip > 0
                    he = dth;
                else
                    Te = front(phi_pad(ip-1:ip+2,jp-1:jp+1),rp,'x');
                    he = phip / (phip - phie) * dth;
                end
                % Check west side.
                if phiw * phip > 0
                    hw = dth;
                else
                    Tw = front(phi_pad(ip-2:ip+1,jp-1:jp+1),rp,'x');
                    hw = phip / (phip - phiw) * dth;
                end

                % Adjust matrix values.

                if phie * phip >= 0 %NO cross to the east 
                    Vp(i,j) = Vp(i,j) + 1/rp^2/dth^2 - 2/(hw+he)/he/rp^2;
                    Ve(i,j) = 2/he/(hw+he)/rp^2;

                else %if we cross to the east
                    if phip >= 0 %and we are in fluid 2
                        dive = eta2*he + eta1*(dth - he);
                        Vp(i,j) = Vp(i,j) + 1/rp^2/dth^2 - 2*eta2/rp^2/(he+hw)/dive;
                        Ve(i,j) = 2*eta2/rp^2/(he+hw)/dive;
                        b(i,j) = b(i,j) - 2*eta2/rp^2/(he+hw)/dive * Te;

                    else %and we are in fluid 1
                        dive = eta1*he + eta2*(dth - he);
                        Vp(i,j) = Vp(i,j) + 1/rp^2/dth^2 - 2*eta1/rp^2/(he+hw)/dive;
                        Ve(i,j) = 2*eta1/rp^2/(he+hw)/dive;
                        b(i,j) = b(i,j) + 2*eta1/rp^2/(he+hw)/dive * Te;
                    end
                end

                if phiw * phip >= 0 %No cross to the west
                    Vp(i,j) = Vp(i,j) + 1/rp^2/dth^2 - 2/(hw+he)/hw/rp^2;
                    Vw(i,j) = 2/hw/(he+hw)/rp^2;

                else %if we cross to the west
                    if phip>=0 %and we are in fluid 2
                        divw = eta2*hw + eta1*(dth - hw);
                        Vp(i,j) = Vp(i,j) + 1/rp^2/dth^2 - 2*eta2/rp^2/(he+hw)/divw;
                        Vw(i,j) = 2*eta2/rp^2/(he+hw)/divw;
                        b(i,j) = b(i,j) - 2*eta2/rp^2/(he+hw)/divw * Tw;

                    else %and we are in fluid 1
                        divw = eta1*hw + eta2*(dth - hw);
                        Vp(i,j) = Vp(i,j) + 1/rp^2/dth^2 - 2*eta1/rp^2/(he+hw)/divw;
                        Vw(i,j) = 2*eta1/rp^2/(he+hw)/divw;
                        b(i,j) = b(i,j) + 2*eta1/rp^2/(he+hw)/divw * Tw;
                    end
                end
            end


            % r direction.

            % Extract neighbours.
            phin = phi_pad(ip,jp+1);
            phis = phi_pad(ip,jp-1);
            % Check if near the interface.
            if phin*phip <= 0 || phis*phip <= 0
                % Check north side.
                if phin * phip > 0
                    hn = dr;
                else
                    Tn = front(phi_pad(ip-1:ip+1,jp-1:jp+2),R_pad(ip,jp:jp+1),'y');
                    hn = phip / (phip - phin)*dr;
                end
                rI_n = (rp + hn);
                rhat_n = 0.5*(rp + rI_n);

                % Check south side.
                if phis * phip > 0
                    hs = dr;
                else
                    Ts = front(phi_pad(ip-1:ip+1,jp-2:jp+1),R_pad(ip,jp-1:jp),'y');
                    hs = phip / (phip - phis)*dr;
                end
                rI_s = (rp - hs);
                rhat_s = 0.5*(rp + rI_s);

                % Adjust matrix values.
                if phin * phip >= 0 %NO cross to the north
                    Vp(i,j) = Vp(i,j) + 1/dr^2 - 2*(rhat_n/hn)/rp/(hn+hs);
                    Vn(i,j) = 2*rhat_n/rp/hn/(hs+hn);

                else %if we cross to the north
                    if phip >= 0 % and if we are on fluid 2
                        divn = eta2*hn + eta1*(dr - hn);
                        Vp(i,j) = Vp(i,j) + 1/dr^2 - 2*rhat_n*eta2/rp/(hn+hs)/divn;
                        Vn(i,j) = 2*rhat_n*eta2/rp/(hn+hs)/divn;
                        b(i,j) = b(i,j) - 2*rhat_n*eta2/rp/(hn+hs)/divn * Tn; 

                    else % and if we are on fluid 1 
                        divn = eta1*hn + eta2*(dr - hn);
                        Vp(i,j) = Vp(i,j) + 1/dr^2 - 2*rhat_n*eta1/rp/(hn+hs)/divn;
                        Vn(i,j) = 2*rhat_n*eta1/rp/(hn+hs)/divn;
                        b(i,j) = b(i,j) + 2*rhat_n*eta1/rp/(hn+hs)/divn * Tn;
                    end
                end

                if phis * phip >= 0 %NO cross to the south
                    Vp(i,j) = Vp(i,j) + 1/dr^2 - 2*(rhat_s/hs)/rp/(hn+hs);
                    Vs(i,j) = 2*rhat_s/rp/hs/(hs+hn);

                else %if we cross to the south
                    if phip >= 0 % and if we are on fluid 2
                        divs = eta2*hs + eta1*(dr - hs);
                        Vp(i,j) = Vp(i,j) + 1/dr^2 - 2*rhat_s*eta2/rp/(hn+hs)/divs;
                        Vs(i,j) = 2*rhat_s*eta2/rp/(hn+hs)/divs;
                        b(i,j) = b(i,j) - 2*rhat_s*eta2/rp/(hn+hs)/divs * Ts;

                    else %and if we are on fluid 1
                        divs = eta1*hs + eta2*(dr - hs);
                        Vp(i,j) = Vp(i,j) + 1/dr^2 - 2*rhat_s*eta1/rp/(hn+hs)/divs;
                        Vs(i,j) = 2*rhat_s*eta1/rp/(hn+hs)/divs;
                        b(i,j) = b(i,j) + 2*rhat_s*eta1/rp/(hn+hs)/divs * Ts;
                    end
                end
            end
%             if min([he,hw,hs,hn]) < dr*dr
%                 Tp = front(phi_pad(ip-1:ip+1,jp-1:jp+1),rp,'p');
%                 Vp(i,j) = 1;
%                 Vw(i,j) = 0;
%                 Ve(i,j) = 0;
%                 Vn(i,j) = 0;
%                 Vs(i,j) = 0;
%                 b(i,j) = Tp;
%             end
%         end
%         % Rescale matrix system.
%         %   Neighbours.
%         Ve(i,j) = Ve(i,j) / Vp(i,j);
%         Vw(i,j) = Vw(i,j) / Vp(i,j);
%         Vn(i,j) = Vn(i,j) / Vp(i,j);
%         Vs(i,j) = Vs(i,j) / Vp(i,j);
%         %   RHS.
%         b(i,j) = b(i,j) / Vp(i,j);
%         %   Current node.
%         Vp(i,j) = 1;
    end
end

%Neumann B.C.
Vp(:,end) = Vp(:,end) + (2/dr^2) - 2*(1 - dr/2/R(end,end))/dr^2; 
Vn(:,end) = 0;
Vs(:,end) = 2*(1 - dr/2/R(end,end))/dr^2;
b(:,end) =  Q*eta2 / (dr*pi*R(end,end));

Vp(:,1) = -(1/dr^2); 
Ve(:,1) = 0;
Vw(:,1) = 0;
Vn(:,1) = 1/dr^2;
Vs(:,1) = 0;




% Reference node
Vp(end/2,5) = 1; 
Vs(end/2,5) = 0;
Vn(end/2,5) = 0;
Ve(end/2,5) = 0;
Vw(end/2,5) = 0;
 b(end/2,5) = 0;

% Form full index matrices.
I_mat = [I I I I I];
J_mat = [Jp Je Jw Jn Js];
V_mat = [Vp Ve Vw Vn Vs];
I_mat2 = vertcat(I_mat(:),I_farf(:));
J_mat2 = vertcat(J_mat(:),J_farf(:));
V_mat2 = vertcat(V_mat(:),V_farf(:));

% For sparse matrix and RHS vector.
A = sparse(I_mat2,J_mat2,V_mat2);
% A = sparse(I_mat(:),J_mat(:),V_mat(:));
b = b(:);

% Solve system.
x = A \ b;
T_pol1 = reshape(x,m,n);
T_pol2 = reshape(x,m,n);
T_pol1(phi_pol>0) = NaN;
T_pol2(phi_pol<0) = NaN;


% Convert back to cart
T1 = griddata(RR, TT, T_pol1, X, Y,'cubic');
% T1 = T;
% T2 = T;
% T1 = zeros(Nodes,Nodes);
T2 = griddata(RR, TT, T_pol2, X, Y,'cubic');
% T1(phi>0) = NaN;
% T2(phi<0) = NaN;

    function Tf = front(phi_local,r_p,d)
        % Compute the temperature on the interface.
        
        switch d
            case 'p'
                % Square box for value on centre node.
                % Derivatives.
                phi_t = (phi_local(3,2) - phi_local(1,2)) / (2*dth);
                phi_tt = (phi_local(1,2) - 2*phi_local(2,2) + phi_local(3,2)) / (dth^2);
                phi_r = (phi_local(2,3) - phi_local(2,1)) / (2*dr);
                phi_rr = (phi_local(2,1) - 2*phi_local(2,2) + phi_local(2,3)) / (dr^2);
                phi_rt = (phi_local(3,3) - phi_local(1,3) - phi_local(3,1) + phi_local(1,1)) / (4*dr*dth);
                % Square firstderivatives.
                phi_t2 = phi_t^2;
                phi_r2 = phi_r^2;
                % Grad.
                grad = sqrt(phi_r2 + phi_t2/(r_p^2));
                % Curvature.
                kappa = (phi_r^3/r_p + phi_tt*phi_r^2/r_p^2 + phi_rr*phi_t^2/r_p^2 ...
                    - 2*phi_r*phi_t*phi_rt/r_p^2 + 2*phi_r*phi_t^2/r_p^3) / (grad^3);
                % Condition on front.
                Tf = - sigma * kappa;
            case 'x'
                % Interpolate from two nodes in x direction.
                % Derivatives.
                phi_t = (phi_local(3:4,2) - phi_local(1:2,2)) / (2*dth);
                phi_tt = (phi_local(1:2,2) - 2*phi_local(2:3,2) + phi_local(3:4,2)) / (dth^2);
                phi_r = (phi_local(2:3,3) - phi_local(2:3,1)) / (2*dr);
                phi_rr = (phi_local(2:3,1) - 2*phi_local(2:3,2) + phi_local(2:3,3)) / (dr^2);
                phi_rt = (phi_local(3:4,3) - phi_local(1:2,3) - phi_local(3:4,1) + phi_local(1:2,1)) / (4*dr*dth);
                % Square first derivatives.
                phi_t2 = phi_t.^2;
                phi_r2 = phi_r.^2;
                % Grad.
                grad = sqrt(phi_r2 + phi_t2/(r_p^2));
                % Curvature.
                kappa = ((phi_r.^3)/r_p + phi_tt.*(phi_r.^2)/r_p^2 + phi_rr.*(phi_t.^2)/r_p^2 ...
                    - 2*phi_r.*phi_t.*phi_rt/r_p^2 + 2*phi_r.*(phi_t.^2)/r_p^3) ./ (grad.^3);
                % New spatial increment.
                h = phi_local(2,2) / (phi_local(2,2) - phi_local(3,2)) * dth;
                % Curvature on front.
                kappaf = (kappa(2,1)*h + kappa(1,1)*(dth-h)) / dth;
                % Condition on front.
                Tf = - sigma * kappaf;
            case 'y'
                % Interpolate from two nodes in y direction.
                % Derivatives.
                phi_t = (phi_local(3,2:3) - phi_local(1,2:3)) / (2*dth);
                phi_tt = (phi_local(1,2:3) - 2*phi_local(2,2:3) + phi_local(3,2:3)) / (dth^2);
                phi_r = (phi_local(2,3:4) - phi_local(2,1:2)) / (2*dr);
                phi_rr = (phi_local(2,1:2) - 2*phi_local(2,2:3) + phi_local(2,3:4)) / (dr^2);
                phi_rt = (phi_local(3,3:4) - phi_local(1,3:4) - phi_local(3,1:2) + phi_local(1,1:2)) / (4*dr*dth);
                % Square first derivatives.
                phi_t2 = phi_t.^2;
                phi_r2 = phi_r.^2;
                % Grad.
                grad = sqrt(phi_r2 + phi_t2./(r_p.^2));
                % Curvature                
                kappa = ((phi_r.^3)./r_p + phi_tt.*(phi_r.^2)./r_p.^2 + phi_rr.*(phi_t.^2)./r_p.^2 ...
                    - 2*phi_r.*phi_t.*phi_rt./r_p.^2 + 2*phi_r.*(phi_t.^2)./r_p.^3) ./ (grad.^3);
                % New spatial increment.
                h = phi_local(2,2) / (phi_local(2,2) - phi_local(2,3)) * dr;
                % Curvature on front.
                kappaf = (kappa(1,2)*h + kappa(1,1)*(dr-h)) / dr;
                % Condition on front.
                Tf = - sigma * kappaf;
        end
        
    end

end










































