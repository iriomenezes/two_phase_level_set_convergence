function CC = driver_HS_twophase
if nargin == 0
    clear, clc
end
tsave = 0;
tinteral = 5;

N = 100; domain = 7.5;
sigma = 1/2000;

eta2 = 1;
eta1 = 1/10;
k_start = 0; time = 0;
tmax = 81;
lbrt = 'n';
yspace = linspace(-domain,domain,N);
dy = yspace(2) - yspace(1);
xspace = linspace(-domain,domain,N)';
dx = dy;
[X, Y] = meshgrid(xspace,yspace);
[TT, RR] = cart2pol(X, Y);
amp = 0.1;
phi = (-1 + RR + amp.*(cos(6.*TT)))';
phi = reinit_crossing_nan(phi,dx,dy,lbrt);

CC = [];
plots = 1;
linestyle = 'b-';
filename = [];


% Initial plot.
if plots ~= 0
    if plots == 1
        figure, h1 = gca;
    else
        h1 = plots;
    end
    box(h1,'on')
end
if isempty(CC)
    plot_sim(1,time)
else
    plot_sim(0,time)
end

plot_sim(0,time)

%% Simulation.

k = k_start;
dt=1e-16;
Q = 1;
while time < tmax
    %First Step
    [T1,T2] = pressure_pol_c(phi,sigma,X,Y,domain,eta1,eta2,Q);
    [dP,F] = velocity_HS(phi,T1,T2,dx,dy,lbrt,eta1,eta2);
    dt = min(0.999*tinteral, dx/(4*max(abs(F(:)))));
    phi1 = phi - dt * dP;    
    [T1,T2] = pressure_pol_c(phi1,sigma,X,Y,domain,eta1,eta2,Q);
    [dP,~] = velocity_HS(phi1,T1,T2,dx,dy,lbrt,eta1,eta2);
    phi2 = phi1 - dt * dP;    
    phi = 0.5*(phi+phi2);    
    time = time + dt; [time, dt]

    %Reinitialization
    if ~mod(k,4)
        for index = 1:5
            [PhiX, PhiY] = gradient(phi, dx, dy);
            PhiGrad = sqrt( PhiX.^2 + PhiY.^2 );
            
            F = phi./sqrt(phi.^2 + (PhiGrad.*dx).^2);
            phiExt = ExtendMatrix(phi);
            phiExt2 = ExtendMatrix(phiExt);
        
            DxP = (-phiExt2(5:end, 3:end-2) + 4.*phiExt2(4:end-1, 3:end-2) - 3.*phiExt2(3:end-2, 3:end-2))./(2.*dx);
            DxM = (3.*phiExt2(3:end-2, 3:end-2) - 4.*phiExt2(2:end-3, 3:end-2) + phiExt2(1:end-4, 3:end-2))./(2.*dx);
            DyP = (-phiExt2(3:end-2, 5:end) + 4.*phiExt2(3:end-2, 4:end-1) - 3.*phiExt2(3:end-2, 3:end-2))./(2.*dy);
            DyM = (3.*phiExt2(3:end-2, 3:end-2) - 4.*phiExt2(3:end-2, 2:end-3) + phiExt2(3:end-2, 1:end-4))./(2.*dy);
        
            Np = sqrt(max(DxM,0).^2+min(DxP,0).^2+max(DyM,0).^2+min(DyP,0).^2);
            Nm = sqrt(min(DxM,0).^2+max(DxP,0).^2+min(DyM,0).^2+max(DyP,0).^2);
        
            DwP   = (-phiExt2(5:end,1:end-4) + 4.*phiExt2(4:end-1,2:end-3) - 3.*phiExt2(3:end-2,3:end-2))/(2.*sqrt(2)*dx);
            DwM   = (3.*phiExt2(3:end-2,3:end-2) - 4.*phiExt2(2:end-3,4:end-1) + 1.*phiExt2(1:end-4,5:end))/(2.*sqrt(2)*dx);
        
            DzP   = (-phiExt2(5:end,5:end) + 4.*phiExt2(4:end-1,4:end-1) -3.*phiExt2(3:end-2,3:end-2))./(2.*sqrt(2)*dx);
            DzM = (3.*phiExt2(3:end-2,3:end-2) -4.*phiExt2(2:end-3,2:end-3) + 1.*phiExt2(1:end-4,1:end-4) )./(2.*sqrt(2)*dx);
        
            Np2 = sqrt(max(DwM,0).^2+min(DwP,0).^2+max(DzM,0).^2+min(DzP,0).^2);
            Nm2 = sqrt(min(DwM,0).^2+max(DwP,0).^2+min(DzM,0).^2+max(DzP,0).^2);
        
            dphi = 0.5.*(max(F,0).*(Np2-1)+min(F,0).*(Nm2-1) + max(F,0).*(Np-1)+min(F,0).*(Nm-1));

            phi = phi - (dx./5).*dphi;
        end
    end
    

    % Plotting.
    if ~mod(k,10)
       plot_sim(0,time)
    end
    if time>tsave
        plot_sim(0,time)
        if ~isempty(filename)
            param.time = time;
            param.N = N;
            param.amp = amp;
            param.sigma = sigma;
            param.domain = domain;
            param.eta1 = eta1;
            param.eta2 = eta2;
            param.X = X;
            param.Y = Y;
            param.phi = phi;
            str = ['MATs/' filename '_' num2str(ceil(time))];
            str(str=='.') = '-';
            save(str,'param')
            str = ['MATs/' filename '_' num2str(ceil(time)) '.mat'];
            str(str=='.') = '-';
            delete(str)
        end
        tsave = tsave + tinteral;
    end
    k = k + 1; 
end


%% Plotting Functions
% Plot function used during simulation.
    function plot_sim(sp,t)
        cc = contourc(X(1,:),Y(:,1),phi,[0 0]);
        if plots ~= 0
            cla(h1), plot_contour_radial([CC cc],h1,linestyle);
            axis equal
            axis([X(1,1) X(end,end) Y(1,1) Y(end,end)]);
            title(string(t))
            drawnow
        end
        if sp
            CC = [CC cc];
        end
    end

% Plot function used after simulation.
    function plot_end
        cla(h1);
        str = [];
        plot_contour_radial(CC,h1,linestyle,str);
        axis equal
        axis([X(1,1) X(end,end) Y(1,1) Y(end,end)]);
        title(string(time));
        drawnow
    end
end




















