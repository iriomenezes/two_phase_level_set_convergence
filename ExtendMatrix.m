function Xext = ExtendMatrix(X)

[Ny, Nx] = size(X);

Xext  = zeros(Ny + 2, Nx + 2);
Xext(2:end-1,2:end-1) = X;

Xext(1,:) = 2*Xext(2,:)-Xext(3,:);
Xext(:,end) = 2*Xext(:,end-1)-Xext(:,end-2);
Xext(:,1) = 2*Xext(:,2)-Xext(:,3);
Xext(end,:) = 2*Xext(end-1,:)-Xext(end-2,:);

end