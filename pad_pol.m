function A_pad = pad_pol(A,k,lbrt)
% Pad matrix.
% Can be padded: symmetrically (s);
%                periodically (p);
%                by extending (e);
%                with zeros (z);
%                or no flux (n).

% Default to symmetric padding.
if nargin < 3
    lbrt = 'ssss';
end

% Ensure lbrt has length 4.
if length(lbrt) == 1
    lbrt = [lbrt lbrt lbrt lbrt];
end

% No flux changes to (e) for phi or (s) otherwise.
if strcmp(inputname(1),'phi') | strcmp(inputname(1),'phi_pol') | strcmp(inputname(1),'R')
    lbrt(lbrt=='n') = 'e';
else
    lbrt(lbrt=='n') = 's';
end

% Size variables.
[m,n] = size(A);

% Pad matrix.
switch lbrt
    % All sides the same.
    case 'xx'
        A_pad = A([k+1:-1:2,1:m,m-1:-1:m-k],[k+1:-1:2,1:n,n-1:-1:n-k]);
    case 'xxx'
        A_pad = A([m-k:m-1,1:m,2:k+1],[n-k+1:n,1:n,2:k+1]);
    case 'xxxx'
        A_pad = zeros(size(A)+2*k);
        A_pad(k+1:k+m,k+1:k+n) = A;
        A_pad(:,1:k) = 3*A_pad(:,k+1)*ones(1,k) ...
            - 3*A_pad(:,k+(k+1:-1:2)) + A_pad(:,k+(2*k+1:-2:3));
        A_pad(k+m+(1:k),:) = 3*ones(k,1)*A_pad(k+m,:) ...
            - 3*A_pad(k+m-(1:k),:) + A_pad(k+m-2*(1:k),:);
        A_pad(:,k+n+(1:k)) = 3*A_pad(:,k+n)*ones(1,k) ...
            - 3*A_pad(:,k+n-(1:k)) + A_pad(:,k+n-2*(1:k));
        A_pad(1:k,:) = 3*ones(k,1)*A_pad(k+1,:) ...
            - 3*A_pad(k+(k+1:-1:2),:) + A_pad(k+(2*k+1:-2:3),:);
    otherwise
        % Sides different.
        A_pad = zeros(size(A)+2*k);
        A_pad(k+1:k+m,k+1:k+n) = A;
        switch lbrt(1)
            case 's'
                A_pad(1:k,:) = A_pad(k+(k+1:-1:2),:);
            case 'p'
                A_pad(1:k,:) = A_pad(k+m+1-(k:-1:1),:);
            case 'e'
                A_pad(1:k,:) = 3*ones(k,1)*A_pad(k+1,:) ...
                    - 3*A_pad(k+(k+1:-1:2),:) + A_pad(k+(2*k+1:-2:3),:);
        end
        switch lbrt(2)
            case 's'
                A_pad(:,1:k) = A_pad(:,k+(k+1:-1:2));
            case 'p'
                A_pad(:,1:k) = A_pad(:,k+n+1-(k:-1:1));
            case 'e'
                A_pad(:,1:k) = 3*A_pad(:,k+1)*ones(1,k) ...
                    - 3*A_pad(:,k+(k+1:-1:2)) + A_pad(:,k+(2*k+1:-2:3));
        end
        switch lbrt(3)
            case 's'
                A_pad(k+m+(1:k),:) = A_pad(k+m-(1:k),:);
            case 'p'
                A_pad(k+m+(1:k),:) = A_pad(k+(1:k),:);
            case 'e'
                A_pad(k+m+(1:k),:) = 3*ones(k,1)*A_pad(k+m,:) ...
                    - 3*A_pad(k+m-(1:k),:) + A_pad(k+m-2*(1:k),:);
        end
        switch lbrt(4)
            case 's'
                A_pad(:,k+n+(1:k)) = A_pad(:,k+n-(1:k));
            case 'p'
                A_pad(:,k+n+(1:k)) = A_pad(:,k+(1:k));
            case 'e'
                A_pad(:,k+n+(1:k)) = 3*A_pad(:,k+n)*ones(1,k) ...
                    - 3*A_pad(:,k+n-(1:k)) + A_pad(:,k+n-2*(1:k));
        end
end





























