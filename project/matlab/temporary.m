Tn = zeros((L+2)*(N+2),1)+T0;
Tn(1:N+2) = T_low;

for n = 1:n_max-1
%     Tn = matrixToVector(T(:,:,n));  % Temperature vector at time n
    % 4.1 COMPUTE DISCRETIZATION COEFFICIENTS
    A = zeros((N+2)*(L+2));
    b = zeros((N+2)*(L+2),1);
    % 4.1.1 INTERNAL NODES
    for j = 2:L+1
        for i = 2:N+1
            % Node P ID
            id = i+(j-1)*(N+2);
            % First version of the discretization coefficients
            aW = lambda_w(j,i)*Sy(j)/dx(i-1);
            aE = lambda_e(j,i)*Sy(j)/dx(i);
            aS = lambda_s(j,i)*Sx(i)/dy(j-1);
            aN = lambda_n(j,i)*Sx(i)/dy(j);
            % Net heat on node P at time n
            QPn = aW*Tn(id-1)+aE*Tn(id+1)+aS*Tn(id-(N+2))+aN*Tn(id+(N+2))-(aW+aE+aS+aN)*Tn(id);
            % Discretization coefficients
            aW = beta*aW;
            aE = beta*aE;
            aS = beta*aS;
            aN = beta*aN;
            dP = rho(j,i)*V(j,i)*cp(j,i)/t_step;
            aP = dP+aW+aE+aS+aN;
            bP = dP*Tn(id)+(1-beta)*QPn;
            % Insert coefficients on matrix
            A(id,id) = aP;
            A(id,id-1) = aW;
            A(id,id+1) = aE;
            A(id,id-(N+2)) = aS;
            A(id,id+(N+2)) = aN;
            b(id,1) = bP;
        end
    end
    % 4.1.2 LOWER WALL AND UPPER WALL NODES
    for i = 2:N+1
        % LOWER WALL
        A(i,i) = 1;
        b(i,1) = T_low;
        % UPPER WALL
        % Node P ID
        id = i+(L+1)*(N+2);
        % First version of the discretization coefficients
        aS = lambda_s(L+2,i)*Sx(i)/dy(L+1);
        % Net heat on node P at time n
        QPn = Q_flow*Sx(i)+aS*(Tn(id-(N+2))-Tn(id));
        % Discretization coefficients
        aS = beta*aS;
        aP = aS;
        bP = beta*Q_flow*Sx(i)+(1-beta)*QPn;
        % Insert coefficients on matrix
        A(id,id) = aP;
        A(id,id-(N+2)) = aS;
        b(id,1) = bP;
    end
    % 4.1.3 LEFT WALL AND RIGHT WALL NODES
    for j = 2:L+1
        % LEFT WALL
        % Node P ID
        id = 1+(j-1)*(N+2);
        % First version of the discretization coefficients
        aE = lambda_e(j,1)*Sy(j)/dx(1);
        % Net heat on node P at time n
        QPn = alpha_g*(Tg-Tn(id))*Sy(j)+aE*(Tn(id+1)-Tn(id));
        % Discretization coefficients
        aE = beta*aE;
        aP = aE + beta*alpha_g*Sy(j);
        bP = beta*alpha_g*Tg*Sy(j) + (1-beta)*QPn;
        % Insert coefficients on matrix
        A(id,id) = aP;
        A(id,id+1) = aE;
        b(id,1) = bP;
        % RIGHT WALL
        % Node P ID
        id = j*(N+2);
        % Insert coefficients on matrix
        A(id,id) = 1;
        b(id,1) = T_fun(n*t_step);
    end
    % 4.2 ERASE SINGULAR NODES
    singularNodes = [1,N+2,1+(L+1)*(N+2),(L+2)*(N+2)];
    A(singularNodes,:) = [];
    A(:,singularNodes) = [];
    b(singularNodes,:) = [];
    Tn(singularNodes,:) = [];
    % 4.3 SOLVE LINEAR SYSTEM
%     [T_new, it] = gaussSeidel(A, b, Tn, itmax, delta);
    T_new = linsolve(A,b);
    % 4.4 COMPUTE SINGULAR NODES TEMPERATURE
    T_new = [0; T_new(1:N,1); 0; T_new(N+1:N+L*(N+2),1); 0; T_new(N+L*(N+2)+1:end,1); 0];
    % 4.5 SAVE TEMPERATURE
    T(:,:,n+1) = vectorToMatrix(T_new,L+2,N+2); 
    % 4.6 UPDATE TEMPERATURE
    Tn = T_new;
end