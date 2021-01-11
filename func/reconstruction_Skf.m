function [k1x,k2x,k3x,k4x,Cx,Cxdata,relerr,iter] = reconstruction_Skf(Cdata,Ca,t,t0,C0,k1,k2,k3,k4)

% blood volume fraction (tumor)
Vb = 0.15;
% interstizial volume fraction (tumor)
Vi = 0.3;

% To sum the compartment concentrations
alpha = [1-Vb,1-Vb-Vi];

% Initialization of kinetic parameters
switch nargin
    case 5 % real data
        k1x = rand(1);
        k2x = rand(1);
        k3x = rand(1);
        k4x = 0;
    case 9 % simulation
        Kx = num2cell([k1,k2,k3,k4]+(rand(1,4)-4/5).*[k1,k2,k3,k4]);
        [k1x,k2x,k3x,k4x] = deal(Kx{:});
end
x = [k1x;k2x;k3x;k4x];

% Solve DIRECT PROBLEM (with initial parameters)
Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
Cx = concentration(k1x,Mx,Ca,t0,C0,t);
Cxdata = (alpha*Cx + Vb*Ca(t))';

% Relative error
relerr = norm(Cdata-Cxdata)/norm(Cdata);

%% GAUSS - NEWTON Regularized method

% Iteration numbers
nit = 0;
iter = 0;
cont = 0;
nit_max = 30;

% Initializations
vect_relerr = zeros(3*nit_max,1);
cell_K = cell(1,3*nit_max);

% Lower bound for the relative error
toll = sqrt(sum(Cdata))/norm(Cdata)*5e-1;

while relerr>toll
    
    % Count iteration's number
    nit = nit+1; % nit = number of iteration for a single choise of the initial parameters
    iter = iter+1; % iter = numer of iteration for the whole algorithm
    
    % Derivative with respect to the parameters
    aux1 = dF_dk1(Mx,Ca,t0,t,alpha)';
    auxM = dF_dM(k1x,Mx,Ca,t0,C0,t,alpha)';
    D = [aux1 auxM];
    D = [D(:,1),-D(:,2),D(:,3)-D(:,2),D(:,4)-D(:,5)];
    
    % Regularization parameter --------------------------------------------
    % values for lambda to compute V_lambda
    switch nargin
        case 5 % real data
            vect_lambda = 5e4:5e4:5e6;
        case 9 % simulation
            vect_lambda = 5e2:5e2:5e4;
    end
    vect_lambda = vect_lambda';
    % Generalized Cross Validation
    [r,~] = GCV(D,Cdata-Cxdata,vect_lambda);
    %----------------------------------------------------------------------
    
    % '\' solve the linear system Dh=Z with Z=Cdata-Cxdata h=increment
    % ==> Tychonov regularization: (rI+D'*D)h=D'*Z
    h = (r*diag([1,1,1,1])+D.'*D)\(D.'*(Cdata-Cxdata));
    
    % Check for positive coefficients
    xph = x+h;
    if any(xph<=0)
        h(xph<=0)=0;
    end
    x = x+h;
    
    % Refresh the parameters
    k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4);
    
    % Solve DIRECT PROBLEM  --> refresh data
    Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
    Cx = concentration(k1x,Mx,Ca,t0,C0,t);
    Cxdata = (alpha*Cx + Vb*Ca(t))';
    
    % Relative error
    [relerrprec,relerr] = deal(relerr,norm(Cdata-Cxdata)/norm(Cdata));
    
    % Store the relative error 'relerr' and the solution 'x' for each iteration 'iter'
    vect_relerr(iter) = relerr;
    cell_K{iter} = x;
    
    %.....................................................................%
    if  ( nit>=15 && relerr>=0.3 ) || ( relerr>=0.3 && abs(relerr-relerrprec)<1e-3 ) || ( nit>=25 && abs(relerr-relerrprec)<1e-4 ) || (nit>=nit_max)
        
        nit = 0;
        cont = cont+1;
        
        if cont == 3
            
            vect_relerr(vect_relerr==0) = [];
            x = cell_K{vect_relerr==min(vect_relerr)};
            k1x = x(1); k2x = x(2); k3x = x(3); k4x = x(4);
            
            % Solve DIRECT PROBLEM (with the parameters leading to the smallest relerr)
            Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
            Cx = concentration(k1x,Mx,Ca,t0,C0,t);
            Cxdata = (alpha*Cx + Vb*Ca(t))';
            
            relerr = min(vect_relerr);
            
            break
        else
            
            % Initialization of kinetic parameters
            switch nargin
                case 5 % real data
                    k1x = rand(1);
                    k2x = rand(1);
                    k3x = rand(1);
                    k4x = 0;
                case 9 % simulation
                    Kx = num2cell([k1,k2,k3,k4]+(rand(1,4)-4/5).*[k1,k2,k3,k4]);
                    [k1x,k2x,k3x,k4x] = deal(Kx{:});
            end
            x = [k1x;k2x;k3x;k4x];
            
            % Solve DIRECT PROBLEM (with initial parameters)
            Mx = [[-(k2x+k3x);k3x],[k4x;-k4x]];
            Cx = concentration(k1x,Mx,Ca,t0,C0,t);
            Cxdata = (alpha*Cx + Vb*Ca(t))';
            
            % Relative error
            relerr = norm(Cdata-Cxdata)/norm(Cdata);
            
            continue
            
        end
        
    end
  
%.........................................................................%
%     figure(1)
%     plot(t,Cdata,'b')
%     hold on
%     plot(t,Cxdata,'c--')
%     
%     pause
%.........................................................................%
    
end

end
