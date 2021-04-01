% A = comp_sol_activity(M,Ac,ei,t0,A0,t) computes the analytic solution 
% of the compartmental ODE's system in terms of activities 
% for the FDG model with nc compartments and known input Ac. 
%..........................................................................
% - ODE's system:
% dA/dt = M*A + Ac*e1 A(t0)=A0
% - Analytic solution:
% A(t) = int_t0^t exp(M(t-u)) * Ac(u) * ei du + exp(M(t-t0)) * A0 
%..........................................................................
% - INPUT:
% M is a nc x nc matrix.
% Ac is a function handle accepting a vector as argument and returning a
% vector of the same size.
% ei is a nc x 1 vector of {0,1}.
% t0 is a scalar.
% A0 is a vector of length nc.
% t is the time vector.
% - OUTPUT:
% A is a matrix of size nc x length(t)

% NB: allow directly the computation of the concentration starting from t0 != 0

function A = f_comp_sol_activity(M,Ac,ei,t0,A0,t)

    ngl = 8; [x,w] = gauss_legendre(ngl);

    nc = length(M);
    nt = length(t);
    ind=find(ei==1);
    
    A = zeros(nc,nt);
    Ac = @(tt)(interp1([0;t],[0;Ac],tt,'linear',0));

    f = @(u)( Ac(u) * getcols(expm((t(1)-u)*M),ind) );
    A(:,1) = expm((t(1)-t0)*M) * A0(:) + quadglv(f,t0,t(1),x,w);

    for n=2:nt 
        
        f = @(u)( Ac(u) * getcols(expm((t(n)-u)*M),ind) );
        A(:,n) = expm((t(n)-t(n-1))*M) * A(:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
end