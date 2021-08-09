function [Z,VN]=CRHC2_57(A,B,N,M,Q,R,Pf,F1,G1,h1,F2,G2,h2,x0)
% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
% Q, R, and Pf are the gains in the cost function
% N is the length of the horizon
% Z is the vector of optimal variables and VN is the cost function 
% F1, G1, h1, F2, G2, h2 are constraint matrices
% x0 is the initial condition


% Batch parameters

    gam = kron(eye(N),B);
    om = A;
    for i=1:N-1
        gam = gam + kron(diag(ones(N-i,1),-i),A^i*B);
        om = [om; A^(i+1)];
    end
    Qb = blkdiag( kron(eye(N-1),Q), Pf );
    Rb = kron(eye(N),R);
    H = 2*(gam'*Qb*gam + Rb);
    f = (2*x0'*om'*Qb*gam)';
%COnstraints:    
    if ~isempty(F2)
        Ai = F2*gam+G2;
        Bi = h2-F2*om*x0;
    else
        Ai = []; Bi= [];
    end
       
    % Post control horizon, the inputs are kept constant
    Aeq = kron([zeros(N-M,M-1) -1*ones(N-M,1) eye(N-M)], eye(size(B,2)));
    beq = kron( zeros(N-M,1), zeros(size(B,2),1));
    if ~isempty(F1)
        Aeq = [ Aeq ; F1*gam+G1];
        beq = [ beq ; h1-F1*om*x0];
    end

    
    options = optimoptions('quadprog','Display','none');
    [Z,VN,~] = quadprog(H,f, Ai,Bi, Aeq,beq, [],[], [],options);
end