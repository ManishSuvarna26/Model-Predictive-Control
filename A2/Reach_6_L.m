function P=Reach_09(A,B,S,U)
% A and B are the system matrices x^+=Ax+Bu
% S is the polytope for set S
% U is the polytope for feasible inputs
% P is the polytope Reach(S)


    P = A*S + B*U;

    % OBS: only works when there exists inv(A)
%     Ain = [S.A/A    -S.A/A*B
%            zeros(size(U.A,1),size(S.A/A,2)) U.A];
%     bin = [S.b; U.b];
% 
%     Pext = Polyhedron('A',Ain, 'b',bin);
% 
%     P = projection(Pext,1:length(A));
    
    


%     a = A*S;
%     
%     b = Polyhedron('V',(A*S.V')')
%     plot(a,'alpha',0.3,b,'alpha',0.3)

end
