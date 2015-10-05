function [A,M] = FEMMat(omega,B)
% Construit la matrice A elements finis du systeme Au=b, tel que :
%   A = \int_{omega} [\epsilon(N)]^T*B*[\epsilon(N)] dx
%   et optionnellement
%   M = \int_{omega} N^T*eye(2)*N dx
%
%   B est le tenseur de Hook, ecrit sous forme matriciel avec les notations
%   de Voigt.

    % Verification des donnees
    assert(isa(omega,'Mesh'),'Parametre #1 invalide : objet maillage type invalide');
    assert(isnumeric(B) && (all(size(B) == [3 3]) || all(size(B) == [2 2])),'Mauvaise representation du tenseur de Hook');
    
    A = matrixAssembly(omega,1,B);
    if nargout == 2
        M = matrixAssembly(omega,0,eye(2));
    end
end

function A = matrixAssembly(omega,order,B)
    % Determinant de la matrice jacobienne
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);
    
    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(omega);
    
    % Creation de la matrice
    A = sparse(2*omega.nbNodes,2*omega.nbNodes);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        map = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)'); % index de l'inconnue
        Xe = omega.nodes(ids,:); % coordonnees de l'element
        
        [M1,J] = shapesFunctions(omega,Xg,Xe,order); % Evaluation des fonctions de formes
        
        D = kron(diag(Wg.*detJ(J)),B); % Matrice pour la quadrature de Gauss
        
        A(map(:),map(:)) = A(map(:),map(:))+(M1'*D*M1);
    end
end