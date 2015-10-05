function A = FEMVec(omega,F)
% Construit le vecteur b elements finis du systeme Au=b, tel que :
%   b = \int_{omega} N^T.F dx
%
%   F est le vecteur d'effort.

    % Verification des donnees
    assert(isa(omega,'Mesh'),'Parametre #1 invalide : objet maillage type invalide');
    assert(isnumeric(F) && numel(F) == 2,'Mauvaise representation de l''effort');
    
    % Determinant de la matrice jacobienne
    if omega.type == 2
        detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);
    elseif omega.type == 1
        detJ = @(J) sqrt(J(:,1).^2+J(:,2).^2);
    else % c'est un element de type Node
        detJ = @(J) 1;
    end
    
    % Integration numerique de Gauss
    [Wg,Xg] = gaussPoints(omega);
    
    % Creation de la matrice
    A = sparse(2*omega.nbNodes,1);
    for i=1:omega.nbElems
        ids = omega.elems(i,:); % ids des noeuds
        map = bsxfun(@(id,j) (id-1)*2+j,ids(:)',(1:2)'); % index de l'inconnue
        Xe = omega.nodes(ids,:); % coordonnees de l'element
        
        [M1,J] = shapesFunctions(omega,Xg,Xe,0); % Evaluation des fonctions de formes
        M2 = repmat(F(:),numel(Wg),1); % Evaluation de l'effort aux points de Gauss
        
        D = kron(diag(Wg.*detJ(J)),eye(2)); % Matrice pour la quadrature de Gauss
        
        A(map(:),1) = A(map(:),1)+(M1'*D*M2);
    end
end