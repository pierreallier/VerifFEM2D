function new_mesh = deform(mesh,u,coeff)
% Créer un maillage déformée à partir d'un champs de déplacement 2D
%
% Paramètres:
%   - le maillage initial
%   - un champs de déplacement 2D
%   - un coefficient d'amplification (optionel)
    
    % Vérification des entrées
    assert(isa(mesh,'Mesh'),'Objet maillage inconnu');
    assert(isnumeric(u) && numel(u) == 2*mesh.nbNodes,'Mauvais champs de déplacement 2D');

    if nargin == 2
        coeff = 1;
    else
        assert(isnumeric(coeff),'Mauvais facteur d''amplification');
        coeff = coeff(1);
    end
    
    if any(size(u) == 1)
        u = reshape(u(:),[],mesh.nbNodes)';
    end
    
    % Création du nouveau maillage
    new_mesh = Mesh(mesh.nodes+coeff*u,mesh.elems);
end