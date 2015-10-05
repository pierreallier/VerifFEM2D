function new_mesh = deform(mesh,u,coeff)
% Creer un maillage deformee a partir d'un champs de deplacement 2D
%
% Parametres:
%   - le maillage initial
%   - un champs de deplacement 2D
%   - un coefficient d'amplification (optionel)
    
    % Verification des entrees
    assert(isa(mesh,'Mesh'),'Objet maillage inconnu');
    assert(isnumeric(u) && numel(u) == 2*mesh.nbNodes,'Mauvais champs de deplacement 2D');

    if nargin == 2
        coeff = 1;
    else
        assert(isnumeric(coeff),'Mauvais facteur d''amplification');
        coeff = coeff(1);
    end
    
    if any(size(u) == 1)
        u = reshape(u(:),[],mesh.nbNodes)';
    end
    
    % Creation du nouveau maillage
    new_mesh = Mesh(mesh.nodes+coeff*u,mesh.elems);
end