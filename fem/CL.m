function [index,u0] = CL(varargin)
% Impose les conditions aux limites en déplacement (Dirichlet).
%
% Paramètres: une sucession de couple de 3 données :
%   - un maillage sur lequel on bloque les noeuds
%   - un scalaire imposant la valeur du déplacement en ces noeuds
%   - un vecteur spécifiant les directions bloquée : [0 1] libre en x, bloqué en y
%
% Sortie:
%   - index : un vecteur logique représentant les inconnues bloquée 
%   - u0 : la valeur imposée aux inconnues bloquées
%
% Exemple :
%   [cl_index,u0] = CL(dleft,0,[1 0],dbottom,1.5,[1 1]);
%   u(cl_index) = u0;

    assert(mod(numel(varargin),3) == 0,'Mauvaise utilisation de la fonction');
    assert(isa(varargin{1},'Mesh'),'Objet de type maillage inconnu');
    
    % Initialisation
    dim = 2;
    index = false(2*varargin{1}.nbNodes,1);
    u0 = zeros(2*varargin{1}.nbNodes,1);
    
    for i=1:numel(varargin)/3
        mesh = varargin{3*(i-1)+1};
        val = varargin{3*(i-1)+2}(1);
        dir = logical(varargin{3*i});
        
        assert(isa(mesh,'Mesh'),'Objet de type maillage inconnu');
        assert(numel(dir) == 2,'Les blocages doivent être donnée en 2D');  
        
        ids = unique(mesh.elems(:)); % les noeuds qui appartiennent au domaine
        map = bsxfun(@(id,j) (id-1)*dim+j,ids(:)',(1:dim)'); % ids de l'inconnue u
        map = map(:).*repmat(dir(:),numel(ids),1); % restriction à la direction bloquée
        map = map(map ~= 0);
        
        index(map) = true;
        u0(map) = val;
    end
    
    u0 = u0(index);
end