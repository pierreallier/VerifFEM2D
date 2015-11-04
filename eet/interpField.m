function u_new = interpField(omega,u,omega_new)
% interpolation un champ de deplacement (aux noeuds) sur un maillage
% augmentee en ordre
%
% u_int = interpField(x,u,x_int) interpole le champ de deplacement (aux
% noeuds) definit sur le maillage x sur le nouveau maillage x_int issu
% d'une augmentation en ordre du maillage x et retourne le champ interpole
% u_int.

    assert(isa(omega,'Mesh') && isa(omega_new,'Mesh'),'Mauvais objet de type maillage');
    assert(omega.type == 2 && omega_new.type == 2,'Seul l''interpolation sur des maillages triangulaires est possible');
    assert(numel(u) == 2*omega.nbNodes,'Seul un champ 2d au noeud peut etre interpole');
    assert(omega.nbElems == omega.nbElems,'Le second maillage n''est pas une simple augmentation d''ordre du premier');
   
    % Retrouver l'id des noeuds sommets des triangles
    sommet_nodes = unique(omega.elems(:,1:3));
    node_matchs = zeros(omega.nbNodes,1);
    for i=1:size(sommet_nodes,1)
        id = find(omega_new.nodes(:,1) == omega.nodes(sommet_nodes(i),1) & omega_new.nodes(:,2) == omega.nodes(sommet_nodes(i),2),1);
        if isempty(id)
            error('Les deux maillages n''ont pas le meme support');
        else
            node_matchs(sommet_nodes(i)) = id;
        end
    end
    
    % Trouver la correspondance entre elements
    elem_matchs = zeros(omega.nbElems,1);
    for i=1:size(elem_matchs,1)
        id = find(omega_new.elems(:,1) == node_matchs(omega.elems(i,1)) & omega_new.elems(:,2) == node_matchs(omega.elems(i,2)) & omega_new.elems(:,3) == node_matchs(omega.elems(i,3)),1);
        if isempty(id)
            error('Les deux maillages n''ont pas le meme support');
        else
            elem_matchs(i) = id;
        end
    end
     
    % maps de correspondances
    map = @(x) kron(2*(x-1),[1;1]) + repmat([1;2],size(x));

    if omega.order ~= omega_new.order
        % Coordonnees locales des noeuds du nouveau maillage dans l'element ancien
        switch omega_new.order
            case 1
               a = [1 0;0 1;0 0]; 
            case 2
               a = [1 0;0 1;0 0;1/2 1/2;0 1/2;1/2 0]; 
            case 3
               a = [1 0;0 1;0 0;2/3 1/3;1/3 2/3;0 2/3;0 1/3;1/3 0;2/3 0;1/3 1/3];  
            case 4
               a = [1 0;0 1;0 0;3/4 1/4;1/2 1/2;1/4 3/4;0 3/4;0 1/2;0 1/4;1/4 0;1/2 0;3/4 0;1/4 1/4;1/2 1/4;1/4 1/2];   
            otherwise
                error('Maillage de degree trop important');
        end

        % Interpoler le deplacement aux nouveau noeuds
        M1 = shapesFunctions(omega,a,[],0);
        u_new(map(omega_new.elems(elem_matchs,:)')) = M1*u(map(omega.elems'));
    else
        u_new(map(node_matchs)) = u;
    end
    u_new = u_new(:);
end

function g = gravityElemCenter(omega)
    x = omega.nodes(:,1);
    y = omega.nodes(:,2);
    
    g = [mean(x(omega.elems(:,1:3)),2) mean(y(omega.elems(:,1:3)),2)];
end