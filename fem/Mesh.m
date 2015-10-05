classdef Mesh
% Classe réprésentant un objet de type maillage.
%
% Variables :
%   - nodes : les coordonnées des noeuds (2D)
%   - elems : la table de connectivité des éléments 
%   - type : le type d'éléments (0 = Nodes ; 1 = Tress ; 2 = Triangle)

    properties (SetAccess = protected)
        type; % le type d'éléments (0 = Nodes ; 1 = Tress ; 2 = Triangle)
        nodes; % les coordonnées des noeuds (2D)
        elems; % la table de connectivité des éléments 
    end

    properties(Constant,Hidden)
        % correspondance entre le n° de type d'élément de GMSH avec celui
        % de notre code (valeur négative pour les éléments non implémenté)
        GMSH2ELEM_TYPE = [1 2 -1 -1 -1 -1 -1 1 2 -1 -1 -1 -1 -1 0 -1 -1 -1 -1 -1 2 -1 2 -1 -1 1 1];
    end
    
    methods
        %% CONSTRUCTEUR %%
        function obj = Mesh(varargin)
        % Constructeur d'un objet de type maillage.  
        %   Mesh(filename) lit un fichier maillage GMSH à partir de son
        %   chemin filename.
        %   Mesh(coor,connec,type) crée un objet malliage à partir des
        %   coordonnées des noeuds (coor) et la table de connectivité des
        %   éléments (connec) de type (type). Par défaut type=2.
        
            assert(nargin <= 3,'Bad argument number in Mesh constructor');
            if nargin == 1 % 1er usage de la fonction
                assert(ischar(varargin{1}) && exist(varargin{1},'file') == 2,'Le fichier spécifié n''existe pas');
                obj = obj.read(varargin{1});
            else % 2nd usage
                if nargin == 2
                    obj.type = 2;
                else
                    obj.type = varargin{3}(1);
                    assert(abs(floor(obj.type)-obj.type) <= eps && obj.type <= 2 && obj.type >= 0,'Le type des éléments est incorrect');
                end
                
                % Vérification des données
                assert(isnumeric(varargin{1}) && isnumeric(varargin{2}),'Mauvais arguments');
                if size(varargin{1},1) == 0
                    obj.nodes = zeros(0,2);
                    warning('Empty mesh');
                else
                    obj.nodes = varargin{1}(:,1:2);
                end
                
                overflow_ids = any(varargin{2} > obj.nbNodes & varargin{2} < 1,2);
                if any(overflow_ids)
                    warning('Some elements are ignored');
                end
                obj.elems = floor(varargin{2}(~overflow_ids,:));
                
                assert(obj.order <= 4,'Seul des éléments d''ordre compris entre 1 et 4 sont autorisé');
            end
        end
        
        %% VARIABLES SUPPLEMENTAIRES %%
        function n = nbNodes(obj)
        % Le nombre de noeuds dans ce maillage
            n = size(obj.nodes,1);
        end
        
        function n = dim(obj)
        % La dimension de ce maillage (2D)
            n = size(obj.nodes,2);
        end
        
        function n = nbElems(obj)
        % Le nombre d'éléments dans ce maillage
            n = size(obj.elems,1);
        end
        
        function n = order(obj)
        % L'ordre des éléments du maillage
            switch obj.type
                case 0 % Nodes
                    n = size(obj.elems,2);
                case 1 % Tress
                    n = size(obj.elems,2) - 1;
                case 2 % Triangles
                    n = (sqrt(1+8*size(obj.elems,2))-1)/2-1;
            end
        end
        
        function n = nbNodesPerElems(obj)
        % le nombre de noeuds par éléments
            ord = obj.order; 
            switch obj.type
                case 0 % Nodes
                    n = 1;
                case 1 % Tress
                    n = ord + 1;
                case 2 % Triangles
                    n = (ord+1)*(ord+2)/2;
            end
        end
        
        function h = elementSize(obj)
        % Cette fonction calcule la taille caractéristique de chaque
        % élément présent dans ce maillage et retourne un vecteur de taille
        % le nombre d'éléments.
            h = zeros(obj.nbElems,1);
            if obj.type > 0
                for i=1:numel(h)
                    x = obj.nodes(obj.elems(i,1:obj.type+1),:);
                    h(i) = max(sqrt(sum(diff([x;x(1,:)],1).^2,2)));
                end
            end
        end
        
        function h = elementArea(obj)
        % Cette fonction calcule l'aire de chaque élément et retourne un
        % vecteur de taille le nombre d'éléments.
            if obj.type ~= 2
                h = obj.elementSize;
            else
                h1 = obj.nodes(obj.elems(:,2),:) - obj.nodes(obj.elems(:,1),:);
                h2 = obj.nodes(obj.elems(:,3),:) - obj.nodes(obj.elems(:,1),:);
                h = 1/2*(h1(:,1).*h2(:,2) - h1(:,2).*h2(:,1));
            end
        end
        
        function [conn,eta,normal] = interfaceElems(obj) % NOT TESTED
        % Calcule la table de connectivité (conn) des interfaces.
        %
        % Calcule également si besoin un tableau (eta) permettant de
        % retrouver la correspondance entre les interface et les éléments
        % et la valeur des normales de chaques interface décrites dans la
        % table de connectivité conn. 
        %
        % eta est de taille nb interface x nb elems contenant tel que:
        %   - eta(interface_i,elem_j) = 0 si l'interface i n'appartient pas
        %   à l'élément j, +/- 1 sinon.
        %   - +1 si l'interface partagé par cet élément est dans le même
        %   sens que celui décrit dans la table de connectivité conn à la
        %   ligne i; - -1 sinon (utile pour déterminer le sens de la
        %   normale)
        %
        % normal est une matrice de taille nb interface x 2 où chaque ligne
        % k contient la normale [nx ny] de l'interface décrite à la ligne k
        % de la table conn. eta(:,j)*normal donne la normale de l'ensemble
        % des interfaces apartenant à l'élément j.
        
            if obj.type ~= 2
                error('La fonction interfaceELems ne marche que pour des éléments de type triangle');
            end
            
            d = omega.order-1;
            conn = [omega.elems(:,[1,2,3+(1:d)]) ; omega.elems(:,[2,3,3+d+(1:d)]) ; ...
                    omega.elems(:,[3,1,3+2*d+(1:d)])]; 

            [~,border_ids,elems_ids] = unique(sort(conn,2),'rows','stable');
            conn = conn(border_ids,:);

            if nargout >= 2
                % Calcule la table eta
                elems_ids = reshape(elems_ids,[],3);
                repeated_index = reshape(1:numel(elems_ids),[],3) == elems_ids;
                eta = sparse(elems_ids,repmat((1:omega.nbElems)',1,3),repeated_index - ~repeated_index);
            end
            
            if nargout >= 3
                % Calcule les normales
                normal = [-diff(obj.nodes(conn(:,1:2),2)) diff(obj.nodes(conn(:,1:2),1))]; 
            end
        end
        
        % TODO Patchs
        function n = nbPatchs(obj)
            
        end
        
        function ids = patch(obj,node_id) % NOT TESTED
        % Retourne l'id des éléments appartenant aux patch liée au noeuds node_id
            tmp = 1:obj.nbElems;
            ids = tmp(any(obj.elems == node_id,2));
        end
       
        %% MANIPULATIONS DE MAILLAGE %%
        
        function submesh = restrict(obj,fct)
        % Restreint un maillage pour vérifier une fonction
        %   restrict(mesh,fct_handle) restreint le maillage à partir d'une
        %   fonction prenant en entrée les coordonnées des noeuds et
        %   retournant vrai si le noeuds appartient à ce domaine, faux
        %   sinon. Seul les éléments dont tous les noeuds appartiennent au
        %   domaine sont conservés.
        
            in_nodes = fct(obj.nodes);
            in_elems = all(in_nodes(obj.elems),2);
            submesh = Mesh(obj.nodes,obj.elems(in_elems,:),obj.type);
        end
        
        function submesh = border(obj,fct)
        % border(mesh) Retourne un sous-maillage décrivant le bord du domaine
        %
        % border(mesh,fct) Retourne un sous-maillage décrivant le bord du
        % domaine dont les noeuds vérifient la fonction fct donnée (Voir
        % fonction restriction)
        %
        % Exemple :
        %   omega.border(@(x) x(:,1) == 0) Le bords tel que x == 0

            ord = obj.order;
            if obj.type == 2
                els = cat(1,obj.elems(:,[1:2 3+(1:ord-1)]), ...
                            obj.elems(:,[2:3 2+ord+(1:ord-1)]), ...
                            obj.elems(:,[3 1 1+2*ord+(1:ord-1)])); % construct all 1D elements
                [els,~,ids] = unique(sort(els,2),'rows','stable');
                submesh = Mesh(obj.nodes,els(accumarray(ids,1) == 1,:),1);
            elseif obj.type == 1
                els = obj.elems(:,1:2);
                border_nodes = accumarray(els(:),1) == 1;
                submesh = Mesh(obj.nodes,find(border_nodes),0);
            else
                error('Les éléments de type Node n''ont pas de bord');
            end
            if nargin > 1
                submesh = submesh.restrict(fct);
            end
        end
        
        %% VISUALISATION %%
        function h = plot(obj,varargin)
        % Affiche le maillage.
            if obj.type > 0
                h = patch('Faces',obj.elems(:,1:1+obj.type),'Vertices',obj.nodes,'FaceAlpha',0,varargin{:});
            else
                h = scatter(obj.nodes(obj.elems(:),1),obj.nodes(obj.elems(:),2),'filled','MarkerFaceColor','k',varargin{:});
            end
        end
        
        function h = plotElemField(obj,data,varargin)
        % Affiche un champs élément par élément moyenné sur les éléments
        % Affiche en 1er la composante \sigma_11, mais un menu acessible
        % par un clic droit permet de choisir les autres composantes
        % \sigma_22 et \sigma_12
            % Plot mesh
            h = plot(obj,varargin{:});
            h.FaceAlpha = 1;
            
            dim = 2;
            if numel(data) == obj.nbElems
                dim = 1;
            end
            d = dim*(dim+1)/2;
            
            N = numel(data)/d/obj.nbElems;
            % re-organise les données
            data2 = zeros(obj.nbElems,d);
            for i=1:d
                data2(:,i) = sum(reshape(data(i:d:end),N,[]),1)/N;
            end
            % mise en place du menu
            c = uicontextmenu;
            if dim > 1
                uimenu(c,'Label','s11','Callback',{@setdata,h,1});
            end
            if dim >= 2
                uimenu(c,'Label','s22','Callback',{@setdata,h,2});
                uimenu(c,'Label','s12','Callback',{@setdata,h,3});
            end
            h.UserData = data2;
            h.FaceVertexCData = data2(:,1);
            h.FaceColor = 'flat';
            h.LineStyle = 'none';
            h.UIContextMenu = c;
        end
        
        function h = plotNodeField(obj,data,varargin)
        % Affiche un champs connu aux noeuds.
        % Affiche en 1er la composante u_1, mais un menu acessible
        % par un clic droit permet de choisir l'autre composantes
        % u_2 (si il y a deux composantes par noeuds).
            % Plot mesh
            h = plot(obj,varargin{:});
            h.FaceAlpha = 1;
            
            c = uicontextmenu;
            N = numel(data)/obj.nbNodes;
            if N > 1
                for i=1:N
                    uimenu(c,'Label',['u' num2str(i)],'Callback',{@setdata,h,i});
                end
            end
            h.UserData = reshape(data,N,[])';
            h.FaceVertexCData = h.UserData(:,1);
            h.FaceColor = 'interp';
            h.LineStyle = 'none';
            h.UIContextMenu = c;
        end
    end
    
    methods(Access = protected)
        function obj = read(obj,filename)
        % Lit un fichier maillage GMSH à partir de son chemin filename et
        % construit un objet Mesh
            fID = fopen(filename,'r');
            while ~feof(fID)
                tline = fgetl(fID);
                if strcmpi(tline,'$Nodes') % définition des noeuds
                    nb = sscanf(fgetl(fID), '%d',1);
                    obj.nodes = fscanf(fID,'%*d %f %f %*f',[2 nb])';
                elseif strcmpi(tline,'$Elements') % définition des éléments
                    nb = sscanf(fgetl(fID), '%d',1);
                    
                    els = zeros(nb,1);
                    obj.type = 0;
                    for i=1:nb
                        data = sscanf(fgetl(fID), '%d');
                        if obj.GMSH2ELEM_TYPE(data(2)) == obj.type
                            els(i,:) = data(data(3)+4:end);
                        elseif obj.GMSH2ELEM_TYPE(data(2)) > obj.type
                            % On ne garde que les éléments de plus grande dimension
                            obj.type = obj.GMSH2ELEM_TYPE(data(2));
                            els = zeros(nb,numel(data)-data(3)-3);
                            els(i,:) = data(data(3)+4:end);
                        end
                    end
                    
                    obj.elems = els(any(els,2),:);% suppression des éléments inexistant
                end
            end
        end
    end
end

function setdata(source,callbackdata,h,i)
% Fonction nécessaire pour le menu déroulant pour la visualisation (voir
% plotElemField et plotNodeField
    if size(h.UserData,2) > 1
        h.FaceVertexCData = h.UserData(:,i);
    end
end