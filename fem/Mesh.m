classdef Mesh
% Classe r�pr�sentant un objet de type maillage.
%
% Variables :
%   - nodes : les coordonn�es des noeuds (2D)
%   - elems : la table de connectivit� des �l�ments 
%   - type : le type d'�l�ments (0 = Nodes ; 1 = Tress ; 2 = Triangle)

    properties (SetAccess = protected)
        type; % le type d'�l�ments (0 = Nodes ; 1 = Tress ; 2 = Triangle)
        nodes; % les coordonn�es des noeuds (2D)
        elems; % la table de connectivit� des �l�ments 
    end

    properties(Constant,Hidden)
        % correspondance entre le n� de type d'�l�ment de GMSH avec celui
        % de notre code (valeur n�gative pour les �l�ments non impl�ment�)
        GMSH2ELEM_TYPE = [1 2 -1 -1 -1 -1 -1 1 2 -1 -1 -1 -1 -1 0 -1 -1 -1 -1 -1 2 -1 2 -1 -1 1 1];
    end
    
    methods
        %% CONSTRUCTEUR %%
        function obj = Mesh(varargin)
        % Constructeur d'un objet de type maillage.  
        %   Mesh(filename) lit un fichier maillage GMSH � partir de son
        %   chemin filename.
        %   Mesh(coor,connec,type) cr�e un objet malliage � partir des
        %   coordonn�es des noeuds (coor) et la table de connectivit� des
        %   �l�ments (connec) de type (type). Par d�faut type=2.
        
            assert(nargin <= 3,'Bad argument number in Mesh constructor');
            if nargin == 1 % 1er usage de la fonction
                assert(ischar(varargin{1}) && exist(varargin{1},'file') == 2,'Le fichier sp�cifi� n''existe pas');
                obj = obj.read(varargin{1});
            else % 2nd usage
                if nargin == 2
                    obj.type = 2;
                else
                    obj.type = varargin{3}(1);
                    assert(abs(floor(obj.type)-obj.type) <= eps && obj.type <= 2 && obj.type >= 0,'Le type des �l�ments est incorrect');
                end
                
                % V�rification des donn�es
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
                
                assert(obj.order <= 4,'Seul des �l�ments d''ordre compris entre 1 et 4 sont autoris�');
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
        % Le nombre d'�l�ments dans ce maillage
            n = size(obj.elems,1);
        end
        
        function n = order(obj)
        % L'ordre des �l�ments du maillage
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
        % le nombre de noeuds par �l�ments
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
        % Cette fonction calcule la taille caract�ristique de chaque
        % �l�ment pr�sent dans ce maillage et retourne un vecteur de taille
        % le nombre d'�l�ments.
            h = zeros(obj.nbElems,1);
            if obj.type > 0
                for i=1:numel(h)
                    x = obj.nodes(obj.elems(i,1:obj.type+1),:);
                    h(i) = max(sqrt(sum(diff([x;x(1,:)],1).^2,2)));
                end
            end
        end
        
        function h = elementArea(obj)
        % Cette fonction calcule l'aire de chaque �l�ment et retourne un
        % vecteur de taille le nombre d'�l�ments.
            if obj.type ~= 2
                h = obj.elementSize;
            else
                h1 = obj.nodes(obj.elems(:,2),:) - obj.nodes(obj.elems(:,1),:);
                h2 = obj.nodes(obj.elems(:,3),:) - obj.nodes(obj.elems(:,1),:);
                h = 1/2*(h1(:,1).*h2(:,2) - h1(:,2).*h2(:,1));
            end
        end
        
        function [conn,eta,normal] = interfaceElems(obj) % NOT TESTED
        % Calcule la table de connectivit� (conn) des interfaces.
        %
        % Calcule �galement si besoin un tableau (eta) permettant de
        % retrouver la correspondance entre les interface et les �l�ments
        % et la valeur des normales de chaques interface d�crites dans la
        % table de connectivit� conn. 
        %
        % eta est de taille nb interface x nb elems contenant tel que:
        %   - eta(interface_i,elem_j) = 0 si l'interface i n'appartient pas
        %   � l'�l�ment j, +/- 1 sinon.
        %   - +1 si l'interface partag� par cet �l�ment est dans le m�me
        %   sens que celui d�crit dans la table de connectivit� conn � la
        %   ligne i; - -1 sinon (utile pour d�terminer le sens de la
        %   normale)
        %
        % normal est une matrice de taille nb interface x 2 o� chaque ligne
        % k contient la normale [nx ny] de l'interface d�crite � la ligne k
        % de la table conn. eta(:,j)*normal donne la normale de l'ensemble
        % des interfaces apartenant � l'�l�ment j.
        
            if obj.type ~= 2
                error('La fonction interfaceELems ne marche que pour des �l�ments de type triangle');
            end
            
            d = obj.order-1;
            conn = [obj.elems(:,[1,2,3+(1:d)]) ; obj.elems(:,[2,3,3+d+(1:d)]) ; ...
                    obj.elems(:,[3,1,3+2*d+(1:d)])]; 

            [~,border_ids,elems_ids] = unique(sort(conn,2),'rows','stable');
            conn = conn(border_ids,:);

            if nargout >= 2
                % Calcule la table eta
                elems_ids = reshape(elems_ids,[],3);
                repeated_index = reshape(1:numel(elems_ids),[],3) == elems_ids;
                eta = sparse(elems_ids,repmat((1:obj.nbElems)',1,3),repeated_index - ~repeated_index);
            end
            
            if nargout >= 3
                % Calcule les normales
                normal = [-(obj.nodes(conn(:,1),2)-obj.nodes(conn(:,2),2)) (obj.nodes(conn(:,1),1)-obj.nodes(conn(:,2),1))]; 
                n = sqrt(normal(:,1).^2+normal(:,2).^2);
                normal = [normal(:,1)./n normal(:,2)./n];
            end
        end
       
        %% MANIPULATIONS DE MAILLAGE %%
        
        function submesh = restrict(obj,fct)
        % Restreint un maillage pour v�rifier une fonction
        %   restrict(mesh,fct_handle) restreint le maillage � partir d'une
        %   fonction prenant en entr�e les coordonn�es des noeuds et
        %   retournant vrai si le noeuds appartient � ce domaine, faux
        %   sinon. Seul les �l�ments dont tous les noeuds appartiennent au
        %   domaine sont conserv�s.
        
            in_nodes = fct(obj.nodes);
            in_elems = all(in_nodes(obj.elems),2);
            submesh = Mesh(obj.nodes,obj.elems(in_elems,:),obj.type);
        end
        
        function submesh = border(obj,fct)
        % border(mesh) Retourne un sous-maillage d�crivant le bord du domaine
        %
        % border(mesh,fct) Retourne un sous-maillage d�crivant le bord du
        % domaine dont les noeuds v�rifient la fonction fct donn�e (Voir
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
                error('Les �l�ments de type Node n''ont pas de bord');
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
        % Affiche un champs �l�ment par �l�ment moyenn� sur les �l�ments
        % Affiche en 1er la composante \sigma_11, mais un menu acessible
        % par un clic droit permet de choisir les autres composantes
        % \sigma_22 et \sigma_12
            % Plot mesh
            h = plot(obj,varargin{:});
            set(h,'FaceAlpha',1);
            
            dim = 2;
            if numel(data) == obj.nbElems
                dim = 1;
            end
            d = dim*(dim+1)/2;
            
            N = numel(data)/d/obj.nbElems;
            % re-organise les donn�es
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
            set(h,'UserData',data2);
            set(h,'FaceVertexCData',data2(:,1));
            set(h,'FaceColor','flat');
            set(h,'LineStyle','none');
            set(h,'UIContextMenu',c);
        end
        
        function h = plotNodeField(obj,data,varargin)
        % Affiche un champs connu aux noeuds.
        % Affiche en 1er la composante u_1, mais un menu acessible
        % par un clic droit permet de choisir l'autre composantes
        % u_2 (si il y a deux composantes par noeuds).
            % Plot mesh
            h = plot(obj,varargin{:});
            set(h,'FaceAlpha',1);
            
            c = uicontextmenu;
            N = numel(data)/obj.nbNodes;
            if N > 1
                for i=1:N
                    uimenu(c,'Label',['u' num2str(i)],'Callback',{@setdata,h,i});
                end
            end
            data = reshape(data,N,[])';
            set(h,'UserData',data);
            set(h,'FaceVertexCData',data(:,1));
            set(h,'FaceColor','interp');
            set(h,'LineStyle','none');
            set(h,'UIContextMenu',c);
        end
    end
    
    methods(Access = protected)
        function obj = read(obj,filename)
        % Lit un fichier maillage GMSH � partir de son chemin filename et
        % construit un objet Mesh
            fID = fopen(filename,'r');
            while ~feof(fID)
                tline = fgetl(fID);
                if strcmpi(tline,'$Nodes') % d�finition des noeuds
                    nb = sscanf(fgetl(fID), '%d',1);
                    obj.nodes = fscanf(fID,'%*d %f %f %*f',[2 nb])';
                elseif strcmpi(tline,'$Elements') % d�finition des �l�ments
                    nb = sscanf(fgetl(fID), '%d',1);
                    
                    els = zeros(nb,1);
                    obj.type = 0;
                    for i=1:nb
                        data = sscanf(fgetl(fID), '%d');
                        if obj.GMSH2ELEM_TYPE(data(2)) == obj.type
                            els(i,:) = data(data(3)+4:end);
                        elseif obj.GMSH2ELEM_TYPE(data(2)) > obj.type
                            % On ne garde que les �l�ments de plus grande dimension
                            obj.type = obj.GMSH2ELEM_TYPE(data(2));
                            els = zeros(nb,numel(data)-data(3)-3);
                            els(i,:) = data(data(3)+4:end);
                        end
                    end
                    
                    obj.elems = els(any(els,2),:);% suppression des �l�ments inexistant
                end
            end
        end
    end
end

function setdata(source,callbackdata,h,i)
% Fonction n�cessaire pour le menu d�roulant pour la visualisation (voir
% plotElemField et plotNodeField
    data = get(h,'UserData');
    if size(data,2) > 1
        set(h,'FaceVertexCData',data(:,i));
    end
end