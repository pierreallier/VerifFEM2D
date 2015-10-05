function msh_filename = geo2msh(geo_filename,varargin)
% Cette fonction permet de g�n�rer un maillage gr�ce � Gmsh depuis Matlab
% et retourne le chemin vers le fichier maillage g�n�r�.
%
% Param�tre obligatoires:
%   - geo_filename : le chemin vers un fichier .geo � mailler
%
% Param�tres optionnels:
%   - msh_filename : le nom du fichier maillage o� �crire (.msh)
%   - pos_filename : un chemin vers un fichier d'adaptation (.pos)
%   - order_val : un entier repr�sentant l'ordre des �l�ments [1-4]
%   - scale_val : un flottant repr�sentant un facteur d'�chelle de la taille des �l�ments
%
%   /!\ L'ordre des param�tres optionnel est quelconque sauf si scal_val
%   est entier. Dans ce cas, il faut pr�ciser l'ordre avant le facteur
%   d'�chelle.
% Exemple:
%   geo2msh('../meshes/plate_hole.geo')
%   geo2msh('../meshes/plate_hole.geo','test.msh')
%   geo2msh('../meshes/plate_hole.geo',1.5) facteur d'�chelle � 1.5
%   geo2msh('../meshes/plate_hole.geo',2) ordre fix� � 2
%   geo2msh('../meshes/plate_hole.geo',1,2) ordre fix� � 1 et facteur d'�chelle � 2
%
    % valeurs par defaut
    scale_val = 1;
    order_val = 1;
    pos_filename = '';
    msh_filename = [geo_filename(1:find(geo_filename == '.',1,'last')) 'msh'];
    
    % lecture des param�tres optionnel
    string_ids = cellfun(@ischar,varargin);
    float_vals = varargin(cellfun(@(x) abs(floor(x)-x) >= eps,varargin(~string_ids)));
    integer_vals = varargin(cellfun(@(x) abs(floor(x)-x) <= eps,varargin(~string_ids)));
    if ~isempty(float_vals)
        scale_val = float_vals{1};
    end
    if ~isempty(integer_vals)
        order_val = integer_vals{1};
        if numel(integer_vals) > 1
            scale_val = integer_vals{2};
        end
    end
    filenames = varargin(string_ids);
    for i=1:numel(filenames)
        [~,~,ext] = fileparts(filenames{i});
        if strcmpi(ext,'.pos')
            pos_filename = filenames{i};
        else
            msh_filename = filenames{i};
        end
    end
    
    % V�rifications des param�tres
    assert(exist(geo_filename, 'file') == 2,['Fichier ' geo_filename ' introuvable.']);
    assert(order_val > 0 & order_val < 4,'L''ordre des �l�ments est incorrect (entier de 1 � 4).');
    
    % G�n�ration de la commande
    command = ['"' findGmshPath() '" "' geo_filename '" -2 -clscale ' num2str(scale_val) ' -order ' num2str(order_val) ' -o "' msh_filename '"'];
    if ~isempty(pos_filename)
        assert(exist(pos_filename, 'file') == 2,['Fichier ' pos_filename ' introuvable.']);
        command = [ command ' -bgm "' pos_filename '"'];
    end
    
    % Execution de la commande
    disp(command);
    if ~(ispc || ismac)
        setenv('LD_LIBRARY_PATH', '');
    end
    [err,cmdout] = system(command);
    if err
        error('geo2msh:BadCommand',['Une erreur s''est produite, abandon.\n' cmdout]);
    end
    
end