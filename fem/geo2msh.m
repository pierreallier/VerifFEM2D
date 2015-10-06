function msh_filename = geo2msh(geo_filename,varargin)
% Cette fonction permet de generer un maillage grace a Gmsh depuis Matlab
% et retourne le chemin vers le fichier maillage genere.
%
% Parametre obligatoires:
%   - geo_filename : le chemin vers un fichier .geo e mailler
%
% Parametres optionnels:
%   - msh_filename : le nom du fichier maillage a ecrire (.msh)
%   - pos_filename : un chemin vers un fichier d'adaptation (.pos)
%   - order_val : un entier representant l'ordre des elements [1-4]
%   - scale_val : un flottant representant un facteur d'echelle de la taille des elements
%
%   /!\ L'ordre des parametres optionnel est quelconque sauf si scal_val
%   est entier. Dans ce cas, il faut preciser l'ordre avant le facteur
%   d'echelle.
% Exemple:
%   geo2msh('../meshes/plate_hole.geo')
%   geo2msh('../meshes/plate_hole.geo','test.msh')
%   geo2msh('../meshes/plate_hole.geo',1.5) facteur d'echelle e 1.5
%   geo2msh('../meshes/plate_hole.geo',2) ordre fixe e 2
%   geo2msh('../meshes/plate_hole.geo',1,2) ordre fixe e 1 et facteur d'echelle e 2
%
    % valeurs par defaut
    scale_val = 1;
    order_val = 1;
    pos_filename = '';
    msh_filename = [geo_filename(1:find(geo_filename == '.',1,'last')) 'msh'];
    
    % lecture des parametres optionnel
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
    
    % Verifications des parametres
    assert(exist(geo_filename, 'file') == 2,['Fichier ' geo_filename ' introuvable.']);
    assert(order_val > 0 & order_val <= 4,'L''ordre des elements est incorrect (entier de 1 a 4).');
    
    % Generation de la commande
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