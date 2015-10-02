function path = findGmshPath()
% Cette fonction retourne le chemin de l'ex�cutable de gmsh. 

    % Modifier user_gmsh_path pour forcer un chemin diff�rent
    user_gmsh_path = '';
    
    if ispc 
        % Chemin par d�faut de windows (les stations des salles CAO du DGM)
        gmsh_path = 'C:\Program Files\GMSH\gmsh.exe';
    elseif ismac
        % Chemin par d�faut pour un mac
        gmsh_path = '/Applications/Gmsh.app/Contents/MacOS/gmsh';
    else
        % Chemin par d�faut pour unix (les stations du CdC au LMT)
        gmsh_path = 'gmsh';
    end
    
    if isempty(user_gmsh_path)
        path = gmsh_path;
    else
        path = user_gmsh_path;
    end
    
    % Teste le chemin
    if ~strcmpi(path,'gmsh') && exist(path,'file') ~= 2
        warning('Impossible de trouver le programme GMSH. Est-il install� ?');
    end
end