function [omega_adm,sigma_adm] = eet(omega,sigma,B,dfomega,Fd,clu1,clu2,k_ord)
% Calcul le champ de contrainte admissible par la methode EET
%
% [x_adm,s_adm] = eet(x,s,K,dfx,F,clu,m) calcul le champ de contrainte
% admissible par la methode EET sur le domaine maille x, dont la contrainte
% element finis est s, tel que s = K*epsilon(u) (K l'operateur de Hook dans
% la notation de Voigt). Cette solution doit verifier les conditions aux
% limites en effort sur le bord dfx, de valeur le vecteur 2D F et u.e1=0
% sur clu1 et u.e2=0 sur cl2. La contrainte admissible s_adm est alors
% calculee sur un maillage admissible x_adm augmente en ordre, de valeur
% p+m ou p est l'ordre du maillage x et m l'ordre choisi par l'utilisateur.

    %% Verifications
    assert(isa(omega,'Mesh') && isa(dfomega,'Mesh') && isa(clu1,'Mesh') && isa(clu2,'Mesh'),'Objet mailalge invalide');
    assert(mod(numel(sigma),omega.nbElems*3) == 0,'Mauvais champ de contraintes');
    assert(isnumeric(B) && size(B,1) == size(B,2) && size(B,1) == 3,'Mauvais tenseur de Hook');
    assert(isnumeric(Fd) && numel(Fd) == 2,'Mauvais chargement');
    assert(abs(ceil(k_ord) - k_ord) < eps && k_ord >= 0 && k_ord <= 3,'Mauvais ordre d''augmentation');
    
    %% Initialisation
    % Recuperation des interfaces
    domega = omega.border;
    [interf,eta,n] = omega.interfaceElems;
    
    % Recuperation du bord ou le chargement est impose
    border_nodes_ids = unique(omega.border.elems(:));
    clu_nodes_ids = unique(union(clu1.elems(:),clu2.elems(:)));
    loading_interf_ids = ismember(sort(interf,2),sort(dfomega.elems,2),'rows');
    border_interf_ids = ismember(sort(interf,2),sort(domega.elems,2),'rows');
    
    cl1_interf_ids = ismember(sort(interf,2),sort(clu1.elems,2),'rows');
    cl2_interf_ids = ismember(sort(interf,2),sort(clu2.elems,2),'rows');
    
    clf_interf_ids = not(cl1_interf_ids | cl2_interf_ids) & border_interf_ids;
    
    %% Calcul de la condition de prolongement \int_{E} \omega_h \grad(\phi_i) dx
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);
    [Wg,Xg] = gaussPoints(omega);
    Ng = numel(Wg);

    Q = zeros(size(omega.elems,2)*2,omega.nbElems); % stocke par [Fx^1 Fy^1 Fx^2 Fy^2 Fx^3 Fy^3 ...]' 1,2,3 l'id local des noeuds des elements pour chaque element
    for i=1:omega.nbElems
        nodes_ids = omega.elems(i,:);
        Xe = omega.nodes(nodes_ids,:);

        M1 = sigma((i-1)*Ng*3 + (1:Ng*3));
        [M2,J] = shapesFunctions(omega,Xg,Xe,1);
        D = kron(diag(Wg.*detJ(J)),eye(3));

        Q(:,i) = M2'*D*M1;
    end
    
    %% Calcul des conditions de Neumann \int_{G} \omega_E.n_E \phi_i dx | \int_{G} Fd \phi_i dx
    mean_sigma = reshape([mean(reshape(sigma(1:3:end),Ng,[]),1)' mean(reshape(sigma(2:3:end),Ng,[]),1)' mean(reshape(sigma(3:3:end),Ng,[]),1)']',[],1);
    detJ = @(J) sqrt(J(:,1).^2+J(:,2).^2);
    [Wg,Xg] = gaussPoints(domega);
    Ng = numel(Wg);

    q = zeros(size(interf,2)*2,size(interf,1));
    for i=1:size(interf,1)
        nodes_ids = interf(i,:);
        Xe = omega.nodes(nodes_ids,:);

        [M2,J] = shapesFunctions(domega,Xg,Xe,0);

        if loading_interf_ids(i) % bord sur lequel on impose un effort
            D = kron(diag(Wg.*detJ(J)),eye(2));
            M1 = repmat(Fd(:),Ng,1);
            q(:,i) = M2'*D*M1;
        else
            D = kron(diag(Wg.*detJ(J)),[n(i,1) 0 n(i,2);0 n(i,2) n(i,1)]);
            elems_ids = find(eta(i,:));
            for j=elems_ids
                M1 = repmat(mean_sigma((j-1)*3 + (1:3)),Ng,1);
                q(:,i) = q(:,i) + M2'*D*M1;
            end
            q(:,i) = q(:,i)/numel(elems_ids);
        end
        
        if clf_interf_ids(i) && not(loading_interf_ids(i))
            q(:,i) = 0;
        elseif cl1_interf_ids(i)
            q(2:2:end,i) = 0;
        elseif cl2_interf_ids(i)
            q(1:2:end,i) = 0;
        end
    end
    
    %% Calcul des flux equilibre pour chaque noeud
    interf_nodes = unique(interf(:));
    
    b = zeros(2*size(interf,2),size(interf,1)); % stocke par [bx^1 by^1 bx^2 by^2 ...]' 1,2 l'id local des noeuds des interfaces pour chaque interface
    for k=1:numel(interf_nodes)
        i = interf_nodes(k);
        patch_elems_ids = any(omega.elems == i,2); % elements qui partagent ce noeud
        patch_borders_ids = any(interf == i,2); % bords qui partagent ce noeud

        Bn = kron(eta(patch_borders_ids,patch_elems_ids)',eye(2));
        Qn = Q(logical(kron((omega.elems == i)',[1;1])));
        
        % blocage du noyau
        if all(clu_nodes_ids ~= i)  
            % On retire 2 equations de compatibilite si noeud interieur ou
            % chargement pur
            Bn = Bn(1:end-2,:);
            Qn = Qn(1:end-2);
        elseif any(cl1_interf_ids(patch_borders_ids))
            Bn = Bn(1:end-1,:);
            Qn = Qn(1:end-1);
        elseif any(cl2_interf_ids(patch_borders_ids))
            Bn = Bn([1:end-2 end],:);
            Qn = Qn([1:end-2 end]);
        end

        % Ajout des conditions de chargement
        if any(border_nodes_ids == i)
            Cn = diag(border_interf_ids(patch_borders_ids));
            Cn = kron(Cn(any(Cn,2),:),eye(2));

            qn = q(logical(kron((interf == i)',[1;1])));
            % Mise a 0 des qn qui ne sont sur le bord du chargement
            qn(~logical(kron(loading_interf_ids(patch_borders_ids),[1;1]))) = 0;
            % restriction aux bords
            qn = qn(logical(kron(border_interf_ids(patch_borders_ids),[1;1])));
            
            % Suppression des relations bloquee en deplacement
            if any(cl1_interf_ids(patch_borders_ids) & border_interf_ids(patch_borders_ids))
                id = cl1_interf_ids(patch_borders_ids);
                id = find(id(border_interf_ids(patch_borders_ids)));
                Cn(2*(id-1)+1,:) = NaN;
                qn(2*(id-1)+1) = NaN;
            end
            if any(cl2_interf_ids(patch_borders_ids) & border_interf_ids(patch_borders_ids))
                id = cl2_interf_ids(patch_borders_ids);
                id = find(id(border_interf_ids(patch_borders_ids)));
                Cn(2*id,:) = NaN;
                qn(2*id) = NaN;
            end
            
            Bn = [Cn(~isnan(Cn(:,1)),:);Bn]; %#ok<AGROW>
            Qn = [qn(~isnan(qn));Qn]; %#ok<AGROW>
        end

        % Ajout de conditions supplementaires
        if size(Bn,1) < size(Bn,2)
            M = eye(size(Bn,2));
            bn = q(logical(kron((interf == i)',[1;1])));
            
            Bn = [M Bn' ; Bn zeros(size(Bn,1))]; %#ok<AGROW>
            Qn = [M*bn;Qn]; %#ok<AGROW>
        else
            Bn = Bn(1:size(Bn,2),:);
            Qn = Qn(1:size(Bn,2),:);
        end
        
        bn = Bn\Qn;

        b(logical(kron(interf' == i,[1;1]))) = bn(1:2*sum(patch_borders_ids));
    end
    
    %% Calcul des flux pour chaque interface
    detJ = @(J) sqrt(J(:,1).^2+J(:,2).^2);
    [Wg,Xg] = gaussPoints(domega);
    
    F = zeros(2*size(interf,2),size(interf,1));% stocke par [bx^1 by^1 bx^2 by^2 ...]' 1,2 l'id local des noeuds des interfaces pour chaque interface
    for i=1:size(interf,1)
        nodes_ids = interf(i,:);
        Xe = omega.nodes(nodes_ids,:);

        [M2,J] = shapesFunctions(domega,Xg,Xe,0);
        D = kron(diag(Wg.*detJ(J)),eye(2));

        F(:,i) = (M2'*D*M2)\b(:,i);
    end
    
    %% Calcul d'un nouveau maillage
    d = omega.order;
    if k_ord
        switch d + k_ord
            case 1
               a = [1 0;0 1;0 0]; 
            case 2
               a = [1 0;0 1;0 0;1/2 1/2;0 1/2;1/2 0]; 
            case 3
               a = [1 0;0 1;0 0;2/3 1/3;1/3 2/3;0 2/3;0 1/3;1/3 0;2/3 0;1/3 1/3];  
            case 4
               a = [1 0;0 1;0 0;3/4 1/4;1/2 1/2;1/4 3/4;0 3/4;0 1/2;0 1/4;1/4 0;1/2 0;3/4 0;1/4 1/4;1/2 1/4;1/4 1/2];   
            otherwise
                error('Degrée d''augmentation trop important');
        end

        M1 = shapesFunctions(omega,a,[],0);
        x = omega.nodes(:,1);
        x = M1(1:2:end,1:2:end)*x(omega.elems');
        y = omega.nodes(:,2);
        y = M1(1:2:end,1:2:end)*y(omega.elems');

        elems = reshape((1:numel(x))',[],omega.nbElems)';
        [nodes,~,ids] = uniquetol([x(:) y(:)],'ByRows',true);
        [node_ids,~,tmp_ids] = unique(ids,'stable');
        omega_adm = Mesh(nodes(node_ids,:),tmp_ids(elems));

        [interf_adm,eta] = omega_adm.interfaceElems;

        tmp_elems = omega.elems(:,1:3);
        inner_nodes = unique(tmp_elems(not(ismember(tmp_elems(:),domega.elems(:)))));
        if not(omega_adm.nbElems*(omega_adm.order+1)*(omega_adm.order+2)/2 - sum(sum(abs(eta),2)>1)*(omega_adm.order+1) + numel(inner_nodes) - omega_adm.nbNodes < eps)
            warning('Tous les noeuds doublons n''ont pas ete retire');
        end
    else
        omega_adm = omega;
        interf_adm = interf;
    end
    
    %% Interpolation du flux equilibre sur les nouvelles interfaces
    
    switch d + k_ord
        case 1
           c = [-1;1]; 
        case 2
           c = [-1;1;0]; 
        case 3
           c = [-1;1;-1/3;1/3];  
        case 4
           c = [-1;1;-1/2;0;1/2];   
        otherwise
           error('Degrée d''augmentation trop important');
    end

    M1 = shapesFunctions(domega,c,[],0);
    F_adm = M1*F; % stocke par [bx^1 by^1 bx^2 by^2 ...]' 1,2 l'id local des noeuds des interfaces pour chaque interface

    warning_error = warning('error','MATLAB:nearlySingularMatrix'); %#ok<CTPCT>
    %% Reconstruction du champ equilibre
    domega_adm = omega_adm.border;
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);
    detJb = @(J) sqrt(J(:,1).^2+J(:,2).^2);

    [Wg,Xg] = gaussPoints(omega_adm);
    [Wgb,Xgb] = gaussPoints(domega_adm,2*omega_adm.order);

    sigma_adm = zeros(numel(Wg)*3,omega_adm.nbElems);
    for i=1:omega_adm.nbElems
        nodes_ids = omega_adm.elems(i,:);
        Xe = omega_adm.nodes(nodes_ids,:);

        % Matrice de rigidite de l'element
        [M1,J] = shapesFunctions(omega_adm,Xg,Xe,1);
        D = kron(diag(Wg.*detJ(J)),B);
        K = M1'*D*M1;

        % Second membre
        Q = zeros(2*numel(nodes_ids),1);
        border_elems = find(eta(:,i));
        
        for j=1:numel(border_elems)
            border_nodes_ids = interf_adm(border_elems(j),:);
            Xe = omega_adm.nodes(border_nodes_ids,:);

            % retrouver les indices des noeuds localement a l'element
            [map,~] = find(bsxfun(@rdivide,nodes_ids(:),border_nodes_ids(:)') == 1);
            map = [(map(:)-1)*2+1 map(:)*2]';

            [N1,J] = shapesFunctions(domega_adm,Xgb,Xe,0);
            D = kron(diag(Wgb.*detJb(J)),eye(2));

            Q(map(:)) = Q(map(:)) + eta(border_elems(j),i)*N1'*D*(N1*F_adm(:,border_elems(j)));    
        end
        
        u = zeros(2*size(omega_adm.elems,2),1);
        % calcul le deplacement du probleme en bloquant les mvts rigide
        try
            u(2:end-2) = K(2:end-2,2:end-2)\Q(2:end-2);
        catch err
            if strcmpi(err.identifier,'MATLAB:nearlySingularMatrix')
                try
                    u(1:end-3) = K(1:end-3,1:end-3)\Q(1:end-3);
                catch err
                    if strcmpi(err.identifier,'MATLAB:nearlySingularMatrix')
                        warning(['Matrices singuliere pour l''element ' num2str(i) '.Resultats possiblement incorrect.']);
                    else
                        rethrow(err);
                    end
                end
            else
                warning(warning_error);
                rethrow(err);
            end
        end

        % calcul le champ equilibre resultant (aux points de gauss)
        D = kron(eye(numel(Wg)),B);
        sigma_adm(:,i) = D*M1*u;
    end
    sigma_adm = sigma_adm(:);
    warning(warning_error);
end