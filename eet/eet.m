function [omega_adm,sigma_adm] = eet(omega,sigma,B,dfomega,Fd,k_ord)

    %% Initialisation
    % Récupération des interfaces
    domega = omega.border;
    [interf,eta,n] = omega.interfaceElems;
    
    % Recuperation du bord où le chargement est impose
    loading_nodes_ids = unique(dfomega.elems(:));
    border_nodes_ids = unique(omega.border.elems(:));
    clu_nodes_ids = setdiff(border_nodes_ids,loading_nodes_ids);
    loading_interf_ids = ismember(sort(interf,2),sort(dfomega.elems,2),'rows');

    %% Calcul de la condition de prolongement \int_{E} \omega_h \grad(\phi_i) dx
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);
    [Wg,Xg] = gaussPoints(omega);
    Ng = numel(Wg);

    Q = zeros(size(omega.elems,2)*2,omega.nbElems); % stocké par [Fx^1 Fy^1 Fx^2 Fy^2 Fx^3 Fy^3 ...]' 1,2,3 l'id local des noeuds des elements pour chaque element
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
    end
    
    %% Calcul des flux equilibre pour chaque noeud
    interf_nodes = unique(interf(:));
    
    b = zeros(2*size(interf,2),size(interf,1)); % stocké par [bx^1 by^1 bx^2 by^2 ...]' 1,2 l'id local des noeuds des interfaces pour chaque interface
    for k=1:numel(interf_nodes)
        i = interf_nodes(k);
        patch_elems_ids = any(omega.elems == i,2); % elements qui partagent ce noeud
        patch_borders_ids = any(interf == i,2); % bords qui partagent ce noeud

        Bn = kron(eta(patch_borders_ids,patch_elems_ids)',eye(2));
        Qn = Q(logical(kron((omega.elems == i)',[1;1])));
        
        cond_prolong = sum([Qn(1:2:end) Qn(2:2:end)],1);
        % verifications condtions noeud interieur
%         if all(border_nodes_ids ~= i)
%             disp(['Condition sur le noeud ' num2str(i) ' interieur : ' num2str(cond_prolong)]);
%         end
        
        % On retire les équations de compatibilité si noeud interieur/bord du chargement
        if all(clu_nodes_ids ~= i)
            Bn = Bn(1:end-2,:);
            Qn = Qn(1:end-2);
        end

        % Ajout des conditions de chargement
        if any(loading_nodes_ids == i)
            Cn = diag(loading_interf_ids(patch_borders_ids));
            Cn = kron(Cn(any(Cn,2),:),eye(2));

            qn = q(logical(kron((interf == i)',[1;1])));
            qn = qn(logical(kron(loading_interf_ids(patch_borders_ids),[1;1])));
            
            % verification conditions
            cond_prolong = cond_prolong - sum([qn(1:2:end) qn(2:2:end)],1);
%             disp(['Condition sur le noeud ' num2str(i) ' du bord : ' num2str(cond_prolong)]);

            Bn = [Cn;Bn]; %#ok<AGROW>
            Qn = [qn;Qn]; %#ok<AGROW>
        end

        % Ajout de conditions supplementaires
        if size(Bn,1) < size(Bn,2)
            M = eye(size(Bn,2));
            bn = q(logical(kron((interf == i)',[1;1])));

            Bn = [M Bn' ; Bn zeros(size(Bn,1))]; %#ok<AGROW>
            Qn = [M*bn;Qn]; %#ok<AGROW>
        else
            Bn = Bn(1:size(Bn,2),:);
            Qn = Qn(1:size(Bn,2));
        end

        % Resolution
        bn = Bn\Qn;
        b(logical(kron(interf' == i,[1;1]))) = bn(1:2*sum(patch_borders_ids));
        
    end
    
    %% Calcul des flux pour chaque interface
    detJ = @(J) sqrt(J(:,1).^2+J(:,2).^2);
    [Wg,Xg] = gaussPoints(domega);
    
    F = zeros(2*size(interf,2),size(interf,1));% stocké par [bx^1 by^1 bx^2 by^2 ...]' 1,2 l'id local des noeuds des interfaces pour chaque interface
    for i=1:size(interf,1)
        nodes_ids = interf(i,:);
        Xe = omega.nodes(nodes_ids,:);

        [M2,J] = shapesFunctions(domega,Xg,Xe,0);
        D = kron(diag(Wg.*detJ(J)),eye(2));

        F(:,i) = (M2'*D*M2)\b(:,i);
    end
    
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
    
    %% % Interpolation du flux equilibre sur les nouvelles interfaces
    
    switch d + k_ord
        case 1
           c = [-1;1]; 
        case 2
           c = [-1;1;0]; 
        case 3
           c = [-1;1;-2/3;2/3];  
        case 4
           c = [-1;1;-1/2;0;1/2];   
        otherwise
           error('Degrée d''augmentation trop important');
    end

    M1 = shapesFunctions(domega,c,[],0);
    F_adm = M1*F; % stocké par [bx^1 by^1 bx^2 by^2 ...]' 1,2 l'id local des noeuds des interfaces pour chaque interface

    warning_error = warning('error','MATLAB:nearlySingularMatrix'); %#ok<CTPCT>
    %% Reconstruction du champ equilibre
    domega_adm = omega_adm.border;
    detJ = @(J) J(1:2:end,1).*J(2:2:end,2) - J(2:2:end,1).*J(1:2:end,2);
    detJb = @(J) sqrt(J(:,1).^2+J(:,2).^2);

    [Wg,Xg] = gaussPoints(omega_adm);
    [Wgb,Xgb] = gaussPoints(domega_adm);

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

            % retrouver les indices des noeuds localement à l'element
            [map,~] = find(bsxfun(@rdivide,nodes_ids(:),border_nodes_ids(:)') == 1);
            map = [(map(:)-1)*2+1 map(:)*2]';

            [N1,J] = shapesFunctions(domega_adm,Xgb,Xe,0);
            D = kron(diag(Wgb.*detJb(J)),eye(2));

            Q(map(:)) = Q(map(:)) + eta(border_elems(j),i)*N1'*D*(N1*F_adm(:,border_elems(j)));    
        end
        
        % Verification
%         eq = [sum(Q(1:2:end)) sum(Q(2:2:end))];
%         if any(eq > 10^-10)
%             disp(['Equilibre non verifié sur l''elem ' num2str(i) ' : ' num2str(eq)]);
%         end
        
        % calcul le deplacement du probleme en bloquant les mvts rigide
        try
            u = zeros(2*size(omega_adm.elems,2),1);
            u(2:end-2) = K(2:end-2,2:end-2)\Q(2:end-2);
        catch err
            if strcmpi(err.identifier,'MATLAB:nearlySingularMatrix')
                warning(['Matrices singulière pour l''element ' num2str(i) '.Resultats possiblement incorrect.']);
            else
                warning(warning_error);
                rethrow(err);
            end
        end

        % calcul le champ equilibre resultant (au points de gauss)
        D = kron(eye(numel(Wg)),B);
        sigma_adm(:,i) = D*M1*u;
    end
    sigma_adm = sigma_adm(:);
    warning(warning_error);
end