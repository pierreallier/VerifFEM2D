function h = plotHighOrderStress(omega,stress)

    % TODO quelque verifs
    stress = stress(:);
    if omega.order == 1
        h = omega.plotElemField(stress);
    else
        s11 = reshape(stress(1:3:end),[],omega.nbElems)';
        s22 = reshape(stress(2:3:end),[],omega.nbElems)';
        s12 = reshape(stress(3:3:end),[],omega.nbElems)';
        switch omega.order
            case 2
                connec = omega.elems(:,[1 4 6 4 5 6 4 2 5 6 5 3]);
                sindex = [1 1 2 3];
            case 3
                connec = omega.elems(:,[1 4 9 4 10 9 4 5 10 5 6 10 5 2 6 10 6 7 9 10 8 10 7 8 8 7 3]);
                sindex = [4 4 5 5 6 1 3 3 2];
            case 4
                connec = omega.elems(:,[1 4 12 4 13 12 4 5 13 5 14 13 5 6 14 6 7 14 6 2 7 12 13 11 13 ...
                                        15 11 13 14 15 14 8 15 14 7 8 11 15 10 15 9 10 15 8 9 10 9 3]);
                sindex = [4 4 11 11 8 8 6 9 9 1 3 10 12 12 7 5];
            otherwise
                error('Maillage d''ordre trop eleve, abandon');
        end
        new_omega = Mesh(omega.nodes,reshape(connec',3,[])');
        
        s11 = s11(:,sindex);
        s22 = s22(:,sindex);
        s12 = s12(:,sindex);
        new_stress = reshape([reshape(s11',[],1) reshape(s22',[],1) reshape(s12',[],1)]',[],1);
        h = new_omega.plotElemField(new_stress);
    end
end