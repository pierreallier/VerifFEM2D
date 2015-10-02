function [N,J] = shapesFunctions(omega,a,x,od)
% Evalue les fonctions de forme (ou ses dérivées) aux points de quadratures
% a sur l'éléments dont les coordonnées des noeuds sont dans x, connaissant
% le maillage support.

    type = omega.type;
    order = omega.order;

    % Shapes functions
    if od == 0
        switch type*10+order
            case 1 % Node
                N  = 1+0*a;
            case 11 % Tress order 1
                N  = [(1-a)/2  (1+a)/2];
            case 12 % Tress order 2
                N  = [a.*(a-1)/2  a.*(1+a)/2 1-a.^2];
            case 13 % Tress order 3
                N  = [(1-a).*(9*a.^2-1)  (a+1).*(9*a.^2-1)  9*(1-3*a).*(1-a.^2)  9*(1+3*a).*(1-a.^2)]/16;
            case 14 % Tress order 4
                N  = [a.*(a-1).*(4*a.^2-1)/6  a.*(a+1).*(4*a.^2-1)/6  4*a.*(2*a-1).*(1-a.^2)/3  (a.^2-1).*(4*a.^2-1)  4*a.*(2*a+1).*(1-a.^2)/3];
            case 21 % Triangle order 1
                N  = [a(:,1)  a(:,2)  1-a(:,1)-a(:,2)];
            case 22 % Triangle order 2
                N  = [a(:,1).*(2*a(:,1)-1)  a(:,2).*(2*a(:,2)-1)  (1-a(:,1)-a(:,2)).*(1-2*a(:,1)-2*a(:,2))  4*a(:,1).*a(:,2)  4*a(:,2).*(1-a(:,1)-a(:,2))  4*a(:,1).*(1-a(:,1)-a(:,2))];
            case 23 % Triangle order 3

            case 24 % Triangle order 4

            otherwise
                error('Element inconnu');
        end
        N = kron(N,eye(2));
    end
    
    % Derivative of shape function
    if od == 1 || nargout > 1
        switch type*10+order
            case 1 % Node
                DN = 0*a;
            case 11 % Tress order 1
                DN = [-1/2+0*a  1/2+0*a];
            case 12 % Tress order 2
                DN = [a-1/2 a+1/2  -2*a];
            case 13 % Tress order 3
                DN = [1-27*a^2+18*a  27*a.^2+18*a-1  81*a.^2-18*a-27  27-81*a.^2-18*a]/16;
            case 14 % Tress order 4
                DN = [(4*a-1).*(4*a.^2-2*a-1)/6  (4*a+1).*(4*a.^2+2*a-1)/6  (-32*a.^3+12*a.^2+16*a-4)/3  2*a.*(8*a.^-5)  (-32*a.^3-12*a.^2+16*a+4)/3];
            case 21 % Triangle order 1
                DN = [1+0*a(:,1)  0*a(:,1)  -1+0*a(:,1)  ...
                      0*a(:,2)  1+0*a(:,2)  -1+0*a(:,2)];
            case 22 % Triangle order 2
                DN = [4*a(:,1)-1  0*a(:,1)  4*(a(:,1)+a(:,2)-1)+1  4*a(:,2)  -4*a(:,2)  4*(1-2*a(:,1)-a(:,2))  ...
                      0*a(:,2)  4*a(:,2)-1  4*(a(:,1)+a(:,2)-1)+1  4*a(:,1)  4*(1-a(:,1)-2*a(:,2))  -4*a(:,1)];
            case 23 % Triangle order 3

            case 24 % Triangle order 4

            otherwise
                error('Element inconnu');
        end
        if type == 2
            DN = reshape(DN',[],size(a,1)*2)';
        end
        J = DN*x;
        
        if od == 1
            M = mat2cell(J,size(J,2)*ones(size(J,1)/size(J,2),1),size(J,2));
            J = blkdiag(M{:});
            N = J\DN;
            N = kron(N(1:2:end,:),[1 0;0 0;0 1]) + kron(N(2:2:end,:),[0 0;0 1;1 0]);
        end
    end
end