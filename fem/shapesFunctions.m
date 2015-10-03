function [N,J] = shapesFunctions(omega,a,x,od)
% Evalue les fonctions de forme (ou ses dérivées) aux points de quadratures
% a sur l'éléments dont les coordonnées des noeuds sont dans x, connaissant
% le maillage support.

    type = omega.type;
    order = omega.order;
    
    % Coefficients polynomiaux (pour les éléments d'ordre 3-4 triangulaires)
    if type == 2 && order > 2
        C1 = [1 0 0 0 0 0 0 0 0 0;-5.5 1 0 9 -4.5 0 0 0 0 0;-5.5 0 1 0 0 0 0 -4.5 9 0;9 -4.5 0 -22.5 18 0 0 0 0 0;9 0 -4.5 0 0 0 0 18 -22.5 0;18 0 0 -22.5 4.5 -4.5 -4.5 4.5 -22.5 27;-4.5 4.5 0 13.5 -13.5 0 0 0 0 0;-4.5 0 4.5 0 0 0 0 -13.5 13.5 0;-13.5 0 0 27 -13.5 13.5 0 0 13.5 -27;-13.5 0 0 13.5 0 0 13.5 -13.5 27 -27];
        C2 = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;-25/3 -1 0 16 -12 16/3 0 0 0 0 0 0 0 0 0;-25/3 0 -1 0 0 0 0 0 0 16/3 -12 16 0 0 0;70/3 22/3 0 -208/3 76 -112/3 0 0 0 0 0 0 0 0 0;70/3 0 22/3 0 0 0 0 0 0 -112/3 76 -208/3 0 0 0;140/3 0 0 -208/3 28 -16/3 16/3 4 16/3 -16/3 28 -208/3 96 -32 -32; ...
              -80/3 -16 0 96 -128 224/3 0 0 0 0 0 0 0 0 0;-80/3 0 -16 0 0 0 0 0 0 224/3 -128 96 0 0 0;-80 0 0 192 -144 32 -32 -16 0 0 -16 96 -224 160 32;-80 0 0 96 -16 0 0 -16 -32 32 -144 192 -224 32 160;32/3 32/3 0 -128/3 64 -42.667 0 0 0 0 0 0 0 0 0;32/3 0 32/3 0 0 0 0 0 0 -128/3 64 -128/3 0 0 0; ...
              128/3 0 0 -128 128 -128/3 128/3 0 0 0 0 -128/3 128 -128 0;64 0 0 -128 64 0 0 64 0 0 64 -128 256 -128 -128;128/3 0 0 -128/3 0 0 0 0 128/3 -128/3 128 -128 128 0 -128];
    end
    % Fonctions de forme
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
                N  = [a(:,1).*(2*a(:,1)-1) a(:,2).*(2*a(:,2)-1) ...
                      (1-a(:,1)-a(:,2)).*(1-2*a(:,1)-2*a(:,2))  4*a(:,1).*a(:,2) ...
                      4*a(:,2).*(1-a(:,1)-a(:,2))  4*a(:,1).*(1-a(:,1)-a(:,2))];
            case 23 % Triangle order 3
                N = [1+0*a(:,1) a a.^2 prod(a,2) a.^3 prod(a,2).*a(:,1) prod(a,2).*a(:,2)]*C1;
            case 24 % Triangle order 4
                N = [1+0*a(:,1) a a.^2 prod(a,2) a.^3 prod(a,2).*a(:,1) prod(a,2).*a(:,2) a.^4 a(:,1).^3.*a(:,2) prod(a,2).^2 a(:,1).*a(:,2).^3]*C2;
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
                DN = cat(2,[0*a(:,1) 1+0*a(:,1) 0*a(:,2) 2*a(:,1) 0*a(:,2) a(:,2) 3*a(:,1).^2 0*a(:,2) 2*prod(a,2) a(:,2).^2]*C1, ...
                           [0*a(:,1) 0*a(:,1) 1+0*a(:,2) 0*a(:,1) 2*a(:,2) a(:,1) 0*a(:,2) 3*a(:,2).^2 a(:,1).^2 2*prod(a,2)]*C1);
            case 24 % Triangle order 4
                DN = cat(2,[0*a(:,1) 1+0*a(:,1) 0*a(:,2) 2*a(:,1) 0*a(:,2) a(:,2) 3*a(:,1).^2 0*a(:,2) 2*prod(a,2) a(:,2).^2 4*a(:,1).^3 0*a(:,2) 3*prod(a,2).*a(:,1) 2*a(:,2)*prod(a,2) a(:,2).^3]*C2, ...
                           [0*a(:,1) 0*a(:,1) 1+0*a(:,2) 0*a(:,1) 2*a(:,2) a(:,1) 0*a(:,1) 3*a(:,2).^2 a(:,1).^2 2*prod(a,2) 0*a(:,1) 4*a(:,2).^3 a(:,1).^3 2*a(:,1)*prod(a,2) 3*prod(a,2).*a(:,2)]*C2);
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