function [N,J] = shapesFunctions(omega,a,x,od)
% Evalue les fonctions de forme (ou ses derivees) aux points de quadratures
% a sur l'elements dont les coordonnees des noeuds sont dans x, connaissant
% le maillage support.

    type = omega.type;
    order = omega.order;
    
    % Coefficients polynomiaux (pour les elements d'ordre 3-4 triangulaires)
    if type == 2 && order > 2
        C1 = [1 0 0 0 0 0 0 0 0 0;-5.5 1 0 9 -4.5 0 0 0 0 0;-5.5 0 1 0 0 0 0 -4.5 9 0;9 -4.5 0 -22.5 18 0 0 0 0 0;9 0 -4.5 0 0 0 0 18 -22.5 0;18 0 0 -22.5 4.5 -4.5 -4.5 4.5 -22.5 27;-4.5 4.5 0 13.5 -13.5 0 0 0 0 0;-4.5 0 4.5 0 0 0 0 -13.5 13.5 0;-13.5 0 0 27 -13.5 13.5 0 0 13.5 -27;-13.5 0 0 13.5 0 0 13.5 -13.5 27 -27];
        C2 = 1/3*[3 0 0 0 0 0 0 0 0 0 0 0 0 0 0; ...
                  -25 -3 0 48 -36 16 0 0 0 0 0 0 0 0 0; ...
                  -25 0 -3 0 0 0 0 0 0 16 -36 48 0 0 0; ...
                  70 22 0 -208 228 -112 0 0 0 0 0 0 0 0 0; ...
                  70 0 22 0 0 0 0 0 0 -112 228 -208 0 0 0; ...
                  140 0 0 -208 84 -16 16 12 16 -16 84 -208 288 -96 -96; ...
                  -80 -48 0 288 -384 224 0 0 0 0 0 0 0 0 0; ...
                  -80 0 -48 0 0 0 0 0 0 224 -384 288 0 0 0; ...
                  -240 0 0 576 -432 96 -96 -48 0 0 -48 288 -672 480 96; ...
                  -240 0 0 288 -48 0 0 -48 -96 96 -432 576 -672 96 480; ...
                  32 32 0 -128 192 -128 0 0 0 0 0 0 0 0 0; ...
                  32 0 32 0 0 0 0 0 0 -128 192 -128 0 0 0; ...
                  128 0 0 -384 384 -128 128 0 0 0 0 -128 384 -384 0; ...
                  192 0 0 -384 192 0 0 192 0 0 192 -384 768 -384 -384; ...
                  128 0 0 -128 0 0 0 0 128 -128 384 -384 384 0 -384];
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
                %N  = [a.*(1-a).*(1-4*a.^2) a.*(1+a).*(4*a.^2-1) 8*a.*(a.^2-1).*(1-2*a) 6*(1+a.^2.*(4*a.^2-5)) 8*a.*(1-a.^2).*(1+2*a)]/6;
                N = bsxfun(@power,a(:),0:order)*[0 0 0 6 0;1 -1 -8 0 8;-1 -1 16 -30 16; -4 4 8 0 -8;4 4 -16 24 -16]/6;
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
                DN = [1-27*a.^2+18*a  27*a.^2+18*a-1  81*a.^2-18*a-27  27-81*a.^2-18*a]/16;
            case 14 % Tress order 4
                DN = [0*a(:) 1+0*a(:) 2*a(:) 3*a(:).^2 4*a(:).^3]*[0 0 0 6 0;1 -1 -8 0 8;-1 -1 16 -30 16; -4 4 8 0 -8;4 4 -16 24 -16]/6;
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
                DN = cat(2,[0*a(:,1) 1+0*a(:,1) 0*a(:,2) 2*a(:,1) 0*a(:,2) a(:,2) 3*a(:,1).^2 0*a(:,2) 2*prod(a,2) a(:,2).^2 4*a(:,1).^3 0*a(:,2) 3*prod(a,2).*a(:,1) 2*a(:,2).*prod(a,2) a(:,2).^3]*C2, ...
                           [0*a(:,1) 0*a(:,1) 1+0*a(:,2) 0*a(:,1) 2*a(:,2) a(:,1) 0*a(:,1) 3*a(:,2).^2 a(:,1).^2 2*prod(a,2) 0*a(:,1) 4*a(:,2).^3 a(:,1).^3 2*a(:,1).*prod(a,2) 3*prod(a,2).*a(:,2)]*C2);
            otherwise
                error('Element inconnu');
        end
        if type == 2
            DN = reshape(DN',[],size(a,1)*2)';
        end
        J = DN*x;
        
        if od == 1
            M = mat2cell(J,size(J,2)*ones(size(J,1)/size(J,2),1),size(J,2));
            J2 = blkdiag(M{:});
            N = J2\DN;
            N = kron(N(1:2:end,:),[1 0;0 0;0 1]) + kron(N(2:2:end,:),[0 0;0 1;1 0]);
        end
    end
end