function [Wg,Xg] = gaussPoints(omega)
% Calcule les poids Wg et les points Xg pour la m�thode de quadrature de Gauss
% (int�gration) � partir d'un type d'�l�ment sp�cifique (donn�e par omega).

    switch omega.type
        case 2
            xw = triangularQuad(omega.order);
            Wg = xw(:,end);
            Xg = xw(:,1:end-1);
        case 1
            xw = linearQuad(omega.order);
            Wg = xw(:,end);
            Xg = xw(:,1:end-1);
        otherwise
            Wg = 1;
            Xg = 0;
    end
end

function xw = linearQuad(order)
%linearQuad Computes the Legendre-Gauss nodes and weights on an
% interval [-1,1] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [-1,1]
% which you can evaluate at any x in [-1,1]. Simply evaluate it at
% all of the values contained in the xw(:,1) vector to obtain a
% vector f. Then compute the definite integral using sum(f.*xw(:,3));
%
% xw is a nx2 matrix where the 1st two columns are coordinates of
% points and third one is the weights.
%
% Values from Analyse des solides d�formables par la MEF (M. Bonnet)
%
    N = 2*ceil((order+1)/2);
    switch N
        case 1
            xw=[0 2];
        case 2
            xw=[-0.57735026918962573106 1
                0.57735026918962573106 1];
        case 3
            xw=[-0.77459666924148340428 0.55555555555555558023
                0                      0.88888888888888883955
                0.77459666924148340428 0.55555555555555558023];
        case 4
            xw=[-0.86113631159405257254 0.34785484513745384971
                -0.33998104358485631282 0.65214515486254620580
                0.33998104358485631282 0.65214515486254620580
                0.86113631159405257254 0.34785484513745384971];
        case 5
            xw=[-0.90617984593866399279 0.23692688505618908751
                -0.53846931010568309103 0.47862867049936646804
                0                       0.56888888888888888889
                0.53846931010568309103  0.47862867049936646804
                0.90617984593866399279  0.23692688505618908751];
        case 6
            xw=[-0.93246951420315202781 0.17132449237917034504
                -0.66120938646626451366 0.36076157304813860756
                -0.23861918608319690863 0.46791393457269104738
                0.23861918608319690863  0.46791393457269104738
                0.66120938646626451366  0.36076157304813860756
                0.93246951420315202781  0.17132449237917034504];
        otherwise
            % An algorithm proposed by Greg von Winckel - 02/25/2004
            N=N-1;
            N1=N+1; N2=N+2;
            xu=linspace(-1,1,N1)';

            % Initial guess
            x=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

            % Compute the zeros of the N+1 Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method
            L=zeros(N1,N2); % Legendre-Gauss Vandermonde Matrix
            Lp=zeros(N1,N2);% Derivative of LGVM
            x0=2;
            % Iterate until new points are uniformly within epsilon of old points
            while max(abs(x-x0))>eps
                L(:,1)=1;
                Lp(:,1)=0;
                L(:,2)=x;
                Lp(:,2)=1;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*x.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-x.*L(:,N2) )./(1-x.^2);
                x0=x;
                x=x0-L(:,N2)./Lp;
            end
            % Compute the weights
            xw=[x 2./((1-x.^2).*Lp.^2)*(N2/N1)^2];
    end
end

function xw = triangularQuad(N)
%triangularQuad Computes the Legendre-Gauss nodes and weights on an
% triangular interval [0,1]x[0,1-x] with truncation order N
%
% Suppose you have a continuous function f(x,y) which is defined on
% [0,1]x[0,1-x] which you can evaluate at any x,y in [0,1]x[0,1-x].
% Simply evaluate it at all of the values contained in the xw(:,1:2) vector
% to obtain a vector f. Then compute the definite integral using
% sum(f.*xw(:,3));
%
% xw is a nx3 matrix where the 1st two columns are coordinates of
% points and third one is the weights.
%
% Values from http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF
%
    switch N
        case 1
            xw=[0.33333333333333 0.33333333333333 1.00000000000000];
        case 2
            xw=[0.16666666666667 0.16666666666667 0.16666666666667
                0.16666666666667 0.66666666666667 0.16666666666667
                0.66666666666667 0.16666666666667 0.16666666666667];
        case 3
            xw=[0.33333333333333 0.33333333333333 -0.56250000000000
                0.20000000000000 0.20000000000000 0.52083333333333
                0.20000000000000 0.60000000000000 0.52083333333333
                0.60000000000000 0.20000000000000 0.52083333333333];
        case 4
            xw=[0.44594849091597 0.44594849091597 0.22338158967801
                0.44594849091597 0.10810301816807 0.22338158967801
                0.10810301816807 0.44594849091597 0.22338158967801
                0.09157621350977 0.09157621350977 0.10995174365532
                0.09157621350977 0.81684757298046 0.10995174365532
                0.81684757298046 0.09157621350977 0.10995174365532];
        case 5
            xw=[0.33333333333333 0.33333333333333 0.22500000000000
                0.47014206410511 0.47014206410511 0.13239415278851
                0.47014206410511 0.05971587178977 0.13239415278851
                0.05971587178977 0.47014206410511 0.13239415278851
                0.10128650732346 0.10128650732346 0.12593918054483
                0.10128650732346 0.79742698535309 0.12593918054483
                0.79742698535309 0.10128650732346 0.12593918054483];
        case 6
            xw=[0.24928674517091 0.24928674517091 0.11678627572638
                0.24928674517091 0.50142650965818 0.11678627572638
                0.50142650965818 0.24928674517091 0.11678627572638
                0.06308901449150 0.06308901449150 0.05084490637021
                0.06308901449150 0.87382197101700 0.05084490637021
                0.87382197101700 0.06308901449150 0.05084490637021
                0.31035245103378 0.63650249912140 0.08285107561837
                0.63650249912140 0.05314504984482 0.08285107561837
                0.05314504984482 0.31035245103378 0.08285107561837
                0.63650249912140 0.31035245103378 0.08285107561837
                0.31035245103378 0.05314504984482 0.08285107561837
                0.05314504984482 0.63650249912140 0.08285107561837];
        case 7
            xw=[0.33333333333333 0.33333333333333 -0.14957004446768
                0.26034596607904 0.26034596607904 0.17561525743321
                0.26034596607904 0.47930806784192 0.17561525743321
                0.47930806784192 0.26034596607904 0.17561525743321
                0.06513010290222 0.06513010290222 0.05334723560884
                0.06513010290222 0.86973979419557 0.05334723560884
                0.86973979419557 0.06513010290222 0.05334723560884
                0.31286549600487 0.63844418856981 0.07711376089026
                0.63844418856981 0.04869031542532 0.07711376089026
                0.04869031542532 0.31286549600487 0.07711376089026
                0.63844418856981 0.31286549600487 0.07711376089026
                0.31286549600487 0.04869031542532 0.07711376089026
                0.04869031542532 0.63844418856981 0.07711376089026];
        otherwise
            xw=[0.33333333333333 0.33333333333333 0.14431560767779
                0.45929258829272 0.45929258829272 0.09509163426728
                0.45929258829272 0.08141482341455 0.09509163426728
                0.08141482341455 0.45929258829272 0.09509163426728
                0.17056930775176 0.17056930775176 0.10321737053472
                0.17056930775176 0.65886138449648 0.10321737053472
                0.65886138449648 0.17056930775176 0.10321737053472
                0.05054722831703 0.05054722831703 0.03245849762320
                0.05054722831703 0.89890554336594 0.03245849762320
                0.89890554336594 0.05054722831703 0.03245849762320
                0.26311282963464 0.72849239295540 0.02723031417443
                0.72849239295540 0.00839477740996 0.02723031417443
                0.00839477740996 0.26311282963464 0.02723031417443
                0.72849239295540 0.26311282963464 0.02723031417443
                0.26311282963464 0.00839477740996 0.02723031417443
                0.00839477740996 0.72849239295540 0.02723031417443];
    end
    xw(:,end) = xw(:,end)/2;
end
