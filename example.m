%% Script d'exemple de resolution du probleme elements finis
addpath('./fem/');
addpath('./eet/');

% Parametres du probleme
E = 200;
nu = 0.3;
B = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2]; % Matrix de Hook (Avec les notations de Voigt)
Fd = 1000*[0 1];
a = 1;

% Charge le maillage
omega = Mesh(geo2msh('./meshes/plate_hole.geo'));

% Extraction de maillage secondaires
H = max(omega.nodes(:,2));
domega = omega.border;
dtop = domega.restrict(@(x) x(:,2) == H);
dleft = domega.restrict(@(x) x(:,1) == 0);
dbottom = domega.restrict(@(x) x(:,2) == 0 & x(:,1) >= a);

% Ecriture du systeme
K = FEMMat(omega,B);
F = FEMVec(dtop,Fd);

% Initialisation de l'inconnue
u = zeros(size(K,2),1);

% Imposition des conditions aux limites
[cl_index,u0] = CL(dleft,0,[1 0],dbottom,0,[0 1]);
u(cl_index) = u0;
   
% Resolution
u(~cl_index) = K(~cl_index,~cl_index)\(F(~cl_index) - K(~cl_index,cl_index)*u(cl_index));
sigma = stress(omega,B,u);

% Visualisation
figure('Name','Solution');
  subplot(1,2,1)
    plotElemField(deform(omega,u,1./max(abs(u))),sigma);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Stresses');
  subplot(1,2,2);
    plot(omega.border);
    hold on
    quiver(omega.nodes(:,1),omega.nodes(:,2),u(1:2:end),u(2:2:end));
    xlabel('x');
    ylabel('y');
    title('Displacement');
    
[omega_adm,sigma_adm] = eet(omega,sigma,B,dtop,Fd,dleft,dbottom,1);
u_adm = interpField(omega,u,omega_adm);

n = max(u);
figure('Name','Interpolation');
  subplot(1,2,1);
    plot(omega.border);
    hold on
    quiver(omega.nodes(:,1),omega.nodes(:,2),u(1:2:end)./n,u(2:2:end)./n,2);
    xlabel('x');
    ylabel('y');
    title('Displacement');
  subplot(1,2,2);
    plot(omega_adm.border);
    hold on
    quiver(omega_adm.nodes(:,1),omega_adm.nodes(:,2),u_adm(1:2:end)./n,u_adm(2:2:end)./n,0);
    xlabel('x');
    ylabel('y');
    title('Interpolated Displacement');

figure('Name','CRE');
  subplot(1,2,1);
    plotElemField(deform(omega,u,1./max(abs(u))),sigma);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Stresses');
  subplot(1,2,2);
    plotAdmField(deform(omega_adm,u_adm,1./max(abs(u_adm))),sigma_adm);
    xlabel('x');
    ylabel('y');
    colorbar;
    title('Admissible Stresses');