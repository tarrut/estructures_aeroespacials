%-------------------------------------------------------------------------%
% ASSIGNMENT 03 - (A)
%-------------------------------------------------------------------------%
% Date: April, 2022
% Author/s: Albert Bosch i Pau Tarrés
% 
% Code developed during the "Estructures Aeroespacials" course of GRETA 
% degree at UPC's ESEIAAT. The code implements finite elements to compute 
% the displacements, shear stresses and moments of a plane's wing.

clear;
close all;

%% INPUT DATA

% Material properties
E = 85e9;

% Cross-section parameters
t1 = 1.5e-3;
t2 = 4e-3;
h1 = 0.5;
h2 = 0.25;
b = 0.775;


% Other data
g = 9.81;
L1 = 5;
L2 = 10;
L = L1+L2;
Me = 2550;
M = 35000;

% Number of elements for each part
n_d = 1;
n_el = 30;
n_nod = 31;
n_ne = 2;
n_i = 2;
n_dof = n_nod*n_i;
n_dof_el = n_ne*n_i;


%% PRECOMPUTATIONS

% Compute section: 
% A  - Section area
alpha = atan((h2-h1)/(2*b));

A = t2*(h1+h2)+b/cos(alpha)*2*t1;

z_c = -h2*b/(h1+h2);
y_c = h1/2;

% Iz - Section inertia
Iz1 = 1/12*h1^3*t2;
Iz2 = 1/12*h2^3*t2;
Iz3 = 1/12*t1*(b/cos(alpha))^3*sin(alpha)^2+t1*b/cos(alpha)*(h2/2+(h1-h2)/4)^2;

Iz = Iz1 + Iz2 + 2*Iz3;

% Compute parameter l:
% l - Equilibrium parameter
l_s = sym('l_s');
x = sym('x');
lam1 = -(M/(4*(L1+L2))+3*M*(L1-x)/(2*L2^2))*g;
lam2 = -(M/(4*(L1+L2)))*g;
q1 = l_s*(0.8-0.2*cos(pi*x/L1));
q2 = l_s*(1-(x-L1)/L2)*(1+(x-L1)/L2);

int_lam1 = int(lam1,x,0,L1);
int_lam2 = int(lam2,x,L1,L1+L2);
massa = -(int_lam1+int_lam2)/g+Me;
int_q1 = int(q1,x,0,L1);
int_q2 = int(q2,x,L1,L1+L2);

eq = int_lam1+int_lam2+int_q1+int_q2-Me*g;
l_d = double(solve(eq));
q_t_1 = lam1 + subs(q1,l_s,l_d);
q_t_2 = lam2 + subs(q2,l_s,l_d);

q_t = piecewise(x<L1,q_t_1, x>=L1,q_t_2);
fplot(q_t)

%%
X = zeros(n_nod,n_d);
for n = 1:n_nod
   X(n) = (L1+L2)/n_el*(n-1); 
end

Tn = zeros(n_el, n_ne);
for e=1:n_el
    for ne = 1:n_ne
       Tn(e,ne) = e+ne-1; 
    end
end

Td = zeros(n_el, n_dof_el);
for element = 1:n_el
    for nod = 1:n_ne
        for dof = 1:n_i          
            Td(element,(nod-1)*n_i+dof) = (Tn(element,nod)-1)*n_i+dof;
        end
    end
end

Tmat = zeros(n_el, 1);
for e = 1:n_el
   Tmat(e) = 1;
end

mat = [% Young M.        Section A.    Inertia 
              E,                A,         Iz;  % Material (1)
    ];


%%Compute Kel

Kel = zeros(4, 4, n_el);
for e = 1:n_el

    x_e = zeros(2, n_d);

    for i = 1:2
        for j = 1:n_d
            x_e(i,j) = X(Tn(e,i),j);
            if i==2
                x_e(i,j) = -x_e(i,j);
            end
        end
    end

    l = sqrt(sum(sum(x_e).^2));

    R = [12, 6*l, -12, 6*l;
        6*l, 4*l^2, -6*l, 2*l^2;
        -12, -6*l, 12, -6*l;
        6*l, 2*l^2, -6*l, 4*l^2];

     Kel(:, :, e) = (mat(Tmat(e),3)*mat(Tmat(e),1)/l^3).*R;
    
end

%% Assembly KG
KG = zeros(n_dof);
    
for e = 1:n_el
    for i = 1:n_dof_el
        I = Td(e, i);
        for j = 1:n_dof_el
            J = Td(e, j);
            KG(I, J) = KG(I, J) + Kel(i, j, e);
        end
    end
end

%% Force Vector
F_el = zeros(n_dof_el,n_el);
for e = 1:n_el
    x_e = zeros(2, n_d);

    for i = 1:2
        for j = 1:n_d
            x_e(i,j) = X(Tn(e,i),j);
            if i==2
                x_e(i,j) = x_e(i,j);
            end
        end
    end
    
    l = sqrt(sum(sum(x_e).^2));
    
    q_m = int(q_t, x_e(1,1),x_e(2,1))/l;
    
    F_e = q_m*l/2.*[1;l/6;1;-l/6];
    for r = 1:n_dof_el
       F_el(r,e) = F_el(r,e) + F_e(r); 
    end
end

%% Compute Fext
Fext = zeros(n_dof,1);
Fext(21) = -Me*g;
for e = 1:n_el
    for i = 1:n_dof_el
        I = Td(e, i);
        Fext(I) = Fext(I) + F_el(i,e);
    end
end

%% Global system of equations
u_R = [0; 0];
v_R = [1; 2];
v_L = (3:n_dof)';

K_LL = KG(v_L, v_L);
K_LR = KG(v_L, v_R);
K_RL = KG(v_R, v_L);
K_RR = KG(v_R, v_R);

F_L = KG(v_L);
F_R = KG(v_R);

u_L = K_LL \ (F_L - K_LR*u_R);
R_R = K_RR*u_R + K_RL*u_L - F_R;

u(v_L,1) = u_L;
u(v_R,1) = u_R;


%% Compute internal distributions

for e = 1:length(n_el)
    x_e = zeros(2, n_d);

    for i = 1:2
        for j = 1:n_d
            x_e(i,j) = X(Tn(e,i),j);
            if i==2
                x_e(i,j) = x_e(i,j);
            end
        end
    end
    
    l = sqrt(sum(sum(x_e).^2));
    u_e = zeros(n_dof_el,1);
    
    for i = 1:n_dof_el
        I = Td(e,i);
        u_e(i,1) = u(I);
    end
    
    F_int_e = Kel(:,:,e)*u_e;
    
    Fy(e, 1) = -F_int_e(1);
    Fy(e, 2) = F_int_e(3);
    Mz(e, 1) = -F_int_e(2);
    Mz(e, 2) = F_int_e(4);
    
    
    A = [2 l -2 l; -3*l -2*l^2 3*l -l^2; 0 l^3 0 0; l^3 0 0 0];
    coef = 1/l^3 * A * u_e;
    pu(e,:) = coef;
    pt(e,:) = [3*coef(1), 2*coef(2), coef(3)];
    
end
for k = 3:3:30
%% POSTPROCESS

% Number of subdivisions and plots
nsub = k;
plotBeams1D(fig,X,Tn,nsub,pu,pt,Fy,Mz)
drawnow;

end

%% Loop through each of the number of elements
for k = 1:length(n_el)

    %% PREPROCESS
    
    % Nodal coordinates
    %  x(a,j) = coordinate of node a in the dimension j
    % Complete the coordinates
    
    % Nodal connectivities  
    %  Tnod(e,a) = global nodal number associated to node a of element e
    
    % Material properties matrix
    %  mat(m,1) = Young modulus of material m
    %  mat(m,2) = Section area of material m
    %  mat(m,3) = Section inertia of material m

    % Material connectivities
    %  Tmat(e) = Row in mat corresponding to the material associated to element e 
    
    
    %% SOLVER
    
    % Compute:
    % u  - Displacements and rotations vector [ndof x 1]
    % pu - Polynomial coefficients for displacements for each element [nel x 4]
    % pt - Polynomial coefficients for rotations for each element [nel x 3]
    % Fy - Internal shear force at each elements's nodes [nel x nne]
    % Mz - Internal bending moment at each elements's nodes [nel x nne]
    
    %% POSTPROCESS
    
    % Number of subdivisions and plots
    nsub = n_el/3;
    plotBeams1D(fig,x,Tnod,nsub,pu,pt,Fy,Mz)
    drawnow;
    
end

% Add figure legends
figure(fig)
legend(strcat('N=',cellstr(string(nel))),'location','northeast');