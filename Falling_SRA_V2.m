% 2022.12.20

clear all; clc;
format short;

% Fixed position in inertia frame
h=1; % The height for the base of the arm is fixed.
% For the x direction position the SRA is fixed, it should be relative to
% the distance between the two walls d.

% The number of link for each side of the SRA.
n_link=2;

L=0.6; % Length of each SRA
m=2.5; % Mass of each SRA
d=1; % Distance between two walls
g=9.81; % Gravity acceleration

k=2; % Stiffness of spring between each two adjacent COMs of links.
c=1; % Damping coefficient of the damper between the COMs of two adjacent links.


L_theta=sym('L_theta',[1,n_link]);
R_theta=sym('R_theta',[1,n_link]);
L_theta_dot=sym('L_theta_dot',[1,n_link]);
R_theta_dot=sym('R_theta_dot',[1,n_link]);
% Temporarily the rotation matrix for left and right SRA are the same, becasue the ball is released from the
%  center of the wall and the two SRAs are equal length, so the motion and animation will be symmetric 
% about the central line of the wall.
R=struct('L',zeros(4,4),'R',zeros(4,4)); % The rotation matrix for left and right SRA are different.
for i=1:n_link
    % The expression of the standand derived transformation/rotation
    % matrix.
    R(i).L=[cos(-L_theta(i)) -sin(-L_theta(i)) 0 0;
        sin(-L_theta(i)) cos(-L_theta(i)) 0 0;
        0 0 1 0;
        0 0 0 1]*[1 0 0 0;
        0 1 0 L/2;
        0 0 1 0;
        0 0 0 1];
    R(i).R=[cos(-R_theta(i)) -sin(-R_theta(i)) 0 0;
        sin(-R_theta(i)) cos(-R_theta(i)) 0 0;
        0 0 1 0;
        0 0 0 1]*[1 0 0 0;
        0 1 0 L/2;
        0 0 1 0;
        0 0 0 1];
end

% Calculate the intermediate rotation matrix.
R_intm=struct('L',zeros(4,4),'R',zeros(4,4));
for j=1:n_link
    if j==1
        R_intm(j).L=R(j).L;
        R_intm(j).R=R(j).R;
    else
        R_intm(j).L=simplify(R_intm(j-1).L*R(j).L);
        R_intm(j).R=simplify(R_intm(j-1).R*R(j).R);
    end
end

% The next step is to work out  the real position of each joint in the
% inertia fram by these rotation matrices calculated above.

% The rotation matrix from F0 to F1
% Assume the origin is located at the center of two walls on the ground.
R_01=struct('L',zeros(4,4),'R',zeros(4,4));
R_01.L=[1 0 0 -d/2;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]*[1 0 0 0;
    0 1 0 h;
    0 0 1 0;
    0 0 1 0]*[0 1 0 0;
    -1 0 0 0;
    0 0 1 0;
    0 0 0 1];
R_01.R=[1 0 0 d/2;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]*[1 0 0 0;
    0 1 0 h;
    0 0 1 0;
    0 0 0 1]*[0 -1 0 0;
    1 0 0 0;
    0 0 1 0;
    0 0 0 1]*[-1 0 0 0;
    0 1 0 0;
    0 0 -1 0;
    0 0 0 1];

R_abs=struct('L',zeros(4,4),'R',zeros(4,4));
for i=1:n_link
    R_abs(i).L=simplify(R_01.L*R_intm(i).L);
    R_abs(i).R=simplify(R_01.R*R_intm(i).R);
end

o_c=struct('L',zeros(3,1),'R',zeros(3,1));
z=struct('L',zeros(3,1),'R',zeros(3,1));
o0=[0;0;0];
for i=1:n_link
    if i==1
        o_c(i).L=1/2*R_abs(i).L(1:3,4);
        o_c(i).R=1/2*R_abs(i).R(1:3,4);
    else
        o_c(i).L=R_abs(i).L(1:3,4)+1/2*(R_abs(i).L(1:3,4)-R_abs(i-1).L(1:3,4));
        o_c(i).R=R_abs(i).R(1:3,4)+1/2*(R_abs(i).R(1:3,4)-R_abs(i-1).R(1:3,4));
    end
end
z0=[0;0;1];
for i=1:n_link
    z(i).L=R_abs(i).L(1:3,3);
    z(i).R=R_abs(i).R(1:3,3);
end

% Jacobian matrix
J_cell=struct('L',zeros(6,n_link),'R',zeros(6,n_link));
% Due to the limitations of struct array, we could not define the Jacobian
% matrix with the format of struct, since the code is not written in
% Simulink, we can try cell.
J_intm_L=cell(2,n_link);
J_intm_R=cell(2,n_link);
for i=1:n_link
    for j=1:2
        J_intm_L{j,i}=[0;0;0];
        J_intm_R{j,i}=[0;0;0];
    end
end
for i=1:n_link
    if i==1
            % Define all the intermediate variables again to avoid the
            % issue that cross function cannot deal with symbolic
            % variables.
%             z_sym=sym(z0);
%             oc_L1=sym(o_c(i).L-o0);
%             oc_R1=sym(o_c(i).R-o0);
%             J(i).L(1:3,j)=cross(z_sym,oc_L1);
%             J(i).R(1:3,j)=cross(z_sym,oc_R1);

            J_intm_L{1,i}=cross(z0,(o_c(i).L-o0));
            J_intm_R{1,i}=cross(z0,(o_c(i).R-o0));
            J_intm_L{2,i}=z0;
            J_intm_R{2,i}=z0;

%                  J(i).L(1:3,j)=cross(z0,(o_c(i).L-o0));
%                  J(i).R(1:3,j)=cross(z0,(o_c(i).R-o0));
    else
        for j=1:n_link
            if j==1
                J_intm_L{1,j}=cross(z0,(o_c(i).L-o0));
                J_intm_R{1,j}=cross(z0,(o_c(i).R-o0));
                J_intm_L{2,j}=z0;
                J_intm_R{2,j}=z0;
            else
%             z_intm_L2=sym(z(j-1).L);
%             z_intm_R2=sym(z(j-1).R);
%             oc_L2=sym((o_c(i).L-R_abs(i-1).L(1:3,4)));
%             oc_R2=sym((o_c(i).R-R_abs(i-1).R(1:3,4)));
%             J(i).L(1:3,j)=cross(z_intm_L2,oc_L2);
%             J(i).R(1:3,j)=cross(z_intm_R2,oc_R2);

            J_intm_L{1,j}=cross(z(j-1).L,(o_c(i).L-R_abs(j-1).L(1:3,4)));
            J_intm_R{1,j}=cross(z(j-1).R,(o_c(i).L-R_abs(j-1).L(1:3,4)));
            J_intm_L{2,j}=z(j-1).L;
            J_intm_R{2,j}=z(j-1).R;

%             J(i).L(1:3,j)=cross(z(j-1).L,(o_c(i).L-R_abs(j-1).L(1:3,4)));
%             J(i).R(1:3,j)=cross(z(j-1).R,(o_c(i).L-R_abs(j-1).L(1:3,4)));
            end
        end
    end
    J_cell(i).L=J_intm_L;
    J_cell(i).R=J_intm_R;
end
J=struct('L',zeros(6,n_link),'R',zeros(6,n_link));
for i=1:n_link
    J(i).L=cell2sym(J_cell(i).L);
    J(i).R=cell2sym(J_cell(i).R);
end

% Moment of inertia
I=1/12*(m/2)*(L/2)^2;

% Inertia matrix (M matrix)
M=struct('L',zeros(n_link,n_link),'R',zeros(n_link,n_link));
for i=1:n_link
    M.L=simplify(M.L+m/2*transpose(J(i).L(1:3,:))*J(i).L(1:3,:)+transpose(J(i).L(4:6,:))...
        *R_abs(i).L(1:3,1:3)*I*transpose(R_abs(i).L(1:3,1:3))*J(i).L(1:3,:));
    M.R=simplify(M.R+m/2*transpose(J(i).R(1:3,:))*J(i).R(1:3,:)+transpose(J(i).R(4:6,:))...
        *R_abs(i).R(1:3,1:3)*I*transpose(R_abs(i).R(1:3,1:3))*J(i).R(1:3,:));
end

% Coriolis-centripedal matrix (C matrix)
C=struct('L',zeros(n_link,n_link),'R',zeros(n_link,n_link));
C_element_L=cell(n_link,n_link);
C_element_R=cell(n_link,n_link);
Element_L=0; Element_R=0;
for i=1:n_link
    for j=1:n_link
        for k=1:n_link
%             C.L=simplify(C.L+1/2*(diff(M.L(i,j),theta(k))+diff(M.L(i,k),theta(j))...
%                 -diff(M.L(k,j),theta(i)))*theta_dot(k));
%             C.R=simplify(C.R+1/2*(diff(M.R(i,j),theta(k))+diff(M.R(i,k),theta(j))...
%                 -diff(M.R(k,j),theta(i)))*theta_dot(k));
            Element_L=Element_L+1/2*(diff(M.L(i,j),L_theta(k))+diff(M.L(i,k),L_theta(j))...
                -diff(M.L(k,j),L_theta(i)))*L_theta_dot(k);
            Element_R=Element_R+1/2*(diff(M.R(i,j),R_theta(k))+diff(M.R(i,k),R_theta(j))...
                -diff(M.R(k,j),R_theta(i)))*R_theta_dot(k);
        end
        C_element_L{i,j}=Element_L;
        C_element_R{i,j}=Element_R;
    end
end
C.L=simplify(cell2sym(C_element_L));
C.R=simplify(cell2sym(C_element_R));

%% Potential Energy Part %%
% Potential energy due to gravity
Pg_total=struct('L',0,'R',0);
Pg_cell_L=cell(n_link,1); Pg_cell_R=cell(n_link,1);
Pg=struct('L',zeros(n_link,1),'R',zeros(n_link,1));
for i=1:n_link
    Pg_total.L=Pg_total.L+m/2*g*o_c(i).L(2);
    Pg_total.R=Pg_total.R+m/2*g*o_c(i).R(2);
end
for i=1:n_link
    Pg_cell_L{i,1}=(diff(Pg_total.L,L_theta(i)));
    Pg_cell_R{i,1}=(diff(Pg_total.R,R_theta(i)));
end
Pg.L=cell2sym(Pg_cell_L);
Pg.R=cell2sym(Pg_cell_R);

% Potential energy due to "spring" stiffness between the COMs of each
% adjacent link.
% We assume that the unstretched length of the "spring" is L, we can adjust
% this according to our meet in the future.
dc=struct('L',0,'R',0); % Distance between COMs.
for i=1:(n_link-1)
    dc(i).L=((o_c(i+1).L(1)-o_c(i).L(1))^2+(o_c(i+1).L(2)-o_c(i).L(2))^2)^(1/2);
    dc(i).R=((o_c(i+1).R(1)-o_c(i).R(1))^2+(o_c(i+1).R(2)-o_c(i).R(2))^2)^(1/2);
end
Pk_total=struct('L',0,'R',0);
Pk_cell_L=cell(n_link,1); Pk_cell_R=cell(n_link,1);
Pk=struct('L',zeros(n_link,1),'R',zeros(n_link,1));
for i=1:(n_link-1)
    Pk_total.L=Pk_total.L+1/2*k*(L-dc(i).L)^2;
    Pk_total.R=Pk_total.R+1/2*k*(L-dc(i).R)^2;
end
for i=1:n_link
    Pk_cell_L{i,1}=simplify(diff(Pk_total.L,L_theta(i)));
    Pk_cell_R{i,1}=simplify(diff(Pk_total.R,R_theta(i)));
end
Pk.L=cell2sym(Pk_cell_L);
Pk.R=cell2sym(Pk_cell_R);

% Potential energy due to damping coefficient
vc_rel=struct('L',zeros(3,1),'R',zeros(3,1)); % The relative velocities between two adjacent COMs.
for i=1:(n_link-1)
    vc_rel(i).L=((jacobian(o_c(i+1).L,L_theta)-jacobian(o_c(i).L,L_theta))...
    *transpose(L_theta_dot));
    vc_rel(i).R=((jacobian(o_c(i+1).R,R_theta)-jacobian(o_c(i).R,R_theta))...
    *transpose(R_theta_dot));
end
Pc_total=struct('L',0,'R',0);
Pc_cell_L=cell(n_link,1); Pc_cell_R=cell(n_link,1);
Pc=struct('L',zeros(n_link,1),'R',zeros(n_link,1));
for i=1:(n_link-1)
    Pc_total.L=Pc_total.L+1/2*c*transpose(vc_rel(i).L)*vc_rel(i).L;
    Pc_total.R=Pc_total.R+1/2*c*transpose(vc_rel(i).R)*vc_rel(i).R;
end
for i=1:n_link
    Pc_cell_L{i,1}=diff(Pc_total.L,L_theta_dot(i));
    Pc_cell_R{i,1}=diff(Pc_total.L,L_theta_dot(i));
end
Pc.L=simplify(cell2sym(Pc_cell_L));
Pc.R=simplify(cell2sym(Pc_cell_R));

%% External Force Calculation %%



