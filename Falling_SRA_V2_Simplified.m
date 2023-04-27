% 2023.01.03
% This modified version is for the special case, when the ball is released
% from an arbitrary point along the central line between the two walls, and
% in this case the motions of the two arms will be exactly the same or
% symmetric about the central line.

% ABOUT EXTERNAL FORCE EXERTED ON THE SRA TIP
% We can treat the two arms and the ball as a whole object, thus the
% interaction force (or external force exerted on the tip) can be ignored.

clear all; clc;
format short;

% Fixed position in inertia frame
par.h=1; % The height for the base of the arm is fixed.
% For the x direction position the SRA is fixed, it should be relative to
% the distance between the two walls d.

% The number of link for each side of the SRA.
par.n_link=2;

par.L=0.6; % Length of each SRA
par.m=2.5; % Mass of each SRA
par.mb=1; % Mass of the ball
par.d=1; % Distance between two walls
par.g=9.81; % Gravity acceleration

par.k=2; % Stiffness of spring between each two adjacent COMs of links.
par.c=1; % Damping coefficient of the damper between the COMs of two adjacent links.

% [M,C,Pg,Pk,Pc,ddtheta]=Matrix_Derivation(par);
[M,ddtheta]=Matrix_Derivation(par);
eom_intm=eom_proc(ddtheta,par);
save('v2s.mat');

%% ode function setting up using implicit method.
% theta=sym('theta',[1,par.n_link]);
% theta_dot=sym('theta_dot',[1,par.n_link]);
% eom=@(t,q) [q(3);
%     q(4);
%     subs(ddtheta(1),[theta,theta_dot],[q(1),q(2),q(3),q(4)]);
%     subs(ddtheta(2),[theta,theta_dot],[q(1),q(2),q(3),q(4)])];

% eom=@(t,theta) [theta_dot(1);
%     theta_  dot(2);
%     ddtheta(1);
%     ddtheta(2)];

% [t,y]=ode23(eom,tspan,y0);


% 2023.01.09
% The ode function set up has been failed, probably we can try another
% method above.
% function dy=eom(t,y,ddtheta)
% theta=sym('theta',[1,par.n_link]);
% theta_dot=sym('theta_dot',[1,par.n_link]);
% dy(1)=y(2);
% dy(2)=subs(ddtheta,[theta,theta_dot],[y(1),y(2)]);    
% end

% 2023.01.10
% Try another method to transfer the symbolic variable expression into the
% format that ode function could accept.
function post_process=eom_proc(ddtheta,par)
post_process=cell(par.n_link,1);
theta=sym('theta',[1,par.n_link]);
theta_dot=sym('theta_dot',[1,par.n_link]);
intm_1ord=cell(par.n_link,1); intm_2ord=cell(par.n_link,1);
for i=1:par.n_link
    intm_1ord{i,1}=char(theta(i));
    intm_2ord{i,1}=char(theta_dot(i));
end
for i=1:par.n_link
    post_process{i,1}=char(ddtheta(i));
end
for i=1:par.n_link
    post_process=replace(post_process,intm_1ord{i,1},['q(',num2str(i),')']);
    post_process=replace(post_process,intm_2ord{i,1},['q(',num2str(i+2),')']);
end
post_process=replace(post_process,'*','.*');
post_process=replace(post_process,'^','.^');
post_process=str2sym(post_process);
end



% function [M_coupled,C_coupled,Pg_coupled,Pk_coupled,Pc_coupled,ddtheta]=Matrix_Derivation(par)
function [M_coupled,ddtheta]=Matrix_Derivation(par)

L_theta=sym('L_theta',[1,par.n_link]);
R_theta=sym('R_theta',[1,par.n_link]);
L_theta_dot=sym('L_theta_dot',[1,par.n_link]);
R_theta_dot=sym('R_theta_dot',[1,par.n_link]);
% Temporarily the rotation matrix for left and right SRA are the same, becasue the ball is released from the
%  center of the wall and the two SRAs are equal length, so the motion and animation will be symmetric 
% about the central line of the wall.
R=struct('L',zeros(4,4),'R',zeros(4,4)); % The rotation matrix for left and right SRA are different.
for i=1:par.n_link
    % The expression of the standand derived transformation/rotation
    % matrix.
    R(i).L=[cos(-L_theta(i)) -sin(-L_theta(i)) 0 0;
        sin(-L_theta(i)) cos(-L_theta(i)) 0 0;
        0 0 1 0;
        0 0 0 1]*[1 0 0 0;
        0 1 0 par.L/2;
        0 0 1 0;
        0 0 0 1];
    R(i).R=[cos(-R_theta(i)) -sin(-R_theta(i)) 0 0;
        sin(-R_theta(i)) cos(-R_theta(i)) 0 0;
        0 0 1 0;
        0 0 0 1]*[1 0 0 0;
        0 1 0 par.L/2;
        0 0 1 0;
        0 0 0 1];
end

% Calculate the intermediate rotation matrix.
R_intm=struct('L',zeros(4,4),'R',zeros(4,4));
for j=1:par.n_link
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
R_01.L=[1 0 0 -par.d/2;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]*[1 0 0 0;
    0 1 0 par.h;
    0 0 1 0;
    0 0 1 0]*[0 1 0 0;
    -1 0 0 0;
    0 0 1 0;
    0 0 0 1];
R_01.R=[1 0 0 par.d/2;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1]*[1 0 0 0;
    0 1 0 par.h;
    0 0 1 0;
    0 0 0 1]*[0 -1 0 0;
    1 0 0 0;
    0 0 1 0;
    0 0 0 1]*[-1 0 0 0;
    0 1 0 0;
    0 0 -1 0;
    0 0 0 1];

R_abs=struct('L',zeros(4,4),'R',zeros(4,4));
for i=1:par.n_link
    R_abs(i).L=simplify(R_01.L*R_intm(i).L);
    R_abs(i).R=simplify(R_01.R*R_intm(i).R);
end

o_c=struct('L',zeros(3,1),'R',zeros(3,1));
z=struct('L',zeros(3,1),'R',zeros(3,1));
o0=[0;0;0];
for i=1:par.n_link
    if i==1
        o_c(i).L=1/2*R_abs(i).L(1:3,4);
        o_c(i).R=1/2*R_abs(i).R(1:3,4);
    else
        o_c(i).L=R_abs(i).L(1:3,4)+1/2*(R_abs(i).L(1:3,4)-R_abs(i-1).L(1:3,4));
        o_c(i).R=R_abs(i).R(1:3,4)+1/2*(R_abs(i).R(1:3,4)-R_abs(i-1).R(1:3,4));
    end
end
z0=[0;0;1];
for i=1:par.n_link
    z(i).L=R_abs(i).L(1:3,3);
    z(i).R=R_abs(i).R(1:3,3);
end

% Jacobian matrix
J_cell=struct('L',zeros(6,par.n_link),'R',zeros(6,par.n_link));
% Due to the limitations of struct array, we could not define the Jacobian
% matrix with the format of struct, since the code is not written in
% Simulink, we can try cell.
J_intm_L=cell(2,par.n_link);
J_intm_R=cell(2,par.n_link);
for i=1:par.n_link
    for j=1:2
        J_intm_L{j,i}=[0;0;0];
        J_intm_R{j,i}=[0;0;0];
    end
end
for i=1:par.n_link
    if i==1
            % Define all the intermediate pariables again to avoid the
            % issue that cross function cannot deal with symbolic
            % pariables.
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
        for j=1:par.n_link
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
            J_intm_R{1,j}=cross(z(j-1).R,(o_c(i).R-R_abs(j-1).R(1:3,4)));
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
J=struct('L',zeros(6,par.n_link),'R',zeros(6,par.n_link));
for i=1:par.n_link
    J(i).L=cell2sym(J_cell(i).L);
    J(i).R=cell2sym(J_cell(i).R);
end

% Moment of inertia
I=1/12*(par.m/2)*(par.L/2)^2;

% Inertia matrix (M matrix)
M=struct('L',zeros(par.n_link,par.n_link),'R',zeros(par.n_link,par.n_link));
for i=1:par.n_link
    M.L=simplify(M.L+par.m/2*transpose(J(i).L(1:3,:))*J(i).L(1:3,:)+transpose(J(i).L(4:6,:))...
        *R_abs(i).L(1:3,1:3)*I*transpose(R_abs(i).L(1:3,1:3))*J(i).L(1:3,:));
    M.R=simplify(M.R+par.m/2*transpose(J(i).R(1:3,:))*J(i).R(1:3,:)+transpose(J(i).R(4:6,:))...
        *R_abs(i).R(1:3,1:3)*I*transpose(R_abs(i).R(1:3,1:3))*J(i).R(1:3,:));
end

% Coriolis-centripedal matrix (C matrix)
C=struct('L',zeros(par.n_link,par.n_link),'R',zeros(par.n_link,par.n_link));
C_element_L=cell(par.n_link,par.n_link);
C_element_R=cell(par.n_link,par.n_link);
Element_L=0; Element_R=0;
for i=1:par.n_link
    for j=1:par.n_link
        for k=1:par.n_link
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
Pg_cell_L=cell(par.n_link,1); Pg_cell_R=cell(par.n_link,1);
Pg=struct('L',zeros(par.n_link,1),'R',zeros(par.n_link,1));
for i=1:par.n_link
    Pg_total.L=Pg_total.L+par.m/2*par.g*o_c(i).L(2);
    Pg_total.R=Pg_total.R+par.m/2*par.g*o_c(i).R(2);
end
for i=1:par.n_link
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
for i=1:(par.n_link-1)
    dc(i).L=((o_c(i+1).L(1)-o_c(i).L(1))^2+(o_c(i+1).L(2)-o_c(i).L(2))^2)^(1/2);
    dc(i).R=((o_c(i+1).R(1)-o_c(i).R(1))^2+(o_c(i+1).R(2)-o_c(i).R(2))^2)^(1/2);
end
Pk_total=struct('L',0,'R',0);
Pk_cell_L=cell(par.n_link,1); Pk_cell_R=cell(par.n_link,1);
Pk=struct('L',zeros(par.n_link,1),'R',zeros(par.n_link,1));
for i=1:(par.n_link-1)
    Pk_total.L=Pk_total.L+1/2*par.k*(par.L-dc(i).L)^2;
    Pk_total.R=Pk_total.R+1/2*par.k*(par.L-dc(i).R)^2;
end
for i=1:par.n_link
    Pk_cell_L{i,1}=simplify(diff(Pk_total.L,L_theta(i)));
    Pk_cell_R{i,1}=simplify(diff(Pk_total.R,R_theta(i)));
end
Pk.L=cell2sym(Pk_cell_L);
Pk.R=cell2sym(Pk_cell_R);

% Potential energy due to damping coefficient
vc_rel=struct('L',zeros(3,1),'R',zeros(3,1)); % The relative velocities between two adjacent COMs.
for i=1:(par.n_link-1)
    vc_rel(i).L=((jacobian(o_c(i+1).L,L_theta)-jacobian(o_c(i).L,L_theta))...
    *transpose(L_theta_dot));
    vc_rel(i).R=((jacobian(o_c(i+1).R,R_theta)-jacobian(o_c(i).R,R_theta))...
    *transpose(R_theta_dot));
end
Pc_total=struct('L',0,'R',0);
Pc_cell_L=cell(par.n_link,1); Pc_cell_R=cell(par.n_link,1);
Pc=struct('L',zeros(par.n_link,1),'R',zeros(par.n_link,1));
for i=1:(par.n_link-1)
    Pc_total.L=Pc_total.L+1/2*par.c*transpose(vc_rel(i).L)*vc_rel(i).L;
    Pc_total.R=Pc_total.R+1/2*par.c*transpose(vc_rel(i).R)*vc_rel(i).R;
end
for i=1:par.n_link
    Pc_cell_L{i,1}=diff(Pc_total.L,L_theta_dot(i));
    Pc_cell_R{i,1}=diff(Pc_total.R,R_theta_dot(i));
end
Pc.L=simplify(cell2sym(Pc_cell_L));
Pc.R=simplify(cell2sym(Pc_cell_R));

%%   Special treatment for symmetric motions of left and right arms   %%
theta=sym('theta',[1,par.n_link]);
theta_dot=sym('theta_dot',[1,par.n_link]);
M.L=subs(M.L,L_theta,theta);
M.R=subs(M.R,R_theta,theta);
M_coupled=M.L+M.R;
C.L=subs(C.L,[L_theta,L_theta_dot],[theta,theta_dot]);
C.R=subs(C.R,[R_theta,R_theta_dot],[theta,theta_dot]);
C_coupled=C.L+C.R;
Pk.L=subs(Pk.L,L_theta,theta); Pk.R=subs(Pk.R,R_theta,theta);
Pk_coupled=Pk.L+Pk.R;
Pc.L=subs(Pc.L,[L_theta,L_theta_dot],[theta,theta_dot]);
Pc.R=subs(Pc.R,[R_theta,R_theta_dot],[theta,theta_dot]);
Pc_coupled=Pc.L+Pc.R;

% For treating the potential energy due to gravity, it's kind of different,
% since there is another term due to the gravity of the ball, we have to
% calculate the total energy again.
Pg_total.L=subs(Pg_total.L,L_theta,theta); Pg_total.R=subs(Pg_total.R,R_theta,theta);
Pg_total=Pg_total.L+Pg_total.R+par.mb*par.g*R_abs(end).L(2,4);
Pg_cell=cell(par.n_link,1);
for i=1:par.n_link
    Pg_cell{i,1}=(diff(Pg_total,theta(i)));
end
Pg_coupled=cell2sym(Pg_cell);
%% Processing End %%

RHS=(-C_coupled*transpose(theta_dot)-Pg_coupled-Pk_coupled-Pc_coupled);
ddtheta=simplify(M_coupled\RHS);

end


