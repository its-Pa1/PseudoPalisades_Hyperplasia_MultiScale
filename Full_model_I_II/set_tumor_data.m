function [q,Q_value,a,b,c,FA,DivDT_x,DivDT_y] = ...
    set_tumor_data(x,y, delta,kappa, center_x_DW,center_y_DW, DT_scale)
% This function deals with everything required for tumor diffusion
% computation

% outputs:
% q = mesoscopic tissue(q)
% Q_value = macroscopic tissue(Q)
% a = DT(1,1) first element of first row of tumor diff tensor
% b = DT(1,2) = DT(2,1) 
% c = DT(2,2)
% FA = fractional anisotropy
% DivDT_x = x component of div of DT
% DivDT_y = y component of div of DT

% inputs:
% x,y : the space variables
% dt = time spacing
% delta = required for computing small q, it combines the two distributions
% kappa = it controls the FA of D_w
% center_x_DW = x position of the cross for D_w(as our D_w makes a cross type sign)
% center_y_DW = y position of the cross for D_w
% DT_scale = s^2/lambda0

%% memory allocations
Lx = length(x);
h = x(2)-x(1);
theta = linspace(0,pi,101);
q = zeros(size(x,2),size(y,2),size(theta,2));
Q_value = zeros(size(x,2),size(y,2));
DC = zeros(2,2,size(x,2),size(y,2));
FA = zeros(size(x,2), size(y,2));

a = zeros(size(x,2),size(y,2));
b = zeros(size(x,2),size(y,2));
c = zeros(size(x,2),size(y,2));

DivDT_x = zeros(Lx-2,Lx-2);
DivDT_y = zeros(Lx-2,Lx-2);

%% Loops to save the Q and q at all points and theta(small q depends of theta)
for jj = 1:length(y)
    for ii = 1:length(x)
        Q_value(ii,jj) = tissue_Q_macro(x(ii),y(jj),center_x_DW,center_y_DW);
        for l = 1:length(theta)
            q(ii,jj,l) = tissue_q_un(x(ii),y(jj),theta(l),delta, kappa, center_x_DW, center_y_DW);
        end
        
    end
end
%% Computation of DT and storing its components in a, b and c
% here trapezoidal method is used for numerical integration
for jj = 1:length(y)
    for ii = 1:length(x)
        sum2 = 0;
        
        for l = 2:length(theta)-1
            sum2 = sum2 + q(ii,jj,l)*([cos(theta(l));sin(theta(l))]*[cos(theta(l));sin(theta(l))]');
        end
        DC(:,:,ii,jj) = DT_scale*2*0.5*(theta(2)-theta(1))*(2*sum2 ...
            + q(ii,jj,1)*[cos(theta(1));sin(theta(1))]*[cos(theta(1));sin(theta(1))]'...
            + q(ii,jj,end)*[cos(theta(end));sin(theta(end))]*[cos(theta(end));sin(theta(end))]');
        
        temp_d = DC(:,:,ii,jj);
        a(ii,jj) = temp_d(1,1);
        b(ii,jj) = temp_d(1,2);
        c(ii,jj) = temp_d(2,2);
    end
end
%% Here FA is stored at all points in the domain
for jj = 1:length(y)
    for ii = 1:length(x)
        lambda = eig(DC(:,:,ii,jj));
        FA(ii,jj) = (abs(lambda(1)-lambda(2)))/(sqrt((lambda(1)^2)+(lambda(2)^2)));
        
    end
end
%% in these loops, the div of DT has been stored
for jj = 2:Lx-1
    for ii = 2:Lx-1
        iii = ii-1;
        jjj = jj-1;
        
        DivDT_x(iii,jjj) = (a(ii,jj) - a(ii-1,jj) + b(ii,jj) - b(ii-1,jj))/h;
        DivDT_y(iii,jjj) = (b(ii,jj) - b(ii,jj-1) + c(ii,jj) - c(ii,jj-1))/h;
    end
end

end