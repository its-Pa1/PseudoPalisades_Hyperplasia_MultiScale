function [q,Q_value,a,b,c,diff_stencil,FA,DivDT_x,DivDT_y] = ...
    set_tumor_diff_explicit(x,y, delta,kappa, center_x_DW,center_y_DW, DT_scale)
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
% theta = linspace(0,pi,100);
theta = linspace(0,pi,101);
q = zeros(size(x,2),size(y,2),size(theta,2));
Q_value = zeros(size(x,2),size(y,2));
DC = zeros(2,2,size(x,2),size(y,2));
FA = zeros(size(x,2), size(y,2));
diff_stencil.alpha1 = zeros(Lx-2,Lx-2);
diff_stencil.alpha2 = zeros(Lx-2,Lx-2);
diff_stencil.alpha3 = zeros(Lx-2,Lx-2);
diff_stencil.alpha4 = zeros(Lx-2,Lx-2);
diff_stencil.alpha5 = zeros(Lx-2,Lx-2);
diff_stencil.alpha6 = zeros(Lx-2,Lx-2);
diff_stencil.alpha7 = zeros(Lx-2,Lx-2);
diff_stencil.alpha8 = zeros(Lx-2,Lx-2);
diff_stencil.alpha9 = zeros(Lx-2,Lx-2);
a = zeros(size(x,2),size(y,2));
b = zeros(size(x,2),size(y,2));
c = zeros(size(x,2),size(y,2));
DivDT_x = zeros(Lx-2,Lx-2);
DivDT_y = zeros(Lx-2,Lx-2);

%% Loops to save the Q and q at all points and theta(small q depends of theta)
for j = 1:length(y)
    for i = 1:length(x)
        Q_value(i,j) = tissue_Q_macro(x(i),y(j),center_x_DW,center_y_DW);
        for l = 1:length(theta)
            q(i,j,l) = tissue_q_un(x(i),y(j),theta(l),delta, kappa, center_x_DW, center_y_DW);
        end
        
    end
end
%% Computation of DT and storing its components in a, b and c
% here trapezoidal method is used for numerical integration
for j = 1:length(y)
    for i = 1:length(x)
        sum2 = 0;
        
        for l = 2:length(theta)-1
            sum2 = sum2 + q(i,j,l)*([cos(theta(l));sin(theta(l))]*[cos(theta(l));sin(theta(l))]');
        end
        DC(:,:,i,j) = DT_scale*2*0.5*(theta(2)-theta(1))*(2*sum2 ...
            + q(i,j,1)*[cos(theta(1));sin(theta(1))]*[cos(theta(1));sin(theta(1))]'...
            + q(i,j,end)*[cos(theta(end));sin(theta(end))]*[cos(theta(end));sin(theta(end))]');
        
        temp_d = DC(:,:,i,j);
        a(i,j) = temp_d(1,1);
        b(i,j) = temp_d(1,2);
        c(i,j) = temp_d(2,2);
    end
end
%% Here FA is stored at all points in the domain
for j = 1:length(y)
    for i = 1:length(x)
        lambda = eig(DC(:,:,i,j));
        FA(i,j) = (abs(lambda(1)-lambda(2)))/(sqrt((lambda(1)^2)+(lambda(2)^2)));
    end
end
%% In this section, all the 9 stencils(as mentioned in Weikart et al, pdf in this folder)
%have been saved as alpha's. alpha1 is the upper-left stencil, alpha2:
%upper-middle and so on. Also, in these loops, the div of DT has been
%stored

for j = 2:Lx-1
    for i = 2:Lx-1
        ii = i-1;
        jj = j-1;
        
        diff_stencil.alpha1(ii,jj) = ((abs(b(i-1,j+1)) - b(i-1,j+1)) + (abs(b(i,j)) - b(i,j)))/(4*h*h);
        
        diff_stencil.alpha2(ii,jj) = ((c(i,j+1)+c(i,j))/(2*h*h)) - ((abs(b(i,j+1))+abs(b(i,j)))/(2*h*h));
        
        diff_stencil.alpha3(ii,jj) = ((abs(b(i+1,j+1)) +  b(i+1,j+1))/(4*h*h)) + (abs(b(i,j))+b(i,j))/(4*h*h);
        
        diff_stencil.alpha4(ii,jj) = ((a(i-1,j)+a(i,j))/(2*h*h)) - ((abs(b(i-1,j))+abs(b(i,j)))/(2*h*h));
        
        diff_stencil.alpha5(ii,jj) =  -((a(i-1,j)+2*a(i,j)+a(i+1,j))/(2*h*h)) - ((abs(b(i-1,j+1))-b(i-1,j+1)+...
            abs(b(i+1,j+1))+b(i+1,j+1))/(4*h*h)) - ((abs(b(i-1,j-1))+b(i-1,j-1)+...
            abs(b(i+1,j-1))-b(i+1,j-1))/(4*h*h)) + ((abs(b(i-1,j))+abs(b(i+1,j))+...
            abs(b(i,j-1))+abs(b(i,j+1))+2*abs(b(i,j)))/(2*h*h)) - ((c(i,j-1)+2*c(i,j)+...
            c(i,j+1))/(2*h*h));
        
        diff_stencil.alpha6(ii,jj) = ((a(i+1,j)+a(i,j))/(2*h*h)) - ((abs(b(i+1,j))+abs(b(i,j)))/(2*h*h));
        
        diff_stencil.alpha7(ii,jj) = ((abs(b(i-1,j-1))+b(i-1,j-1))/(4*h*h)) + (abs(b(i,j))+b(i,j))/(4*h*h);
        
        diff_stencil.alpha8(ii,jj) = ((c(i,j-1)+c(i,j))/(2*h*h)) - ((abs(b(i,j-1))+abs(b(i,j)))/(2*h*h));
        
        diff_stencil.alpha9(ii,jj) = ((abs(b(i+1,j-1))-b(i+1,j-1))/(4*h*h)) + ((abs(b(i,j)) - b(i,j))/(4*h*h));
       
        DivDT_x(ii,jj) = (a(i,j) - a(i-1,j) + b(i,j) - b(i-1,j))/h;
        DivDT_y(ii,jj) = (b(i,j) - b(i,j-1) + c(i,j) - c(i,j-1))/h;
    end
end

end