function diff_stencil = set_diff_st(x,a1,b1,c1,M_old, alpha, gamma2 ,DT_scale,C1)
%Using this function, all the 9 stencils(as mentioned in Table 2.3 of the thesis)
%have been saved as alpha's. alpha1 is the upper-left stencil, alpha2:
%upper-middle and so on.

%% Memory allocation
Lx = length(x);
h = x(2)-x(1);
diff_stencil.alpha1 = zeros(Lx-2,Lx-2);
diff_stencil.alpha2 = zeros(Lx-2,Lx-2);
diff_stencil.alpha3 = zeros(Lx-2,Lx-2);
diff_stencil.alpha4 = zeros(Lx-2,Lx-2);
diff_stencil.alpha5 = zeros(Lx-2,Lx-2);
diff_stencil.alpha6 = zeros(Lx-2,Lx-2);
diff_stencil.alpha7 = zeros(Lx-2,Lx-2);
diff_stencil.alpha8 = zeros(Lx-2,Lx-2);
diff_stencil.alpha9 = zeros(Lx-2,Lx-2);
temp = alpha*M_old;

%% Loop to compute the non-linear diffusion( for explicit method)
for jj = 2:Lx-1
    for ii = 2:Lx-1
        temp(ii,jj) = (gamma2*alpha*M_old(ii,jj))/(sqrt(C1+ (((M_old(ii+1,jj)-M_old(ii-1,jj))/(2*h))^2) ...
                     + (((M_old(ii,jj+1)-M_old(ii,jj-1))/(2*h))^2)));
      
    end
end
%% diffusion tensor update
a = a1 + temp.*(DT_scale-a1);
b = b1 + temp.*(0-b1);
c = c1 + temp.*(DT_scale-c1);

%% Comupation of 9 stencils
for jj = 2:Lx-1
    for ii = 2:Lx-1
        iii = ii-1;
        jjj = jj-1;
        
        diff_stencil.alpha1(iii,jjj) = ((abs(b(ii-1,jj+1)) - b(ii-1,jj+1)) + (abs(b(ii,jj)) - b(ii,jj)))/(4*h*h);
        
        diff_stencil.alpha2(iii,jjj) = ((c(ii,jj+1)+c(ii,jj))/(2*h*h)) - ((abs(b(ii,jj+1))+abs(b(ii,jj)))/(2*h*h));
        
        diff_stencil.alpha3(iii,jjj) = ((abs(b(ii+1,jj+1)) +  b(ii+1,jj+1))/(4*h*h)) + (abs(b(ii,jj))+b(ii,jj))/(4*h*h);
        
        diff_stencil.alpha4(iii,jjj) = ((a(ii-1,jj)+a(ii,jj))/(2*h*h)) - ((abs(b(ii-1,jj))+abs(b(ii,jj)))/(2*h*h));
        
        diff_stencil.alpha5(iii,jjj) =  -((a(ii-1,jj)+2*a(ii,jj)+a(ii+1,jj))/(2*h*h)) - ((abs(b(ii-1,jj+1))-b(ii-1,jj+1)+...
            abs(b(ii+1,jj+1))+b(ii+1,jj+1))/(4*h*h)) - ((abs(b(ii-1,jj-1))+b(ii-1,jj-1)+...
            abs(b(ii+1,jj-1))-b(ii+1,jj-1))/(4*h*h)) + ((abs(b(ii-1,jj))+abs(b(ii+1,jj))+...
            abs(b(ii,jj-1))+abs(b(ii,jj+1))+2*abs(b(ii,jj)))/(2*h*h)) - ((c(ii,jj-1)+2*c(ii,jj)+...
            c(ii,jj+1))/(2*h*h));
        
        diff_stencil.alpha6(iii,jjj) = ((a(ii+1,jj)+a(ii,jj))/(2*h*h)) - ((abs(b(ii+1,jj))+abs(b(ii,jj)))/(2*h*h));
        
        diff_stencil.alpha7(iii,jjj) = ((abs(b(ii-1,jj-1))+b(ii-1,jj-1))/(4*h*h)) + (abs(b(ii,jj))+b(ii,jj))/(4*h*h);
        
        diff_stencil.alpha8(iii,jjj) = ((c(ii,jj-1)+c(ii,jj))/(2*h*h)) - ((abs(b(ii,jj-1))+abs(b(ii,jj)))/(2*h*h));
        
        diff_stencil.alpha9(iii,jjj) = ((abs(b(ii+1,jj-1))-b(ii+1,jj-1))/(4*h*h)) + ((abs(b(ii,jj)) - b(ii,jj))/(4*h*h));
       
        
    
    end
end

end