function A = set_up_const_diff(x,dt,Ds)

h = x(2)-x(1);
Lx = length(x);

% acidity diff matrix
A  = diag((1/dt) + (2*Ds/(h*h)) * ones(Lx-2,1), 0) + diag(-Ds/(h*h)* ones(Lx-3,1), -1) ...
    + diag(- (Ds/(h*h))* ones(Lx-3,1), 1);

A(1,1) = A(1,1)-Ds/(h*h); % BC
A(end,end) = A(end,end)-Ds/(h*h); %BC


end