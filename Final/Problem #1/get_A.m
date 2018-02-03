function A = get_A(Nx,Ny,dt,dx,dy,lambda)
    BOUNDARY_MASK = get_border(Nx,Ny);
    BOUNDARY_MASK = reshape(BOUNDARY_MASK.',[],1); %creates a vector that
                                                   %marks the positions
                                                   %where boundary
                                                   %conditions apply
    Identity_Matrix = speye(Nx*Ny,Nx*Ny);          %sparse identity matrix
    
    a_b = -lambda*dt/(dy^2);
    a_l = -lambda*dt/(dx^2);
    a_c = (1 + 2*lambda*dt*(1/(dx^2) + 1/(dy^2)));
    a_r = a_l;
    a_t = a_b;

    %matrix B will contain the values to be used to fill diagonals
    B = zeros(Nx*Ny,5);
    B(:,1) = a_b;
    B(:,2) = a_l;
    B(:,3) = a_c;
    B(:,4) = a_r;
    B(:,5) = a_t;
    
    %diagonals tells spdiags where to put the values from each column of B
    %along the diagonals of the A matrix. 0 correspsonds to the largest diagonal.
    %+- values are the diagonals above/below the largest diagonal.
    diagonals = [-Nx -1 0 1 Nx];
    Internal = spdiags(B, diagonals, Nx*Ny, Nx*Ny);
    
    %Where boundary conditions apply, we simply want the identity matrix row.
    %Elsewhere, we want to keep the values imposed on A by spdiags.
    A = BOUNDARY_MASK.*Identity_Matrix + ~BOUNDARY_MASK.*Internal;
end