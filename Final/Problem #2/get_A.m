function A = get_A(Nx,Ny,dt,dx,dy,lambda,vy)
  a_b = -lambda*dt/(dy^2);
  a_l = -lambda*dt/(dx^2);
  a_c = 1 + 2*lambda*dt*(1/dx^2 + 1/dy^2);
  a_r = a_l;
  a_t = a_b;
  
  %robin boundary condition values for lhs
  a_l_rob = -lambda*dt/(dx^2);
  a_c_rob = a_c + 2*vy*dt/dy;
  a_r_rob = a_l_rob;
  a_t_rob = a_t + a_b;

  %ROBIN_BC marks where Robin boundary conditions apply
  ROBIN_BC = zeros(Ny,Nx);
  ROBIN_BC(1,2:(Nx-1)) = 1;
  
  %marks the (non-Robin) boundary condition locations.
  BOUNDARY_MASK = get_border(Nx,Ny);
  BOUNDARY_MASK = xor(BOUNDARY_MASK,ROBIN_BC);
  
  ROBIN_BC = reshape(ROBIN_BC.',[],1);
  BOUNDARY_MASK = reshape(BOUNDARY_MASK.',[],1);
  
  %sparse identity matrix
  Identity_Matrix = speye(Nx*Ny,Nx*Ny);

  %insert internal grid point coefficients into matrix
  B = zeros(Nx*Ny,5);
  B(:,1) = a_b;
  B(:,2) = a_l;
  B(:,3) = a_c;
  B(:,4) = a_r;
  B(:,5) = a_t;
  diagonals = [-Nx -1 0 1 Nx];
  Internal = spdiags(B, diagonals, Nx*Ny, Nx*Ny);

  %insert robinson boundary condition coefficients into matrix
  B = zeros(Nx*Ny,4);
  B(:,1) = a_l_rob;
  B(:,2) = a_c_rob;
  B(:,3) = a_r_rob;
  B(:,4) = a_t_rob;
  diagonals = [-1 0 1 Nx];
  Robinson = spdiags(B, diagonals, Nx*Ny, Nx*Ny);

  A = BOUNDARY_MASK.*Identity_Matrix...
      + ROBIN_BC.*Robinson...
      + ~(BOUNDARY_MASK | ROBIN_BC).*Internal;
end