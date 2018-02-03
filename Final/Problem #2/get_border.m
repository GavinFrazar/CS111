function B = get_border(Nx,Ny)
    B = zeros(Ny,Nx);
    B(:,1) = 1; %left
    B(:,Nx) = 1; %right
    B(1,:) = 1; %bottom
    B(Ny,:) = 1; %top
end

