function [psi, omega] = flow_around_cylinder_steady
Re=10; 
%%%%% define the grid %%%%%
n=101; m=101; % number of grid points
N=n-1; M=m-1; % number of grid intervals
h=pi/M; % grid spacing based on theta variable
xi=(0:N)*h; theta=(0:M)*h; % xi and theta variables on the grid
%%%%% Initialize the flow fields %%%%%
psi=zeros(n,m);
omega=zeros(n,m);
psi(n,:)=exp(xi(n)) * sin(theta(:));% Write the free stream bc here
%%%%% Set relax params, tol, extra variables %%%%%
r_psi=1.8; % Set the relaxation parameter here, psi equation
r_omega=0.9; % Set the relaxation parameter here, omega equation
delta=1.e-08; % error tolerance
error=2*delta; % initialize error variable
%%%%% Add any additional variable definitions here %%%%%
vec_const = h^2*exp(2*xi);
ReDiv8 = Re/8;
%%%%% Main SOR Loop %%%%%
while (error > delta)
    psi_old = psi; omega_old = omega;
    for i=2:n-1
        for j=2:m-1
            psi(i,j)=(1-r_psi)*psi(i,j)+.25*r_psi*(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+...
                vec_const(i) * omega(i,j)); % Write psi equation here
        end
    end
    error_psi=max(abs(psi(:)-psi_old(:)));
    omega(1,:)=(psi(3,:) - 8*psi(2,:))/(2*h^2); % Write the boundary condition here
    for i=2:n-1
        for j=2:m-1
             f=(psi(i+1,j)-psi(i-1,j))*(omega(i,j+1)-omega(i,j-1))-...
              (psi(i,j+1)-psi(i,j-1))*(omega(i+1,j)-omega(i-1,j));
            omega(i,j)=(1-r_omega)*omega(i,j)+.25*r_omega*(omega(i+1,j)+omega(i-1,j)+omega(i,j+1)+omega(i,j-1)+...
                ReDiv8*f); % Write omega equation here
        end
    end
    error_omega=max(abs(omega(:)-omega_old(:)));
    error=max(error_psi, error_omega);
end
plot_Re10(psi);

