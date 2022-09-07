function omega = flow_around_cylinder_unsteady
Re=60;
%%%%% define the grid %%%%%
n=101; m=202; % number of grid points
N=n-1; M=m-2; % number of grid intervals: 2 ghost points, theta=-h,2*pi
h=2*pi/M; % grid spacing based on theta variable
xi=(0:N)*h; theta=(-1:M)*h; % xi and theta variables on the grid
%%%%% time-stepping parameters %%%%%
t_start=0; t_end=0.5; %vortex street starts at around t=1000
tspan=[t_start t_end]; 
%%%%% Initialize vorticity field %%%%%
omega=zeros(n,m);
%%%%% Construct the matrix A for psi equation %%%%%
% boundary_index = [bottom, left, top, right]
boundary_index_bottom = 1:n;
boundary_index_left = 1:n:1+(m-1)*n;
boundary_index_top = 1+(m-1)*n:m*n;
boundary_index_right = n:n:m*n;

diagonals1 = [4*ones(m*n,1), -1*ones(m*n,4)];
A=spdiags(diagonals1,[0 -1 1 -n n], m*n, m*n); %interior points
I=speye(m*n);
A(boundary_index_left,:)=I(boundary_index_left,:);
A(boundary_index_right,:)=I(boundary_index_right,:);
diagonals2 = [ones(m*n,1), -1*ones(m*n,2)];
Aux = spdiags(diagonals2,[0 m*(n-2) -m*(n-2)], m*n, m*n);
A(boundary_index_top,:)=Aux(boundary_index_top,:);
A(boundary_index_bottom,:)=Aux(boundary_index_bottom,:);

%%%%% Find the LU decomposition %%%%%
[L,U]=lu(A); clear A;
%%%%% Compute any time-independent constants %%%%%

%%%%% advance solution using ode23 %%%%%
options=odeset('RelTol', 1.e-03);
omega=omega(2:n-1,2:m-1); % strip boundary values for ode23
omega=omega(:); % make a column vector
[t,omega]=ode23...
  (@(t,omega)omega_eq(omega,L,U,theta,xi,h,Re,n,m),tspan, omega, options);
%%%%% expand omega to include boundaries %%%%%
temp=zeros(n,m);
temp(2:n-1,2:m-1)=reshape(omega(end,:),n-2,m-2);
omega=temp; clear temp;
%%%%% compute stream function (needed for omega boundary values) %%%%%
omega_tilde = zeros(n,m);

omega_tilde(1,:)=0;%left boundary 
omega_tilde(n,:)=exp(xi(n))*sin(theta);%right boundary
omega_tilde(:,1)=0;%bottom boundary
omega_tilde(:,m)=0;%top boundary
for i=2:n-1
    for j=2:m-1                   
        omega_tilde(i,j) = h^2*exp(2*xi(i))*omega(i,j);%interior points     
    end    
end  

omega_tilde = omega_tilde(:);

psi = reshape(U\(L\omega_tilde),n,m);

%%%%% set omega boundary conditions %%%%%
omega(1,:)= (psi(3,:) - 8*psi(2,:)) * 1/(2*h^2);
omega(n,:)=0;    
omega(:,1)=omega(:,m-1);           
omega(:,m)=omega(:,2);  

%%%%% plot scalar vorticity field %%%%%
%plot_Re60(omega);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_omega_dt=omega_eq(omega,L,U,theta,xi,h,Re,n,m)
%%%%% expand omega to include boundary points %%%%%
temp=zeros(n,m);
index1=2:n-1; index2=2:m-1;
temp(index1,index2)=reshape(omega,n-2,m-2);
omega=temp; clear temp;
%%%%% compute stream function %%%%%
omega_tilde = zeros(n,m);

omega_tilde(1,:)=0;%left  
omega_tilde(n,:)=exp(xi(n))*sin(theta);%right
omega_tilde(:,1)=0;%bottom
omega_tilde(:,m)=0;%top
for i=2:n-1
    for j=2:m-1                   
        omega_tilde(i,j) = h^2*exp(2*xi(i))*omega(i,j);%interior      
    end    
end  

omega_tilde = omega_tilde(:);

psi = reshape(U\(L\omega_tilde),n,m);

%%%%% compute derivatives of omega %%%%%

omega(1,:)= (psi(3,:) - 8*psi(2,:)) * 1/(2*h^2);
omega(n,:)=0;    
omega(:,1)=omega(:,m-1);           
omega(:,m)=omega(:,2);   

d_omega_dt=zeros(n-2,m-2);
for i=2:n-1            
    for j=2:m-1         
        d_omega_dt(i-1,j-1)=...
        2*exp(-2*xi(i))* 1/(h^2*Re)*...
        (omega(i+1,j)+omega(i-1,j)+...
         omega(i,j+1)+omega(i,j-1)-4*omega(i,j))...                 
        +exp(-2*xi(i))/(4*h^2)*...
        ((psi(i+1,j)-psi(i-1,j))*(omega(i,j+1)-omega(i,j-1))-...
         (psi(i,j+1)-psi(i,j-1))*(omega(i+1,j)-omega(i-1,j)));              
   end     
end  
  d_omega_dt = d_omega_dt(:);
end

