function H=Hamiltonian(N,theta,phi,F,A,lam) 

%Define the Hamiltonian here
% N is the number of lattice sites.
% theta is the angle of polarization, Phi is the phase of microwave, F is
% the strength of linear gradient and A is the curvature of quadratic
% potential. Omega is the strength of microwave which is set to 1. lam is the Lamb-Dicke
% parameter.

    dw=[1,0;0,0];
    up=[0,0;0,1];

    H1=kron(diag((1:N)-(N+1)/2),dw) + kron(diag((1:N)-(N+1)/2 +theta/pi),up); %Linear Gradiant
    H2=kron(diag((1:N)-(N+1)/2).^2,dw) + kron(diag((1:N)-(N+1)/2 +theta/pi).^2,up); %Quadratic potential

    Sp=[0;1]*[1,0];
    Sm=Sp';

    Hmv1= kron(eye(N),(Sp*exp(-1i*phi)+ Sm*exp(1i*phi)));  %Microwave 1 which connects l to l
    Hmv2= kron(circshift(eye(N),-1),Sp*exp(-1i*phi)) + kron(circshift(eye(N),1),Sm*exp(1i*phi)); %Microwave 1 which connects l to l-1

    omega=exp((pi/(4*lam))^2);

% Return the hamiltonian value

    H=F*H1 +  A* H2 + omega*( exp(-(theta/(2*lam))^2) *  Hmv1 + exp(-(((theta-pi)/(2*lam)))^2) * Hmv2);
end