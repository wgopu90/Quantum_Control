function dH_dphi=dHamiltonian_dphi(N,theta,phi,F,A,lam)
%This function calculates the derivative of Hamiltonian with respect to parameter phi in a matrix form.
    Sp=[0;1]*[1,0];
    Sm=Sp';

    dHmv1= -1i * kron(eye(N),(Sp*exp(-1i*phi)- Sm*exp(1i*phi)));  %Microwave 1 which connects l to l
    dHmv2= -1i * (kron(circshift(eye(N),-1),Sp*exp(-1i*phi)) - kron(circshift(eye(N),1),Sm*exp(1i*phi))); %Microwave 1 which connects l to l-1

    omega=exp((pi/(4*lam))^2);


    dH_dphi=omega*(exp(-(theta/(2*lam))^2) * dHmv1 + exp(-(((theta-pi)/(2*lam)))^2) * dHmv2);
end