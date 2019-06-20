function dH_dtheta=dHamiltonian_dtheta(N,theta,phi,F,A,lam)
%This function calculates the derivative of Hamiltonian with respect to parameter theta in a matrix form.

    up=[0,0;0,1];

    Sp=[0;1]*[1,0];
    Sm=Sp';

    dH1=(F/pi)* (kron(eye(N),up));
    dH2=((2*A)/pi)* kron(diag((1:N)-(N+1)/2 +theta/pi),up);

    Hmv1= kron(eye(N),(Sp*exp(-1i*phi)+ Sm*exp(1i*phi)));  %Microwave 1 which connects l to l
    Hmv2= kron(circshift(eye(N),-1),Sp*exp(-1i*phi)) + kron(circshift(eye(N),1),Sm*exp(1i*phi)); %Microwave 1 which connects l to l-1

    omega=exp((pi/(4*lam))^2);


    dH_dtheta=dH1 + dH2 +  omega*((-theta/(2*lam^2))*exp(-(theta/(2*lam))^2) * Hmv1 + ((pi-theta)/(2*lam^2))* exp(-(((theta-pi)/(2*lam)))^2) * Hmv2);
end