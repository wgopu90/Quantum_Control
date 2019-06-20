function H=Hamiltonian(N,Jx,Jz) 

%Define the Hamiltonian here
% N is the number of lattice sites.
% Jx is the array of Jx values for nearest neighbors. Jz is the array of Jz
% values.
    Hx=circshift(diag(Jx),1)+circshift(diag(Jx),1)';
    Hz=diag(Jz);
    H=Hx+Hz;

end