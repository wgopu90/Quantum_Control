function [Fidel,Gradient]=fidelity(initial_state,target_state,target_unitary,dt,N,phi_theta,F,A,lam,casefidelity)

%This function calculates both fidelity and gradient of fidelity with
%hamiltonian given by the function Hamiltonian(...).
num=length(phi_theta)/2;
phi=phi_theta(1:num);
theta=phi_theta(num+1:end);


%H=zeros(2*N,2*N,num);
U1=zeros(2*N,2*N,num);
U=eye(2*N);

    for ii=num:-1:1
        H=Hamiltonian(N,theta(ii),phi(ii),F,A,lam);
        [V,L]=eig(H);
        [L1,L2]=meshgrid(diag(L));
        prefac=(-1i*dt) .* exp(-1i*(dt/2)* (L1+L2)) .* sinc((dt/(2*pi))* (L1-L2));
        
        dH_dphi=dHamiltonian_dphi(N,theta(ii),phi(ii),F,A,lam);
       
        
        dU_dphi(:,:,ii)=V * ((V' * dH_dphi * V).*prefac) * V';    
        
      
        U1(:,:,ii)=expm(-1i*dt*H);
        U=U*U1(:,:,ii);
    end
    
if casefidelity==1
    Fidel=1-((target_state'*U*initial_state)' *(target_state'*U*initial_state));
elseif casefidelity==2
    a=trace(target_unitary'*U);
    Fidel=1-(1/(4*N^2))*(abs(a)^2);
end

if nargout>1
    Gradient=zeros(num,1);
    for ii=1:num   
        UU1=eye(2*N);
       
        for jj=num:-1:1
            if jj==ii
                UU1=UU1*dU_dphi(:,:,jj);
                
            else
                UU1=UU1*U1(:,:,jj);
               
            end
        end
        if casefidelity==1
            Gradient(ii,1)= -2 * real(((target_state'*U*initial_state)')  * (target_state' * UU1 * initial_state));            
            
        elseif casefidelity==2
            Gradient(ii,1)= -2 * (1/(4*N^2))*(real(trace(target_unitary'*UU1)*conj(a)));
            
        end
    end
end
end  


        
