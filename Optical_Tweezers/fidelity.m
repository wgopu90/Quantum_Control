function [Fidel,Gradient]=fidelity(initial_state,target_state,target_unitary,dt,N,J,casefidelity)

    %This function calculates both fidelity and gradient of fidelity with
    %hamiltonian given by the function Hamiltonian(...).
    %casefidelity=1 ---> state to state mapping
    %Casefidelity=2 ---> Unitary mapping.


    Jx=J(:,:,1);
    Jz=J(:,:,2);
    num=size(Jx,2);

    U1=zeros(N,N,num);
    U=eye(N);

    for ii=num:-1:1
        H=Hamiltonian(N,Jx(:,ii)',Jz(:,ii)');
        [V,L]=eig(H);
        [L1,L2]=meshgrid(diag(L));
        prefac=(-1i*dt) .* exp(-1i*(dt/2)* (L1+L2)) .* sinc((dt/(2*pi))* (L1-L2));
        Id=eye(N);
        for q=1:size(Jx,1)
            Ids=Id(:,q)*Id(:,q)';
            dH_dJx=circshift(Ids,1)+circshift(Ids,1)';
            dU_dJx(:,:,q,ii)=V * ((V' * dH_dJx * V).*prefac) * V';
        end
        for q=1:N
            dH_dJz=Id(:,q)*Id(:,q)';
            dU_dJz(:,:,q,ii)=V * ((V' * dH_dJz * V).*prefac) * V';
        end
        U1(:,:,ii)=expm(-1i*dt*H);
        U=U*U1(:,:,ii);
    end
    
    if casefidelity==1
        Fidel=1-((target_state'*U*initial_state)' *(target_state'*U*initial_state));
    elseif casefidelity==2
        a=trace(target_unitary'*U);
        Fidel=1-(1/(N^2))*(abs(a)^2);
    end

    if nargout>1
        Gradient=zeros(N,num,2);
        for q=1:N
        for ii=1:num   
            UU1=eye(N);
            UU2=eye(N);
            for jj=num:-1:1
                if jj==ii
                     UU1=UU1*dU_dJx(:,:,q,jj);
                     UU2=UU2*dU_dJz(:,:,q,jj);
                else
                    UU1=UU1*U1(:,:,jj);
                    UU2=UU2*U1(:,:,jj);
                end
            end
            if casefidelity==1
                Gradient(q,ii,1)= -2 * real(((target_state'*U*initial_state)')  * (target_state' * UU1 * initial_state));            
                Gradient(q,ii,2)= -2 * real((target_state'*U*initial_state)' * target_state'* UU2*initial_state);                
            elseif casefidelity==2
                Gradient(q,ii,1)= -2 * (1/(N^2))*(real(trace(target_unitary'*UU1)*conj(a)));
                Gradient(q,ii,2)= -2 * (1/(N^2))*(real(trace(target_unitary'*UU2)*conj(a)));
            end
        end
        end
    end

end  


        
