clear
N=5;
n_steps=4*N^2;
F=0;
A=1;
lam=0.4;
dt=0.1;




theta=randn(1,n_steps)*2*pi;
phi=randn(1,n_steps)*2*pi;

w0=10000;
dtt=0.001;
dw=0.01;
r=dt/dtt;
t=0:dtt:n_steps*dt;
w=0:dw:w0;
theta1=zeros(1,length(t));
phi1=zeros(1,length(t));
for i=1:n_steps
    for j=1:r
        theta1((i-1)*r+j)=theta(i);
        phi1(((i-1)*r+j))=phi(i);
    end
end
disp('done 1')
theta1(end)=theta(end);
phi1(end)=phi(end);
P1=pws(theta,w,dt);
P2=pws(phi,w,dt);
j=0:n_steps-1;
disp('done2')
thetay=zeros(1,length(t));
phiy=zeros(1,length(t));
for k=1:length(t)
    f1=dt*(j+1)-t(k);
    f2=dt*j - t(k);
    q=(1/(pi))*(sinint(w0*f1)-sinint(w0*f2));
    
    thetay(k)=theta*q';
    phiy(k)=phi*q';
    k
end
plot(t,[theta1;thetay])


U=eye(2*N);
U1=eye(2*N);

    for ii=length(t):-1:1
        H=Hamiltonian(N,theta1(ii),phi1(ii),F,A,lam);
        H1=Hamiltonian(N,thetay(ii),phiy(ii),F,A,lam);
        U=U*expm(-1i*dtt*H);
        U1=U1*expm(-1i*dtt*H1);
        ii
    end
    
    (1/n_steps)*abs(trace(U'*U1)).^2