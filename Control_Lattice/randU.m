function U = randU(n)

% This function generates a Haar random unitary matrix
%
%    U = RANDU(N) generates a random N x N unitary matrix,
%    distributed uniformly according to the Haar measure.


    X = (randn(n) + 1i*randn(n))/sqrt(2);
    [Q,R] = qr(X);
    R = diag(diag(R)./abs(diag(R)));
    U = Q*R;
end