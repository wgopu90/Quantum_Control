function trial_opt(N_min,N_max,num_avg,savefile)
    % This function does the optimization and find the corresponding  waveform 'x'.
    % Fidelity and the waveform is saved in the file savefile.
    % It runs for all system sizes from N_min to N_max.
    % It calls function 'Hamiltonian' and 'fidelity'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---Options for fminunc Optimization function
options = optimset(...  % these are the settings from unitary control search
    'TolX',            1e-16,...    % related to minimum tolerance for change in x
    'TolFun',          1e-8,...    % minumum tolerance for change in the function
    'MaxIter',         2000,...   % maximum number of iterations
    'DerivativeCheck', 'off',...    % compare analytic gradient to numerical estimation (off)
    'GradObj',         'on',...         % tells matlab whether the gradient is supplied to fminunc
    'LargeScale',       'off', ... % when turned off will increase calculation speed
    'Display',         'off',...          % output type in matlab command window, USE 'iter' to turn on
    'MaxFunEvals',     10^6,...         % maximum number of function calls'MaxFunEvals', 500);
    'ObjectiveLimit',  0.00001);    % once objective function reaches this value, fminunc ...
                                      %stops (ONLY works for quasi-Newton method)

mean_fidelity=zeros(N_max-N_min+1,1);
timetoc=zeros(N_max-N_min+1,1);


step=4; 

F=0;
A=1;
lam=0.4;

parfor Nprime=1:N_max-N_min+1
    tic
      Casefidelity=2;   % 1 for state to state mapping and 2 for Unitary
    N=Nprime+N_min-1;
    n_steps=step*N;  %this line is just to initialize the temporary variable, so that it doesn't create an error in parfor
    
    
    dt=2*pi;
    if casefidelity==1
        n_steps=step*N;
    elseif casefidelity==2
        n_steps=step*(N^2);
    end
    T=n_steps*dt;



    a=eye(N);
    lat=1; %which lattice site is the atom initally in.
    initial_state=kron(a(:,lat),[1;0]);

    fidelity1=zeros(num_avg,1);
    
    for j=1:num_avg
        rng shuffle;
        phi0=rand(n_steps,1)*2*pi;
        theta0=(pi/2)+ (2*rand(n_steps,1)-1)*lam;
        %theta0=(pi/2)*ones(n_steps,1);
        %theta=pi/2;
        par0=[phi0;theta0];

        %target_state=kron(randn_target_state(N),[1;0]);
        target_unitary=randU(2*N);
         target_state=(randn_target_state(2*N));
        % target_state=kron(a(:,lat+1),[1;0]);
        
        fun=@(par) fidelity(initial_state,target_state,target_unitary,dt,N,par,F,A,lam,casefidelity);

        [x,fval] = fminunc(fun,par0,options);
        fidelity1(j)=1-fval;
    end
    mean_fidelity(Nprime)=mean(fidelity1)
    N
    timetoc(Nprime)=toc
end

save(savefile,'mean_fidelity','N_min','N_max','timetoc')
end
