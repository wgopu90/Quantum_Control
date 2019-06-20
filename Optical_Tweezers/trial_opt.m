function trial_opt(N,savefile)

    % This function does the optimization and find the corresponding  waveform 'x'.
    % Fidelity and the waveform is saved in the file savefile.
    % It calls function 'Hamiltonian' and 'fidelity'
    % N is the number of lattice sites ( or system size).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    num_avg=100;

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

                                      

    tic
    fidelityk=0;
    
    parfor kk=1:num_avg
        
        
        casefidelity=2;
         % 1 for state to state mapping and 2 for Unitary
        par0=zeros(N,N,2);
        %this line is just to initialize the temporary variable, so that it doesn't create an error in parfor
    
    
        dt=2*pi;


        a=eye(N);
        lat=1; %which lattice site is the atom initally in.
        initial_state=kron(a(:,lat),[1;0]);
    
    
        rng shuffle;
        Jx=zeros(N,N);
        Jz=zeros(N,N);
      
        
        target_unitary=randU(N);
         target_state=(randn_target_state(N));
        par0(:,:,1)=Jx;
        par0(:,:,2)=Jz;
        
        fun=@(par) fidelity(initial_state,target_state,target_unitary,dt,N,par,casefidelity);

        [x,fval] = fminunc(fun,par0,options);
        
        fidelityk=fidelityk+1-fval;
    end
    
    fidelity1=fidelityk/num_avg;                    % Average fidelity
    time_simulation=toc;                          % time taken for simulation
    save(savefile,'fidelity1','x','time_simulation')
end
