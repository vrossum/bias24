par=struct();  % structure to keep all settings

par.circQ=1;
par.one_dQ=0;
if (par.circQ)
    par.n = 4; % # neurons
    par.phi_ar =([1:par.n]-1)/par.n*2*pi+pi/4; % preferred angles of the neurons
else % not checked
    par.phi_ar=-10+pi/4:2*pi/4:10;
   %  phi_ar=pi/4:2*pi/4:2*pi;
   %  phi_ar =([1:n]-1)/n*2*pi+pi/4;
    par.n=numel(phi_ar);
end

par.b = +0.; % offset, if positive-> narrower tuning
par.poissonQ    =0;
par.vonmisesQ   =0;

par.wid=0.4; % for vonMises and Gaussian
par.si=0.1; % std of added Gaussian noise

% this was confusing...Now we alwyas set Tpois indirectly via tpoisfun(si)
#par.Tpois=1/(par.si)^2; % ^2 since 22/6/2024


par.doBayesvarQ =1;
par.ampl =1;
