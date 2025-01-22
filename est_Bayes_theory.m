function [mu,varth]=est_Bayes_theory(th0, par, dovarQ=1)
% mean of bayesian estimator, used to return bias.

    dth     = 0.03; %0.02; % numerical integration step.
    th_ar   = dth:dth:2*pi-dth; % integrate var
    ft      = ftun(1:par.n,th0, par); % response at true angle  1xN
    f       = ftun(1:par.n,th_ar, par); % response at candidate nth x N

    if par.poissonQ==0
        logP    = sum((ft-f).^2,2);  %sum over neurons
        P       = exp(-logP/(4*par.si^2));
    else %poissonQ
        Tpois=tpoisfun(par.si);
        ft *= Tpois; f *= Tpois;
        %P=exp(-ft-f).*besseli(0,2*sqrt(ft.*f));
        % prevent large besseli and complex P for long Tpois, see help(bessel)
        P=exp(2*sqrt(ft.*f)-ft-f).*besseli(0,2*sqrt(ft.*f),1);

        P=prod(P,2);
    end
    P = P'/sum(P);

    if par.circQ
        % formally the angular mean equ. below is wrong, see notes
        mu  = atan2mvr(sum(sin(th_ar).*P),sum(cos(th_ar).*P));
    else
        mu  = sum(th_ar.*P); % in Euclidean space
    end

    MLQ=0;
    if MLQ
        printf("Doing ML in Bayes!\n")
        [~,im]=max(P);
        mu=th_ar(im);
    end

    if par.doBayesvarQ ==0
        varth=0;
        return
    end


    % for the variance we need P(th1,th2,TH), see notes
    % make (th1,th2) grid
    nth = numel(th_ar);
    f1  = repmat(f',[1,1,nth]); % 4 x nth x nth, identical in last dim.
    f2  = permute(f1,[1,3,2]); %4 x nth x nth, identical in 2nd dim.

    if par.poissonQ
        Q = ones (nth,nth);
        for i=1:par.n
            prodmat =squeeze(ft(i)*f1(i,:,:).*f2(i,:,:));
            ee =squeeze(-ft(i)-f1(i,:,:)-f2(i,:,:));
            Q .*= hypgeom_02(prodmat).* exp(ee);  % Pii=hyp.* exp(ee); P.*= Pii;
        end
    else
        logQ    = squeeze(sum((ft'-f1).^2,1)+ sum((ft'-f2).^2,1)+ sum((f1-f2).^2,1) );
        Q       = exp(-logQ/(3*2*par.si^2));  % need to fix factor n
    end

    Q    =  Q/sum(sum(Q));
    varth  = sum(sum((th_ar-mu)'*(th_ar-mu).*Q));
    % this fits much better than
    % varth =sum(sum(th_ar'*th_ar.*Q))-mu^2
end

% marx own implementation of hypgeom. See MMA hypgeom02.nb
% rapidly converging sum.
% with Stirling approx, sum bound could be made more precise
function out = hypgeom_02(inmat)
    mm=max(max(inmat));
    if mm >1e5
        error('hypgeom_02: adjust upper bound')
    end
    kmax=100;

    if mm<50000
        kmax = 50;
    end

    if mm<1000
        kmax = 20;
    end
    if mm<500
        kmax = 10;
    end
    if mm<40
        kmax = 5;
    end
    if mm<5
        kmax = 3;
    end

    out=0*inmat;
    for k = 0:kmax
       out += inmat.^k/factorial(k)^3;
    end
end
