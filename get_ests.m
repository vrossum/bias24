function [biasML sdML biasPV sdPV biasBA sdBA biasBAth sdBAth sd_fi]=get_ests(th0, ntr, par)
% given true value th0, estimate theta with various estimators.

    fth     = ftun(1:par.n, th0, par); % av tuning curves at th0
    if par.poissonQ==0
        for itr=1:ntr % run over trials and estimate each time
            r=fth+par.si*randn(1,par.n);  % could rectify after noise, but that makes ML hard
            th_est_ml(itr) = th_mlgrid(r, th0, par); % about 10x faster
           % th_est_ml(itr) = th_ml(r, th0); % use true value as start for ML search.
            th_est_pv(itr) = th_pv(r,par);
            th_est_BA(itr) = est_Bayes(r,th0,par);
        end

    else %poissonQ
        Tpois=tpoisfun(par.si);
        for itr=1:ntr % run over trials and estimate each time
            %r=poissrnd(fth*Tpois);     % or randp...
            r= randp(fth*Tpois);     % or
            th_est_ml(itr) = th_mlgrid(r, th0, par);
            th_est_pv(itr) = th_pv(r, par);
            th_est_BA(itr) = est_Bayes(r, th0, par);
        end
    end

    [thmean_BAth, thvar_BAth] = est_Bayes_theory(th0,par);

    #keyboard()

    biasML  = mean_fun(th_est_ml,par)-th0; % bias
    sdML    = std_fun(th_est_ml,par);

    biasPV  = mean_fun(th_est_pv,par)-th0;
    sdPV    = std_fun(th_est_pv,par);

    biasBA  = mean_fun(th_est_BA,par)-th0;
    sdBA    = std_fun(th_est_BA,par);


    biasBAth  = thmean_BAth-th0;
    sdBAth = sqrt(abs(thvar_BAth)).*sign(thvar_BAth);
    % negative var means there is a bug
    % check if this correct for circular.

    sd_fi   = 1/sqrt(get_fi(th0,par)); % get min. sd from Fisher Info, not bias corrected
end


function th_est=th_ml(r,th_init,par) % maxlikelihood estimator
	if (par.poissonQ==0)
        th_est = fminsearch(@(th)sum((r-ftun(1:n,th)).^2),th_init);
    else
        Tpois=tpoisfun(par.si);
        eps=1e-5;
        th_est = fminsearch(@(th)sum(-r.*log(eps+ftun(1:par.n,th))+Tpois*ftun(1:par.n,th)),th_init);
    end
end


function th_est=th_mlgrid(r, th_init,  par, ngrid=10000) % maxlikelihood estimator

	%persistent
	th_candidates =linspace(0,2*pi,ngrid); % static for efficiency
	if par.poissonQ==0
        mse_cand = sum((r-ftun(1:par.n, th_candidates,par)).^2,2);
        [~, i] = min(mse_cand);
        th_est = th_candidates(i);
    else
        Tpois=tpoisfun(par.si);
        eps=1e-5;
        mse_cand = sum(-r.*log(eps+Tpois*ftun(1:par.n,th_candidates,par))+Tpois*ftun(1:par.n,th_candidates, par),2);
        %keyboard()
        [~, i] = min(mse_cand);
        th_est = th_candidates(i);
    end
end

function th_est=th_pv(r, par) % pop vector
	if par.circQ==0
        error('stop')   % should do CoM here.
    else
        th_est = atan2mvr( sum(sin(par.phi_ar).*r) , sum(cos(par.phi_ar).*r));
    end
end

function thBA = est_Bayes(r, th0, par) % bayesian estimator for given trial
    dth     = 0.005;
    th_ar   = 0:dth:2*pi-dth;
    if (par.poissonQ==0)
        P = exp(sum( -(r-ftun(1:par.n,th_ar,par)).^2,2)/2/par.si^2);
    else
        Tpois=tpoisfun(par.si);
        ft = Tpois*ftun(1:par.n,th_ar,par);
        P=0*th_ar'+1;
        for i=1:par.n
           % logP += r(i)*log(ft(:,i))-ft(:,i)-lgamma(r(i)+1) % gives log(0) trouble
           P .*= ft(:,i).^r(i).*exp(-ft(:,i))/gamma(r(i)+1);
        end
    end

    P=P'/sum(P);
    if par.circQ  % normalization does not matter
        thBA  =  atan2mvr(sum(sin(th_ar).*P),sum(cos(th_ar).*P));
    else
        thBA   = sum(th_ar.*P);
    end
end

% utilities for circular stats.

function mu=mean_fun(data_ar,par) % circular mean of an array
    if par.circQ
        mu=atan2mvr(mean(sin(data_ar)),mean(cos(data_ar)));
    else
        mu=mean(data_ar);
    end
end

function varout = var_fun(data_ar,par)
    if par.circQ
        varout = 1 -norm( mean(cos(data_ar)), mean(sin(data_ar)) );
        varout = 2*varout; % make var same as linear one in low noise limit.!
    else
        varout= var(data_ar);
    end
end

function std=std_fun(data_ar,par)
    if par.circQ
    % see Wikipedia, circular std != sqrt circ var2
    % this one has the right limit % for finite n, std is limited...
        std=norm( [mean(cos(data_ar)), mean(sin(data_ar))]);
        std=sqrt(-2*log(std));
    else
        std=std(data_ar);
    end
end


function out=df(i,th, par) % derivative of tuning curve w.r.t. th.
	if par.circQ==0
        error('stop')
    else
        if par.vonmisesQ==1
            out = -exp( (cos(th'- par.phi_ar(i))-1)/par.wid) .* sin(th'-par.phi_ar(i))/par.wid;
        else
        % CHEKC THIS
            out = -sin(th'-par.phi_ar(i)).*( cos(th'-par.phi_ar(i))-par.b > 0)/(1-par.b);
        end
    end
end

function fi=get_fi(th, par) % get Fisher info given stimulus angle and noise
	if par.circQ==0
        error('stop')
    end

    if par.poissonQ==0 % BUG there was a square missing until 16/3/20
        fi = sum(arrayfun(@df,1:par.n,th,par).^2)/par.si^2;
    else
        eps=1e-5; % prevent df/f = 0/0 error
        Tpois=tpoisfun(par.si);
        fi = Tpois* sum(arrayfun(@df,1:par.n,th,par).^2./(ftun(1:par.n, th, par)+eps));
    end
end
