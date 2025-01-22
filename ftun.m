function out = ftun(i,th,par)
% mean tuning curves, normalized to peak=1
% df() is not automatically updated!

    persistent firstrunQ=1;

	if par.circQ==1
        if par.vonmisesQ==1
            out = exp((cos(th'-par.phi_ar(i))-1)/par.wid);
        else
            out = max(cos(th'-par.phi_ar(i))-par.b,0)/(1-par.b);
        end
    else % linear space, gaussian coding
        out = exp(-(th'-par.phi_ar(i)).^2/ (2*par.wid^2) );
    end
    out *= par.ampl;

    if par.one_dQ
        if firstrunQ
            printf("ampl override! 1D\n")
        ampl=[1, 0, 0, 0]; % by using individual ampl we can do 1dcase simple
        out.*= ampl(i);

        end
        firstrunQ=0;
    end

    %out=0.3+0.2*th'+par.phi_ar(i)*0; % linear model, for debugging
end
