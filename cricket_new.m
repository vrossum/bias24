# % version for github
#
# Salinas'94 uses cos tuning with b= -0.14 (this increases bias wrt to b=0, but reduces variance)
# noise from TEunissen & Miller sigma= a0+a1*rate
#

#TODO:
# -check circular stats (and does it matter?)
# - bayes theory
# - poisson noise

# von Mises: attractive (bump-valley); cos repulsive
# but also b parameter in [cos-b]  matters.

close all;
clear all;
pkg load parallel

function plot_popr(th0,par)
    ft = ftun(1:par.n,th0);
    figure(121)
    plot(ft)
end

function plot_ftun(par)
    dth=0.01;
    if par.circQ
        th  = 0:dth:2*pi-dth;
    else
        th = -5:dth:10;
    end    
    
    figure(100)
    hold on
    for i=1:par.n
        ff=ftun(i,th,par);
        str=strcat(";",num2str(i),";");
        plot(th,ff,str)
    end
    xlabel("th");ylabel("tuning curves")
    axis([0 2*pi])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_par();

simcase=1
switch simcase
case(0)
    par.b= 0.707;
    par.vonmisesQ   =0;
    plot_ftun(par);
    
case(1)
    # Here: max and min of bias vs offset. Bayes theory approx only.

    ntr= 100  % # trials, should be big. >1000
    th = pi/2: 0.05: pi*3/4; % half side
    res=[];
    par.poissonQ=0;
    par.doBayesvarQ=1;
    par.vonmisesQ=0;
    for par.si=[0.001 0.05 0.1 0.15 0.2]
        par.si

        res=[]; %good for plot , bad for save...
        for par.b = -1.1:0.02:0.1
            par.b;
            nproc=100;
            [biasML sdML biasPV sdPV biasBA sdBA biasBAth sdBAth sd_fi]= pararrayfun(nproc,@(th) get_ests(th,ntr,par), th,"VerboseLevel",0);
            % [biasML sdML biasPV sdPV biasBA sdBA biasBAth sdBAth sd_fi]= arrayfun(@(th) get_ests(th,ntr,par), th);

            biasres = [th; biasML; biasPV;  biasBA;  biasBAth]';
            sdres = [th; sdML; sdPV ;sdBA ;sdBAth]';

            res=[res [par.b, max(biasML) min(biasML) max(biasPV) min(biasPV) max(biasBA) min(biasBA) max(biasBAth) min(biasBAth)]'];
        end %parameter loop
        figure(1)
        plot(res(1,:),res(2,:),';[cos];')
        hold on
        plot(res(1,:),res(3,:),';[cos];')

        xlabel('offset')
        ylabel('max and min of bias')
    end

    rest=res';
    save bias_minmax.dat rest

endswitch
