%Ramp function, largely to be used as an argument to the strainfunc
%argument in solvers for applied strain such as for the stress relaxation
%experiment. 
function k = ramp(strainmax,reps,reptime)
    function ms = rampfn(xs)
    for i = 1:length(xs)
        x = xs(i);
        scalepos = floor(x/reptime);
            if x>reps*reptime
                if mod(reps,2)==0
                    m = 1;
                else
                    m = strainmax;
                end
            else
                x = x-scalepos*reptime;
                x = x/reptime;
                if mod(scalepos,2) ==1
                    m = 1+(strainmax-1)*(1-x);
                else
                    m = 1+(strainmax-1)*x;
                end
            end
            ms(i)=m;
    end
    end
    k = @rampfn;
end




