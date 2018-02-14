%Similar function to ramp, but adds a minimum strain instead of simply 1
function k = ramp2(strainmin,strainmax,reps,reptime)


    function ms = rampfn(xs)
    for i = 1:length(xs)
        x = xs(i);
        if x < reptime*(strainmin-1)/(strainmax-strainmin)
            m = 1+(strainmax-strainmin)/reptime*x;
        else
            x = x-reptime*(strainmin-1)/(strainmax-strainmin);
            scalepos = floor(x/reptime);
            if x>reps*reptime
                if mod(reps,2)==0
                    m = strainmin;
                else
                    m = strainmax;
                end
            else
                x = x-scalepos*reptime;
                x = x/reptime;
                if mod(scalepos,2) ==1
                    m = strainmin+(strainmax-strainmin)*(1-x);
                else
                    m = strainmin+(strainmax-strainmin)*x;
                end
            end
        end
        ms(i)=m;
    end
    end
    k = @rampfn;
end




