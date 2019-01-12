function value = SplineBase(i, k, t, T)
if(k == 1)
    if(t >= T(i) && t < T(i+1))
        value = 1;
    else
        if(T(i+1) == T(length(T)) && T(i) < T(i+1) && t == T(i+1))
            value = 1;
        else
            value = 0;
        end
    end
else
    dt1 = T(i+k-1) - T(i);
    dt2 = T(i+k)-T(i+1);
    if(dt1 == 0)
        w1 = 0;
    else
        w1 = (t-T(i))/dt1;
    end
    if(dt2 == 0)
        w2 = 0;
    else
        w2 = (T(i+k)-t)/dt2;
    end
    if(w1 == 0 && w2 == 0)
        value = 0;
    else
        if w2 == 0
            value = w1 * SplineBase(i, k-1, t, T);
        else
            if w1 == 0
                value =  w2 * SplineBase(i+1, k-1, t, T);
            else
                value = w1 * SplineBase(i, k-1, t, T) + w2 * SplineBase(i+1, k-1, t, T);
            end
        end
    end
end
