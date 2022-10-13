function pq = Pq(q)
    pq = exp(1/4)/sqrt(8*pi) * (Kv((q+1)/2) + Kv((q-1)/2));
end

function kv = elemental_Kv(v)
    TOLERANCE = 1e-4;

    if abs(v - round(v)) < TOLERANCE
        

    if abs(v) < 1e-4
        kv = 1.541506751248303 + 1.491844425771660 * v^2;
    elseif abs(v-1) < 1e-4
        kv = 3.747025974440712 + 6.166027005560792 * (v-1);
    else
        kv = pi/2 * (Iv(-v) - Iv(v)) ./ sin(v*pi);
    end
end

function kv = Kv(V)
    kv = arrayfun(@(v) elemental_Kv(v), V);
end


function iv = Iv(v)
    term = 1 ./ (8.^v .* gamma(v + 1));
    iv  = term;

    for m = 1:6
        term = term ./ (64 * m * (v+m));
        iv = iv + term;
    end
end