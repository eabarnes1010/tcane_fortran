function pq = Pq(q)
    pq = exp(1/4)/sqrt(8*pi) * (Kv((q+1)/2) + Kv((q-1)/2));
end

function kv = Kv(v)
    I = (abs(v) < 1e-4);        
    kv(I) = 1.541506751248303 + 1.491844425771660 * v(I).^2;

    J = (abs(v-1) < 1e-4);
    kv(J) = 3.747025974440712 + 6.166027005560792 * (v(J)-1);

    K = ~I & ~J;
    kv(K) = pi/2 * (Iv(-v(K)) - Iv(v(K))) ./ sin(v(K)*pi);
end

function iv = Iv(v)
    term = 1 ./ (8.^v .* gamma(v + 1));
    iv  = term;

    for m = 1:6
        term = term ./ (64 * m * (v+m));
        iv = iv + term;
    end
end