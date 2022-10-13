function kv = Kvq(V)
    kv = arrayfun(@(v) elemental_Kv(v), V);
end

function kv = elemental_Kv(v)
    TOLERANCE = 1e-4;

    if abs(v - round(v)) > TOLERANCE
        kv = pi/2 * (Iv(-v) - Iv(v)) ./ sin(v*pi);
    elseif abs(v) < TOLERANCE
        kv = 1.541506751248303 + 1.491844425771660 * v^2;
    elseif abs(v-1) < TOLERANCE
        kv = 3.747025974440712 + 6.166027005560792 * (v-1);
    else
        nu = 1 + v - round(v);

        kv_minus1 = 1.541506751248303 + 1.491844425771660 * (v - round(v))^2;
        kv        = 3.747025974440712 + 6.166027005560792 * (v - round(v));

        for j = 1 : int16(round(v) - 1)
            kv_plus1 = kv_minus1 + 8.0*nu*kv;

            nu = nu + 1.0;
            kv_minus1 = kv;
            kv = kv_plus1;
        end
    end
end
