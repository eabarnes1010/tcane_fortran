function ivm = Ivm(v)
    term = 8^v/gamma(1 - v);
    ivm  = term;

    for m = 1:10
        term = term/(64 * m * (m-v));
        ivm = ivm + term;
    end
end