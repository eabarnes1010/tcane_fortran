function ivq = Ivq(v)
    term = 1 ./ (8.^v .* gamma(v + 1));
    ivq  = term;

    for m = 1:10
        term = term./(64 * m * (v+m));
        ivq = ivq + term;
    end
end