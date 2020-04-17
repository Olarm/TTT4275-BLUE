function s = sigma(snr_db)
    [~, w] = size(snr_db);
    s = zeros(1,w);
    
    for i=1:w
        s(i) = 10^(-1/2*log10(2) - 1/20*snr_db(i));
    end
end