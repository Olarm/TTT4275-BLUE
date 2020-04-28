function est = BLUE_33(signal, H, C)
    M = ((H' * (C \ H)) \ H') / C;
    est = M * unwrap(angle(signal));
end