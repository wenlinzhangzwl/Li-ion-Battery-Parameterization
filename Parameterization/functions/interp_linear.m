function y = interp_linear(x, y, x0)
    ind = find(x>=x, 1); 
    m = (y2-y1) / (x2-x1); 
    y = y1 + m * (x-x1);
end

