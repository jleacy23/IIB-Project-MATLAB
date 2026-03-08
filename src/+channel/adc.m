function y = adc(x, ENOBits)
    % reduce resolution to ENOBits
    y = round(x * 2^ENOBits) / 2^ENOBits;
end