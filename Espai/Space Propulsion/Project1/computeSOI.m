function rSOI = computeSOI(body, primary, a)
    rSOI = a * (body.mass / primary.mass)^(2/5);
end
