function rSOI = computeSOI(m_planet, m_primary, a)
    rSOI = a * (m_planet / m_primary)^(2/5);
end
