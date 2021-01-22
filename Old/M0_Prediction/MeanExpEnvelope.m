function mean = MeanExpEnvelope(amplitude,relaxTime,t_0,t_end)
    mean = -relaxTime*amplitude*(exp(-t_end/relaxTime)-exp(-t_0/relaxTime));
end

