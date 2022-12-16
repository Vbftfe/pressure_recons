function rms = calcRms(pivData, reconsField, trueField)
    delta_p_max = max(trueField) - min(trueField);
    error = ((reconsField(pivData.domain) - trueField(pivData.domain))/delta_p_max).^2;
    len = sum(pivData.domain);
    rms = sqrt(sum(error(:))/len)*100;