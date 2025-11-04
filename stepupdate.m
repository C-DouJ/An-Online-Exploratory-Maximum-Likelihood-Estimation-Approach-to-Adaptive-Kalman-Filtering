function QRestimator = stepupdate(QRestimator)
    for i = 1:length(QRestimator.rate)
        if(QRestimator.rate(i)>1)
            QRestimator.rate(i) = max([QRestimator.minrate,QRestimator.ratedr*QRestimator.rate(i)]);
        else
            QRestimator.rate(i) = min([1/QRestimator.minrate,1/QRestimator.ratedr*QRestimator.rate(i)]);
        end    
    end
end