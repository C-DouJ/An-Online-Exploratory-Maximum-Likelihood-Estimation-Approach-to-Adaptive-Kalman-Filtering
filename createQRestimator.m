function QRestimator = createQRestimator(mode,updatemode,winlen,delaylen,xf_SCDed,Pf_SCDed)
    QRestimator.rate = [1.5 1.5 1.5 1.5];
    QRestimator.minrate = 1.01;
    QRestimator.ratedr = 0.999;
    QRestimator.adjidx = 1;
    QRestimator.winlen = winlen;
    QRestimator.measurenum = 0;
    QRestimator.mode = mode;
    QRestimator.updatemode = updatemode;
    if  strcmp(QRestimator.mode,'fix')
        QRestimator.xf_history = xf_SCDed;
        QRestimator.Pf_history = Pf_SCDed;
        QRestimator.zhistorys = zeros(2,QRestimator.winlen);
        QRestimator.uhistorys = zeros(1,QRestimator.winlen);
        QRestimator.totallen = QRestimator.winlen;
    elseif strcmp(QRestimator.mode,'delay')
        QRestimator.delaylen = delaylen;
        statedim = length(xf_SCDed);
        QRestimator.xf_historys = zeros(statedim,QRestimator.winlen+QRestimator.delaylen);
        QRestimator.Pf_historys = zeros(statedim,statedim,QRestimator.winlen+QRestimator.delaylen);
        QRestimator.zhistorys = zeros(2,QRestimator.winlen+QRestimator.delaylen);
        QRestimator.uhistorys = zeros(1,QRestimator.winlen+QRestimator.delaylen);
        QRestimator.totallen = QRestimator.winlen+QRestimator.delaylen-1;
    end
end