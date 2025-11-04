function [Qest,Rest,QRestimator] =  estimateQR(QRestimator,F,H,B,G,Qcur,Rcur)
if QRestimator.measurenum > QRestimator.totallen, display("len error!!!!!!!!!!");end
uhistorys = QRestimator.uhistorys;
zhistorys = QRestimator.zhistorys;
adjidx = QRestimator.adjidx;
rate = QRestimator.rate;
Qexp = Qcur;
Rexp = Rcur;
if adjidx<=2
    Qexp(adjidx,adjidx) = rate(adjidx)*Qexp(adjidx,adjidx);
else
    Rexp(adjidx-2,adjidx-2) = rate(adjidx)*Rexp(adjidx-2,adjidx-2);
end

losscur = 0;
lossexp = 0;
if strcmp(QRestimator.mode,'fix')
    xf_history = QRestimator.xf_history;
    Pf_history = QRestimator.Pf_history;
    losscur = calloss(xf_history,Pf_history,uhistorys,zhistorys,F,H,B,G,Qcur,Rcur,QRestimator.mode);
    lossexp = calloss(xf_history,Pf_history,uhistorys,zhistorys,F,H,B,G,Qexp,Rexp,QRestimator.mode);
elseif strcmp(QRestimator.mode,'delay')
    delaylen = QRestimator.delaylen;
    xf_historys = QRestimator.xf_historys;
    Pf_historys = QRestimator.Pf_historys;
    for i = 1:QRestimator.winlen
        losscur = losscur + calloss(xf_historys(:,i),Pf_historys(:,:,i),uhistorys(:,i:i+delaylen-1),zhistorys(:,i:i+delaylen-1),F,H,B,G,Qcur,Rcur,QRestimator.mode);
        lossexp = lossexp + calloss(xf_historys(:,i),Pf_historys(:,:,i),uhistorys(:,i:i+delaylen-1),zhistorys(:,i:i+delaylen-1),F,H,B,G,Qexp,Rexp,QRestimator.mode);
    end
else
    display("mode error!!!!");
end

if lossexp < losscur
    Qest = Qexp;
    Rest = Rexp;
else
    Qest = Qcur;
    Rest = Rcur;
    QRestimator.rate(adjidx) = 1/rate(adjidx);
end
QRestimator.adjidx = mod(adjidx,length(rate))+1;

QRestimator = stepupdate(QRestimator);
if strcmp(QRestimator.updatemode,'slidewindow')
    if strcmp(QRestimator.mode,'fix')
        if QRestimator.measurenum == QRestimator.winlen
            [QRestimator.xf_history,QRestimator.Pf_history,~,~,~] = kf_scd(xf_history,Pf_history,F,H,B,uhistorys(:,1),G,zhistorys(:,1),Qest,Rest);
            QRestimator.zhistorys(:,1:end-1) = QRestimator.zhistorys(:,2:end);
            QRestimator.uhistorys(:,1:end-1) = QRestimator.uhistorys(:,2:end);
            QRestimator.measurenum = QRestimator.measurenum-1;
        end
    elseif strcmp(QRestimator.mode,'delay')
        if QRestimator.measurenum == QRestimator.totallen
            QRestimator.xf_historys(:,1:end-1) = QRestimator.xf_historys(:,2:end);
            QRestimator.Pf_historys(:,:,1:end-1) = QRestimator.Pf_historys(:,:,2:end);
            QRestimator.zhistorys(:,1:end-1) = QRestimator.zhistorys(:,2:end);
            QRestimator.uhistorys(:,1:end-1) = QRestimator.uhistorys(:,2:end);
            QRestimator.measurenum = QRestimator.measurenum-1;
        end
    end
elseif strcmp(QRestimator.updatemode,'intermittent')
    if strcmp(QRestimator.mode,'fix')
        QRestimator.measurenum = 0;
    elseif strcmp(QRestimator.mode,'delay')
        num = QRestimator.winlen;
        QRestimator.measurenum = QRestimator.measurenum-num;
        QRestimator.xf_historys(:,1:end-num) = QRestimator.xf_historys(:,num+1:end);
        QRestimator.Pf_historys(:,:,1:end-num) = QRestimator.Pf_historys(:,:,num+1:end);
        QRestimator.zhistorys(:,1:end-num) = QRestimator.zhistorys(:,num+1:end);
        QRestimator.uhistorys(:,1:end-num) = QRestimator.uhistorys(:,num+1:end);
        
    end
else
    display("updatemode error!!!!!!!!!!")
end
end