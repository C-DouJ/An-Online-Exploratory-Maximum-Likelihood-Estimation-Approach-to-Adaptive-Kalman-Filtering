function loss = calloss(xf_history,Pf_history,uhistorys,zhistorys,F,H,B,G,Q,R,mode)
    loss = 0;
    xf = xf_history;
    Pf = Pf_history;
    if strcmp(mode,'fix')
        for t = 1:size(zhistorys,2)
            u = uhistorys(:,t);
            z = zhistorys(:,t);
            [xf,Pf,~,rk,Pzz]=kf_scd(xf,Pf,F,H,B,u,G,z,Q,R);
            loss = loss + rk'*inv(Pzz)*rk + log(det(Pzz));
        end
    elseif strcmp(mode,'delay')
        for t = 1:size(zhistorys,2)
            u = uhistorys(:,t);
            z = zhistorys(:,t);
            [xf,Pf,~,rk,Pzz]=kf_scd(xf,Pf,F,H,B,u,G,z,Q,R);
        end
        loss = rk'*inv(Pzz)*rk + log(det(Pzz));
    else
        display("mode error!!!!");
    end

end
