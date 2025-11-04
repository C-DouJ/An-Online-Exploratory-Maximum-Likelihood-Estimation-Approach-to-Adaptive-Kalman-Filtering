%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% 蒙特卡洛仿真程序 %%%%%%%%%
%%%%%%%%   目标跟踪模型  %%%%%%%%%
clear all;
close all;
clc;
%profile on
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));
% randn('state',sum(100*clock));
format long;
% tic;
%%%%%Model parameters
nxp=100;                %%MC次数
nx=4;
nz=2;
T=1;
q1=100;               %%调整名义Q和真实Q之间的倍数关系
q2=100;
r1=0.01;               %%调整名义Q和真实Q之间的倍数关系
r2=0.01;
ts=3000;              %%单次仿真总时间

%%%%%State space model
F=[eye(2) T*eye(2);zeros(2) eye(2)];
H=[eye(2) zeros(2)];
B=zeros(4,1);
u=0;
G=[0.5*T^2,0;0,0.5*T^2;T,0;0,T];

%True noise covariance matrices


Q=zeros(2,2,ts);
Q(1,1,1) = 10;
Q(2,2,1) = 10;
for i=1:ts-1
    Q(1,1,i+1) = 0*max(Q(1,1,i)+0.1*randn,1)+10;
    Q(2,2,i+1) = 0*max(Q(2,2,i)+0.1*randn,1)+10;
end
R=zeros(2,2,ts);
R(1,1,1) = 1;
R(2,2,1) = 1;
for i=1:ts-1
    R(1,1,i+1) = 1+0*max(R(1,1,i)+0.01*randn,0.1);
    R(2,2,i+1) = 1+0*max(R(2,2,i)+0.01*randn,0.1);
end
% load vari_simA_data.mat

%%%%Nominal noise covariance matrices
Q0=diag([Q(1,1,1)/q1;Q(2,2,1)/q2]);   %名义Q
R0=diag([R(1,1,1)/r1;R(2,2,1)/r2]);   %名义Q

saveale = zeros(ts, 4);
pos_rmse_tkf = zeros(ts,1);
pos_rmse_nkf = zeros(ts,1);
pos_rmse_SCDfixim = zeros(ts,1);
pos_rmse_SCDfixsw = zeros(ts,1);
pos_rmse_SCDdelayim = zeros(ts,1);
pos_rmse_SCDdelaysw = zeros(ts,1);

vel_rmse_tkf = zeros(ts,1);
vel_rmse_nkf = zeros(ts,1);
vel_rmse_SCDfixim = zeros(ts,1);
vel_rmse_SCDfixsw = zeros(ts,1);
vel_rmse_SCDdelayim = zeros(ts,1);
vel_rmse_SCDdelaysw = zeros(ts,1);

t_SCDfixim=0;
t_SCDfixsw=0;
t_SCDdelayim=0;
t_SCDdelaysw=0;

Q1ale = zeros(ts, 4);
Q2ale = zeros(ts, 4);
R1ale = zeros(ts, 4);
R2ale = zeros(ts, 4);
for expt=1:nxp
    
    fprintf('MC Run in Process=%d\n',expt);
    
    %%%%变量初始化
    save_qA1 = zeros(ts,4);
    save_qA2 = zeros(ts,4);
    save_rA1 = zeros(ts,4);
    save_rA2 = zeros(ts,4);
    
    %%%%%Set the system initial value%%%%%%%%%%%
    x=[0;0;10;10];                     %%%True initial state value
    P=diag([1 1 1 1]);             %%%Initial estimation error covariance matrix
    Skk=utchol(P);                         %%%Square-root of initial estimation error covariance matrix
    
    
    %%%%Initial state estimate of standard KF    (Kalman filter)
    xf=x+Skk*randn(nx,1);
    Pf=P;
    
    xf_t=xf;
    Pf_t=P;
    xf_n=xf;
    Pf_n=P;
    xf_SCDfixim=xf;
    Pf_SCDfixim=P;
    xf_SCDdelayim=xf;
    Pf_SCDdelayim=P;
    xf_SCDfixsw=xf;
    Pf_SCDfixsw=P;
    xf_SCDdelaysw=xf;
    Pf_SCDdelaysw=P;
    
    Q_SCDfixim = Q0;
    Q_SCDdelayim = Q0;
    Q_SCDfixsw = Q0;
    Q_SCDdelaysw = Q0;
    R_SCDfixim = R0;
    R_SCDdelayim = R0;
    R_SCDfixsw = R0;
    R_SCDdelaysw = R0;
    %%%%%%%%%%%%%%%%%%%%%%
    QRestimatorfixim = createQRestimator('fix','intermittent',5,0,xf_SCDfixim,Pf_SCDfixim); % PICD
    QRestimatordelayim = createQRestimator('delay','intermittent',5,3,xf_SCDdelayim,Pf_SCDdelayim); % ICD
    QRestimatorfixsw = createQRestimator('fix','slidewindow',10,0,xf_SCDfixsw,Pf_SCDfixsw); % PSCD
    QRestimatordelaysw = createQRestimator('delay','slidewindow',5,3,xf_SCDdelaysw,Pf_SCDdelaysw); %SCD
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%
    xA=[];
    ZA=[];
    zA=[];
    xkk_A=[];
    Pkk_A=[];
    for t=1:ts
        
        
        %%%%Square-root of noise covariance matrices
        SQ=utchol(Q(:,:,t));
        SR=utchol(R(:,:,t));
        
        %%%%Simulate true state and measurement
        x=F*x+B*u+G*SQ*randn(2,1);
        z=H*x+SR*randn(nz,1);
        
        
        %%%Save data
        xA=[xA x];
        ZA=[ZA z];
        
        
        %%%%%%%%%%%%%%%%%%%%%%%             Run TKF                  %%%%%%%%%%%%%%%%%%%%%%%
        [xf_t,Pf_t,~]=kf_std(xf_t,Pf_t,F,H,B,u,G,z,Q(:,:,t),R(:,:,t));
        
        %%%%%%%%%%%%%%%%%%%%%%%             Run NKF                  %%%%%%%%%%%%%%%%%%%%%%%
        [xf_n,Pf_n,~]=kf_std(xf_n,Pf_n,F,H,B,u,G,z,Q0,R0);
        
        %%%%%%%%%%%%%%%%%%%%%%%             Run SCDfixsw (PSCD)                %%%%%%%%%%%%%%%%%%%%%%%
        tic;
        QRestimatorfixsw.zhistorys(:,QRestimatorfixsw.measurenum+1) = z;
        QRestimatorfixsw.uhistorys(:,QRestimatorfixsw.measurenum+1) = u;
        if strcmp(QRestimatorfixsw.mode,'delay')
            QRestimatorfixsw.xf_historys(:,QRestimatorfixsw.measurenum+1) = xf_SCDfixsw;
            QRestimatorfixsw.Pf_historys(:,:,QRestimatorfixsw.measurenum+1) = Pf_SCDfixsw;
        end
        QRestimatorfixsw.measurenum = QRestimatorfixsw.measurenum+1;
        if QRestimatorfixsw.measurenum > QRestimatorfixsw.totallen, display("error!!!!!!!!!!");end
        if QRestimatorfixsw.totallen == QRestimatorfixsw.measurenum
            [Qest,Rest,QRestimatorfixsw] = estimateQR(QRestimatorfixsw,F,H,B,G,Q_SCDfixsw,R_SCDfixsw);
            Q_SCDfixsw = Qest;
            R_SCDfixsw = Rest;
        end
        [xf_SCDfixsw,Pf_SCDfixsw,~,~,~]=kf_scd(xf_SCDfixsw,Pf_SCDfixsw,F,H,B,u,G,z,Q_SCDfixsw,R_SCDfixsw);
        t_SCDfixsw = t_SCDfixsw + toc/nxp/ts;
        
        %%%%%%%%%%%%%%%%%%%%%%%             Run SCDdelaysw (SCD)                 %%%%%%%%%%%%%%%%%%%%%%%
        tic;
        QRestimatordelaysw.zhistorys(:,QRestimatordelaysw.measurenum+1) = z;
        QRestimatordelaysw.uhistorys(:,QRestimatordelaysw.measurenum+1) = u;
        if strcmp(QRestimatordelaysw.mode,'delay')
            QRestimatordelaysw.xf_historys(:,QRestimatordelaysw.measurenum+1) = xf_SCDdelaysw;
            QRestimatordelaysw.Pf_historys(:,:,QRestimatordelaysw.measurenum+1) = Pf_SCDdelaysw;
        end
        QRestimatordelaysw.measurenum = QRestimatordelaysw.measurenum+1;
        if QRestimatordelaysw.measurenum > QRestimatordelaysw.totallen, display("error!!!!!!!!!!");end
        if QRestimatordelaysw.totallen == QRestimatordelaysw.measurenum
            [Qest,Rest,QRestimatordelaysw] = estimateQR(QRestimatordelaysw,F,H,B,G,Q_SCDdelaysw,R_SCDdelaysw);
            Q_SCDdelaysw = Qest;
            R_SCDdelaysw = Rest;
        end
        [xf_SCDdelaysw,Pf_SCDdelaysw,~,~,~]=kf_scd(xf_SCDdelaysw,Pf_SCDdelaysw,F,H,B,u,G,z,Q_SCDdelaysw,R_SCDdelaysw);
        t_SCDdelaysw = t_SCDdelaysw + toc/nxp/ts;
        
        %%%%%%%%%%%%%%%%%%%%%%%             Run SCDfixim (PICD)                %%%%%%%%%%%%%%%%%%%%%%%
        tic;
        QRestimatorfixim.zhistorys(:,QRestimatorfixim.measurenum+1) = z;
        QRestimatorfixim.uhistorys(:,QRestimatorfixim.measurenum+1) = u;
        if strcmp(QRestimatorfixim.mode,'delay')
            QRestimatorfixim.xf_historys(:,QRestimatorfixim.measurenum+1) = xf_SCDfixim;
            QRestimatorfixim.Pf_historys(:,:,QRestimatorfixim.measurenum+1) = Pf_SCDfixim;
        elseif strcmp(QRestimatorfixim.mode,'fix') && strcmp(QRestimatorfixim.updatemode,'intermittent')
            if QRestimatorfixim.measurenum == 0
                QRestimatorfixim.xf_history = xf_SCDfixim;
                QRestimatorfixim.Pf_history = Pf_SCDfixim;
            end
        end
        QRestimatorfixim.measurenum = QRestimatorfixim.measurenum+1;
        if QRestimatorfixim.measurenum > QRestimatorfixim.totallen, display("error!!!!!!!!!!");end
        if QRestimatorfixim.totallen == QRestimatorfixim.measurenum
            [Qest,Rest,QRestimatorfixim] = estimateQR(QRestimatorfixim,F,H,B,G,Q_SCDfixim,R_SCDfixim);
            Q_SCDfixim = Qest;
            R_SCDfixim = Rest;
        end
        [xf_SCDfixim,Pf_SCDfixim,~,~,~]=kf_scd(xf_SCDfixim,Pf_SCDfixim,F,H,B,u,G,z,Q_SCDfixim,R_SCDfixim);
        t_SCDfixim = t_SCDfixim + toc/nxp/ts;
        
        %%%%%%%%%%%%%%%%%%%%%%%             Run SCDdelayim (ICD)                %%%%%%%%%%%%%%%%%%%%%%%
        tic;
        QRestimatordelayim.zhistorys(:,QRestimatordelayim.measurenum+1) = z;
        QRestimatordelayim.uhistorys(:,QRestimatordelayim.measurenum+1) = u;
        if strcmp(QRestimatordelayim.mode,'delay')
            QRestimatordelayim.xf_historys(:,QRestimatordelayim.measurenum+1) = xf_SCDdelayim;
            QRestimatordelayim.Pf_historys(:,:,QRestimatordelayim.measurenum+1) = Pf_SCDdelayim;
        elseif strcmp(QRestimatorfixim.mode,'fix') && strcmp(QRestimatordelayim.updatemode,'intermittent')
            if QRestimatordelayim.measurenum == 0
                QRestimatordelayim.xf_history = xf_SCDdelayim;
                QRestimatordelayim.Pf_history = Pf_SCDdelayim;
            end
        end
        QRestimatordelayim.measurenum = QRestimatordelayim.measurenum+1;
        if QRestimatordelayim.measurenum > QRestimatordelayim.totallen, display("error!!!!!!!!!!");end
        if QRestimatordelayim.totallen == QRestimatordelayim.measurenum
            [Qest,Rest,QRestimatordelayim] = estimateQR(QRestimatordelayim,F,H,B,G,Q_SCDdelayim,R_SCDdelayim);
            Q_SCDdelayim = Qest;
            R_SCDdelayim = Rest;
        end
        [xf_SCDdelayim,Pf_SCDdelayim,~,~,~]=kf_scd(xf_SCDdelayim,Pf_SCDdelayim,F,H,B,u,G,z,Q_SCDdelayim,R_SCDdelayim);
        t_SCDdelayim = t_SCDdelayim + toc/nxp/ts;
        
        
        save_qA1(t,1) = Q_SCDfixim(1,1); % PICD
        save_qA2(t,1) = Q_SCDfixim(2,2);
        save_rA1(t,1) = R_SCDfixim(1,1);
        save_rA2(t,1) = R_SCDfixim(2,2);
        
        save_qA1(t,2) = Q_SCDfixsw(1,1); % PSCD
        save_qA2(t,2) = Q_SCDfixsw(2,2);
        save_rA1(t,2) = R_SCDfixsw(1,1);
        save_rA2(t,2) = R_SCDfixsw(2,2);
        
        save_qA1(t,3) = Q_SCDdelayim(1,1); % ICD
        save_qA2(t,3) = Q_SCDdelayim(2,2);
        save_rA1(t,3) = R_SCDdelayim(1,1);
        save_rA2(t,3) = R_SCDdelayim(2,2);
        
        save_qA1(t,4) = Q_SCDdelaysw(1,1); % SCD
        save_qA2(t,4) = Q_SCDdelaysw(2,2);
        save_rA1(t,4) = R_SCDdelaysw(1,1);
        save_rA2(t,4) = R_SCDdelaysw(2,2);
        
        %%%%计算估计误差
        x_err_tkf(t,:)=(xf_t-x)';
        x_err_nkf(t,:)=(xf_n-x)';
        x_err_SCDfixsw(t,:)=(xf_SCDfixsw-x)';
        x_err_SCDdelaysw(t,:)=(xf_SCDdelaysw-x)';
        x_err_SCDfixim(t,:)=(xf_SCDfixim-x)';
        x_err_SCDdelayim(t,:)=(xf_SCDdelayim-x)';
        
        %%%%% 计算位置误差
        pos_rmse_tkf(t,1) = sqrt(pos_rmse_tkf(t,1)^2 + (x_err_tkf(t,1)^2 + x_err_tkf(t,2)^2)/nxp);
        pos_rmse_nkf(t,1) = sqrt(pos_rmse_nkf(t,1)^2 + (x_err_nkf(t,1)^2 + x_err_nkf(t,2)^2)/nxp);
        pos_rmse_SCDfixim(t,1) = sqrt(pos_rmse_SCDfixim(t,1)^2 + (x_err_SCDfixim(t,1)^2 + x_err_SCDfixim(t,2)^2)/nxp);
        pos_rmse_SCDfixsw(t,1) = sqrt(pos_rmse_SCDfixsw(t,1)^2 + (x_err_SCDfixsw(t,1)^2 + x_err_SCDfixsw(t,2)^2)/nxp);
        pos_rmse_SCDdelayim(t,1) = sqrt(pos_rmse_SCDdelayim(t,1)^2 + (x_err_SCDdelayim(t,1)^2 + x_err_SCDdelayim(t,2)^2)/nxp);
        pos_rmse_SCDdelaysw(t,1) = sqrt(pos_rmse_SCDdelaysw(t,1)^2 + (x_err_SCDdelaysw(t,1)^2 + x_err_SCDdelaysw(t,2)^2)/nxp);
        
        pos_RMSE = [pos_rmse_tkf pos_rmse_nkf pos_rmse_SCDfixim pos_rmse_SCDfixsw...
            pos_rmse_SCDdelayim pos_rmse_SCDdelaysw];
        
        pos_armse_tkf = rmse(pos_rmse_tkf,0);
        pos_armse_nkf = rmse(pos_rmse_nkf,0);
        pos_armse_SCDfixim = rmse(pos_rmse_SCDfixim,0);
        pos_armse_SCDfixsw = rmse(pos_rmse_SCDfixsw,0);
        pos_armse_SCDdelayim = rmse(pos_rmse_SCDdelayim,0);
        pos_armse_SCDdelaysw = rmse(pos_rmse_SCDdelaysw,0);
        
        pos_ARMSE = [pos_armse_tkf pos_armse_nkf pos_armse_SCDfixim pos_armse_SCDfixsw...
            pos_armse_SCDdelayim pos_armse_SCDdelaysw]; % PICD PSCD ICD SCD
        
        %%%%% 计算速度误差
        vel_rmse_tkf(t,1) = sqrt(vel_rmse_tkf(t,1)^2 + (x_err_tkf(t,3)^2 + x_err_tkf(t,4)^2)/nxp);
        vel_rmse_nkf(t,1) = sqrt(vel_rmse_nkf(t,1)^2 + (x_err_nkf(t,3)^2 + x_err_nkf(t,4)^2)/nxp);
        vel_rmse_SCDfixim(t,1) = sqrt(vel_rmse_SCDfixim(t,1)^2 + (x_err_SCDfixim(t,3)^2 + x_err_SCDfixim(t,4)^2)/nxp);
        vel_rmse_SCDfixsw(t,1) = sqrt(vel_rmse_SCDfixsw(t,1)^2 + (x_err_SCDfixsw(t,3)^2 + x_err_SCDfixsw(t,4)^2)/nxp);
        vel_rmse_SCDdelayim(t,1) = sqrt(vel_rmse_SCDdelayim(t,1)^2 + (x_err_SCDdelayim(t,3)^2 + x_err_SCDdelayim(t,4)^2)/nxp);
        vel_rmse_SCDdelaysw(t,1) = sqrt(vel_rmse_SCDdelaysw(t,1)^2 + (x_err_SCDdelaysw(t,3)^2 + x_err_SCDdelaysw(t,4)^2)/nxp);
        vel_RMSE = [vel_rmse_tkf vel_rmse_nkf vel_rmse_SCDfixim vel_rmse_SCDfixsw...
            vel_rmse_SCDdelayim vel_rmse_SCDdelaysw];
        
        vel_armse_tkf = rmse(vel_rmse_tkf,0);
        vel_armse_nkf = rmse(vel_rmse_nkf,0);
        vel_armse_SCDfixim = rmse(vel_rmse_SCDfixim,0);
        vel_armse_SCDfixsw = rmse(vel_rmse_SCDfixsw,0);
        vel_armse_SCDdelayim = rmse(vel_rmse_SCDdelayim,0);
        vel_armse_SCDdelaysw = rmse(vel_rmse_SCDdelaysw,0);
        vel_ARMSE = [vel_armse_tkf vel_armse_nkf vel_armse_SCDfixim vel_armse_SCDfixsw...
            vel_armse_SCDdelayim vel_armse_SCDdelaysw];
        
        Q1ale(t,:) = Q1ale(t,:) + abs(log10(save_qA1(t,:)/Q(1,1,t))) / nxp;
        Q2ale(t,:) = Q2ale(t,:) + abs(log10(save_qA1(t,:)/Q(2,2,t))) / nxp;
        R1ale(t,:) = R1ale(t,:) + abs(log10(save_rA1(t,:)/R(1,1,t))) / nxp;
        R2ale(t,:) = R2ale(t,:) + abs(log10(save_rA1(t,:)/R(2,2,t))) / nxp;
    end
    %profile viewer
    
    
    saveale = saveale + 0.25 * (Q1ale + Q2ale + R1ale + R2ale) / nxp; 
    aale = mean(saveale);
end
% toc;

% plotRes2
t_run = ts*[t_SCDfixim t_SCDfixsw t_SCDdelayim t_SCDdelaysw];

plot__result
