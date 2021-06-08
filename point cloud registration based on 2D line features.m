%基于二维线特征的点云配准
pcloud1=pcread('D:\compile document\matlab\data\Indoor and outdoor dataset\apartment1.ply');
PC1=pcloud1.Location;
pcloud2=pcread('D:\compile document\matlab\data\Indoor and outdoor dataset\apartment2.ply');
PC2=pcloud2.Location;
% plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.b','MarkerSize',1);
% hold on;
% plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.r','MarkerSize',1);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off
pr=0.044106349;       %apartment1  apartment2
% pr=0.022431135;          %apartment4  apartment6 
% pr=0.025664071;          %boardroom0  boardroom3 
% pr=0.16526765;          %castle1  castle2 
% pr=0.12308180;          %City1  City2 
tic;
[ Mline1,Slinepoint1] = line2D( PC1,pr);
% plot(Slinepoint1(:,1),Slinepoint1(:,2),'.b','MarkerSize',1);
% hold on;
[ Mline2,Slinepoint2] = line2D( PC2,pr);
t=toc;
tic;
Mline2=[Mline2;-Mline2];
% plot(Slinepoint2(:,1),Slinepoint2(:,2),'.b','MarkerSize',1);
% set(gca,'DataAspectRatio',[1 1 1]);
[n1 m1]=size(Mline1);
Mcostheta1=[];
Mindex1=[];
for i=1:(n1-1)
    for j=i+1:n1
        costheta=Mline1(i,1:2)*Mline1(j,1:2)';
        index=[i j];
        if abs(costheta)<cos(10*pi/180)
        Mcostheta1=[Mcostheta1;costheta];
        Mindex1=[Mindex1;index];
        end
    end
end
[n2 m2]=size(Mline2);
Mcostheta2=[];
Mindex2=[];
for i=1:(n2-1)
    for j=i+1:n2
        costheta=Mline2(i,1:2)*Mline2(j,1:2)';
        index=[i j];
        if abs(costheta)<cos(10*pi/180)
            Mcostheta2=[Mcostheta2;costheta];
            Mindex2=[Mindex2;index];
        end
    end
end
[n3 m3]=size(Mcostheta1);
[n4 m4]=size(Mcostheta2);
overlap00=0;
for i=1:n3
    for j=1:n4
        if abs(Mcostheta1(i)-Mcostheta2(j))<0.2*pr    %参数
            match=[Mindex1(i,:);Mindex2(j,:)];
            n11=Mline1(match(1,1),1:2);
            c11=Mline1(match(1,1),3);
            n12=Mline1(match(1,2),1:2);
            c12=Mline1(match(1,2),3);
            n21=Mline2(match(2,1),1:2);
            c21=Mline2(match(2,1),3);
            n22=Mline2(match(2,2),1:2);
            c22=Mline2(match(2,2),3);
            [U1 S1 V1]=svd([n11' n12']*[n21' n22']');
            r1=V1*U1';
            [U5 S5 V5]=svd([n11' n12']*[n22' n21']');
            r5=V5*U5';
            t1=inv([n21;n22])*[c11-c21;c12-c22];
            t5=inv([n22;n21])*[c11-c21;c12-c22];
            Mr=[r1 r5];
            Mt=[t1 t5];            
            overlap0=0;
            for k=1:2
                r=Mr(:,2*k-1:2*k);
                t=Mt(:,k);
                tMline1=[Mline1(:,1:2)*r' Mline1(:,3)-Mline1(:,1:2)*r'*t];
                [idx dist]=knnsearch(Mline2,tMline1,'k',1);
                idoverlap=find(dist<3*pr);           %参数
                overlap=length(idoverlap)/min(n1,0.5*n2);
                if overlap>overlap0
                    tt=[r t];
                    overlap0=overlap;
                end
            end
            if overlap0>overlap00
            ttt=tt;
            overlap00=overlap0;
            end
        end
    end
end
r=ttt(1:2,1:2);
t=ttt(1:2,3);
tSlinepoint1=Slinepoint1*r'+ones(length(Slinepoint1),1)*t';
% plot(tSlinepoint1(:,1),tSlinepoint1(:,2),'.b','MarkerSize',1);
% hold on;
% plot(Slinepoint2(:,1),Slinepoint2(:,2),'.r','MarkerSize',1);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off;
[idx1 dist1]=knnsearch(Slinepoint2,tSlinepoint1,'k',1);
idoverlap1=find(dist1<pr);
overlappoint=tSlinepoint1(idoverlap1,:);
% plot(overlappoint(:,1),overlappoint(:,2),'.r','MarkerSize',1);
% set(gca,'DataAspectRatio',[1 1 1]);
Mtz=[];
for j=1:50
    Sidx=randperm(length(overlappoint),1);
    [idxx1 distt1]=rangesearch(PC1(:,1:2),Slinepoint1(idoverlap1(Sidx),:),2*pr);
    KNN1=PC1(idxx1{1},:);
%     Sop1=Slinepoint1(idoverlap1(Sidx),:);
%     subplot(1,2,1);
%     plot(Sop1(:,1),Sop1(:,2),'.b','MarkerSize',10);
%     hold on;
%     plot(Slinepoint1(:,1),Slinepoint1(:,2),'.r','MarkerSize',1);
%     set(gca,'DataAspectRatio',[1 1 1]);
%     axis off
%     subplot(1,2,2);
%     plot3(KNN1(:,1),KNN1(:,2),KNN1(:,3),'.b','MarkerSize',5);
%     hold on;
%     plot3(PC1(:,1),PC1(:,2),PC1(:,3),'.r','MarkerSize',1);
%     set(gca,'DataAspectRatio',[1 1 1]);
%     axis off
    [idxx2 distt2]=rangesearch(PC2(:,1:2),overlappoint(Sidx,:),2*pr);
    KNN2=PC2(idxx2{1},:);
%     Sop2=overlappoint(Sidx,:);
%     subplot(1,2,1);
%     plot(Sop2(:,1),Sop2(:,2),'.b','MarkerSize',10);
%     hold on;
%     plot(Slinepoint2(:,1),Slinepoint2(:,2),'.r','MarkerSize',1);
%     set(gca,'DataAspectRatio',[1 1 1]);
%     axis off
%     subplot(1,2,2);
%     plot3(KNN2(:,1),KNN2(:,2),KNN2(:,3),'.b','MarkerSize',5);
%     hold on;
%     plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.r','MarkerSize',1);
%     set(gca,'DataAspectRatio',[1 1 1]);
%     axis off
    zmin1=min(KNN1(:,3));
    zmin2=min(KNN2(:,3));
    tz=zmin2-zmin1;
    Mtz=[Mtz;tz 0];
end
MStz=[];
for i=1:1000
    seed=Mtz(1,:);
    Mtz(1,:)=[];
    for j=1:1000
        [idxtz disttz]=knnsearch(Mtz,seed,'k',1);
        [mi id]=min(disttz);
        if mi<0.2*pr         %参数
            seed=[seed;Mtz(idxtz(id),:)];
            Mtz(idxtz(id),:)=[];
        else
            break;
        end
    end
    [h l]=size(seed);
    Stz=[mean(seed(:,1)) h];
    MStz=[MStz;Stz];
    if isempty(Mtz)
        break
    end
end
[ma idStz]=max(MStz(:,2));
tz=MStz(idStz,1);
R=[r zeros(2,1);zeros(1,2) 1];
T=[t;tz];
TT=[R T;zeros(1,3) 1];
PC1t=PC1*TT(1:3,1:3)'+ones(length(PC1),1)*TT(1:3,4)';
t=toc;
% plot3(PC1t(:,1),PC1t(:,2),PC1t(:,3),'.b','MarkerSize',1);
% hold on;
% plot3(PC2(:,1),PC2(:,2),PC2(:,3),'.r','MarkerSize',1);
% set(gca,'DataAspectRatio',[1 1 1]);
axis off;
tic;
[Tg,PC1t,RMS] = partialoverlapICP3(PC2,PC1t,pr);
t=toc;
TT=Tg*TT;
Rzhen=TT(1:3,1:3);
Tzhen=TT(1:3,4);
errorR=real(acos((trace(Rzhen*inv(R))-1)/2)*(180/pi));
errorTh=norm(Tzhen(1:2)-T(1:2));
errorTv=norm(Tzhen(3)-T(3));
zma=max(PC2(:,3));
zmi=min(PC2(:,3));
interval=(zma-zmi)/5;
idPC1t=find(PC1t(:,3)<zma-interval);
PC1tR=PC1t(idPC1t,:);
idPC2=find(PC2(:,3)<zma-interval);
PC2R=PC2(idPC2,:);
plot3(PC1tR(:,1),PC1tR(:,2),PC1tR(:,3),'.b','MarkerSize',1);
hold on;
plot3(PC2R(:,1),PC2R(:,2),PC2R(:,3),'.r','MarkerSize',1);
set(gca,'DataAspectRatio',[1 1 1]);
axis off