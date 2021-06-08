function [ Mline,Slinepointt] = line2D( PC,pr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
PC2D=PC(:,1:2);
% plot(PC2D(:,1),PC2D(:,2),'.b','MarkerSize',1);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off;
[n m]=size(PC2D);
[idx dist]=rangesearch(PC2D,PC2D,0.3*pr);
for i=1:n
    density(i)=length(idx{i});
end
numden=5;
iddensity=find(density>numden);
linepoint=PC2D(iddensity,:);
% plot(linepoint(:,1),linepoint(:,2),'.b','MarkerSize',5);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off;
pclinepoint=pointCloud([linepoint zeros(length(linepoint),1)]);
Spclinepoint=pcdownsample(pclinepoint,'gridAverage',pr);
Slinepoint=Spclinepoint.Location;
Slinepointt=Slinepoint(:,1:2);
Slinepoint=Slinepoint(:,1:2);
% plot(Slinepoint(:,1),Slinepoint(:,2),'.b','MarkerSize',2);
% set(gca,'DataAspectRatio',[1 1 1]);
% axis off;
[n m]=size(Slinepoint);
[idx dist]=knnsearch(Slinepoint,Slinepoint,'k',10);
Mvar=[];
for i=1:n
    KNN=Slinepoint(idx(i,:),:);
    [h l]=size(KNN);
    M=(KNN-ones(h,1)*mean(KNN))'*(KNN-ones(h,1)*mean(KNN));
    [U S V]=svd(M);
    Mvar=[Mvar;S(2,2)/S(1,1)];
end
Mline=[];
for i=1:10000
    [mi, Sidx]=min(Mvar);
    [idx1 dist1]=knnsearch(Slinepoint,Slinepoint(Sidx,:),'k',10);
    KNN=Slinepoint(idx1,:);
    SG=KNN;
    Slinepoint(idx1,:)=[];
    Mvar(idx1,:)=[];
    for count=1:10000
        [h0 l0]=size(SG);
        [idx2 dist2]=rangesearch(Slinepoint,KNN,3*pr);
        [h l]=size(KNN);
        idKNNN=[];
        for j=1:h
            idKNNN=[idKNNN idx2{j}];
        end
        idKNNN=unique(idKNNN);
        KNNN=Slinepoint(idKNNN,:);
        [h1 l1]=size(KNNN);
        for k=1:h1
            CSG=[SG;KNNN(k,:)];
            [h2 l2]=size(CSG);
            M1=(CSG-ones(h2,1)*mean(CSG))'*(CSG-ones(h2,1)*mean(CSG));
            [U1 S1 V1]=svd(M1);
            nvector=U1(:,2);
            c=-mean(CSG)* nvector;
            d=abs(KNNN(k,:)*nvector+c);
            if d<0.5*pr            %²ÎÊý
                SG=CSG;
                idcc=find(Slinepoint(:,1)==KNNN(k,1) & Slinepoint(:,2)==KNNN(k,2));
                Slinepoint(idcc,:)=[];
                Mvar(idcc,:)=[];
            end
        end
        [h3 l3]=size(SG);
        if h3==h0
            break;
        end
        KNN=SG(h0+1:h3,:);
    end
    if h3>15
        Mline=[Mline;nvector' c];
%         figure(i);
%         plot(SG(:,1),SG(:,2),'.r','MarkerSize',10);
%         hold on;
%         plot(Slinepoint(:,1),Slinepoint(:,2),'.b','MarkerSize',1);
%         set(gca,'DataAspectRatio',[1 1 1]);
    end
    if mi>0.3*pr
        break;
    end
end
end

