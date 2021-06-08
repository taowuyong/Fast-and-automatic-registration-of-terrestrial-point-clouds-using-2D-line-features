function [Tg,PC2t,RMS,overlap,iteration] = partialoverlapICP3(PC1,PC2,pr)
%UNTITLED 此处显示有关此函数的摘要
Rg=eye(3);
tg=zeros(1,3);
[n,m]=size(PC2);
Mdl = KDTreeSearcher(PC1);
for i=1:500
    [idx,dist]=rangesearch(Mdl,PC2,3*pr);
    id2=[];
    id1=[];
    for j=1:n
        if ~isempty(idx{j})
            id2=[id2 j];
            idd1=idx{j};
            id1=[id1 idd1(1)];
        end
    end
    A=PC2(id2,:);
    Y=PC1(id1,:);
    [n1,m1]=size(A);
    uA=[mean(A(:,1)) mean(A(:,2)) mean(A(:,3))];
    uY=[mean(Y(:,1)) mean(Y(:,2)) mean(Y(:,3))];
    H=zeros(3);
    for j=1:n1
        H=H+(A(j,:)-uA)'*(Y(j,:)-uY);
    end
    [U S V]=svd(H);
    D=diag([1 1 det(U*V')]);
    R=V*D*U';
    t=uY-uA*R';
    PC2t=PC2*R'+ones(n,1)*t;
    T=[R t';zeros(1,3) 1];
    if norm(T-eye(4,4))<0.0001
           break
    end
    PC2=PC2t;
    Rg=R*Rg;   
    tg=tg*R'+t;
end
iteration=i;
for j=1:n1
    dd(j)=norm(Y(j,:)-A(j,:));
end
RMS=norm(dd)/sqrt(n1);
Tg=[R*Rg (tg*R'+t)';zeros(1,3) 1];
[n2,m2]=size(PC1);
overlap=length(A)/min(n,n2);
end
