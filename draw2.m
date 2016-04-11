function [  ] = draw2( A,B,distForGroup,distForGroupT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure(11);
subplot(211);
a=zeros(1,10);
plot(a);
[~,idx] = sort(distForGroup(A,:),'ascend');%竹排排序
NoBtoA=find(idx==B);
truthA=find(idx==A);
[~,idxT] = sort(distForGroupT(A,:),'ascend');%竹排排序
NoBtoAT=find(idxT==B);
truthAT=find(idxT==A);
text(2,0.8,'A ');
text(2,0.6,['第 ',num2str(A),'号组群']);
text(2,0.4,['B在这里排第 ',num2str(NoBtoA),'名','(T',num2str(NoBtoAT),')']);
text(2,0.2,['镜像排第 ',num2str(truthA),'名','(T',num2str(truthAT),')']);
text(4,0.8,['相似度',num2str(distForGroup(A,B))]);
text(4,0.6,['T相似度',num2str(distForGroupT(A,B))]);
[~,idx] = sort(distForGroup(:,B),'ascend');%竹排排序
NoAtoB=find(idx==A);
truthB=find(idx==B);
[~,idxT] = sort(distForGroupT(:,B),'ascend');%竹排排序
NoAtoBT=find(idxT==A);
truthBT=find(idxT==B);
text(6,0.8,'B ');
text(6,0.6,['第 ',num2str(B),'号组群']);
text(6,0.4,['A在这里排第 ',num2str(NoAtoB),'名','(T',num2str(NoAtoBT),')']);
text(6,0.2,['镜像排第 ',num2str(truthB),'名','(T',num2str(truthBT),')']);
subplot(212);
spem=ones(1,size(distForGroup,2));
spem(NoBtoA)=2;
spem(truthA)=3;
plot(spem,'Marker','*');
end

