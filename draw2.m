function [  ] = draw2( A,B,distForGroup,distForGroupT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure(11);
subplot(211);
a=zeros(1,10);
plot(a);
[~,idx] = sort(distForGroup(A,:),'ascend');%��������
NoBtoA=find(idx==B);
truthA=find(idx==A);
[~,idxT] = sort(distForGroupT(A,:),'ascend');%��������
NoBtoAT=find(idxT==B);
truthAT=find(idxT==A);
text(2,0.8,'A ');
text(2,0.6,['�� ',num2str(A),'����Ⱥ']);
text(2,0.4,['B�������ŵ� ',num2str(NoBtoA),'��','(T',num2str(NoBtoAT),')']);
text(2,0.2,['�����ŵ� ',num2str(truthA),'��','(T',num2str(truthAT),')']);
text(4,0.8,['���ƶ�',num2str(distForGroup(A,B))]);
text(4,0.6,['T���ƶ�',num2str(distForGroupT(A,B))]);
[~,idx] = sort(distForGroup(:,B),'ascend');%��������
NoAtoB=find(idx==A);
truthB=find(idx==B);
[~,idxT] = sort(distForGroupT(:,B),'ascend');%��������
NoAtoBT=find(idxT==A);
truthBT=find(idxT==B);
text(6,0.8,'B ');
text(6,0.6,['�� ',num2str(B),'����Ⱥ']);
text(6,0.4,['A�������ŵ� ',num2str(NoAtoB),'��','(T',num2str(NoAtoBT),')']);
text(6,0.2,['�����ŵ� ',num2str(truthB),'��','(T',num2str(truthBT),')']);
subplot(212);
spem=ones(1,size(distForGroup,2));
spem(NoBtoA)=2;
spem(truthA)=3;
plot(spem,'Marker','*');
end

