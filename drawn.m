function [h,w] = drawn(imgA,imgB,relationsPa,relationsPe,peopleNumA2,peopleNumB2,option,wBlank,wBlank2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%imA包括pic,no,person,feature,edgefeat。person包括data{n,2}，num，patch{n,pN,2},patchfeature{n,6}
global haxis;
pN=imgA.person.patchNum;
locA=cell(imgA.person.num*pN,2);
locB=cell(imgB.person.num*pN,2);
counter=1;
color='rgbcmykw';
for i=1:imgA.person.num
    for j=1:pN
        locA{counter,1}=imgA.person.patch{i,j,1};
        locA{counter,2}=imgA.person.patch{i,j,2};
        counter=counter+1;
    end
end
counter=1;
for i=1:imgB.person.num
    for j=1:pN
        locB{counter,1}=imgB.person.patch{i,j,1};
        locB{counter,2}=imgB.person.patch{i,j,2};
        counter=counter+1;
    end
end
imb=imgB.pic;
ima=imgA.pic;


[ah,aw,~]=size(ima);
[bh,bw,~]=size(imb);
h=max(ah,bh);
w=aw+bw;
imO=imread('blank.png');
imO=imresize(imO,[h,w]);
imO(1:ah,1:aw,:)=ima;
imO(1:bh,(1+aw):(aw+bw),:)=imb;
figure(10);
imshow(imO);
%%画heatmap
heat=zeros(ah,aw);
for i=1:imgA.person.num*pN
    a=locA{i,1}(2);b=a+locA{i,2}(2)-1;
    c=locA{i,1}(1);d=c+locA{i,2}(1)-1;
    heat(a:b,c:d)=repmat(imgA.wN(i),b-a+1,d-c+1);
end
cc=imgA.person.num;
if option==0
    for i=1:cc
        for j=1:pN
            x=locA{(i-1)*pN+j,1}(1)+round(locA{(i-1)*pN+j,2}(1)/2);
            y=locA{(i-1)*pN+j,1}(2)+round(locA{(i-1)*pN+j,2}(2)/2);
            for k=i+1:cc
                for b=1:pN
                    heat=myLine(heat,[x,y],[locA{(k-1)*pN+b,1}(1)+round(locA{(k-1)*pN+b,2}(1)/2),...
                        locA{(k-1)*pN+b,1}(2)+round(locA{(k-1)*pN+b,2}(2)/2)],...
                        imgA.wE( ((i-1)*cc-i+1-(i-1)*(i-2)/2)*pN^2+(k-i-1)*pN^2+(j-1)*pN+b ));
                end
            end
        end
    end
end
hold on;
colormap('hot');
haxis=imagesc(heat);
alpha(haxis,0.5);
if option~=3
    colorbar; 
end
hold off;
%%结束
if option==1||option==3||option==4
    according=relationsPa{imgA.no,imgB.no};
    for k=1:size(according,1)%k代表A中的序号
        co=color(ceil(k/pN));
        figure(10);
        if option==3||option==4
            tmpWValue=sum(imgA.wN(k,:))/size(imgA.wN(k,:),2);
            text(round(locA{k,1}(1)+locA{k,2}(1)/2),...
                    round(locA{k,1}(2)+locA{k,2}(2)/2),...
                    num2str(tmpWValue));
            if wBlank(k)==1
                rectangle('Position',[locA{k,1}(1)+5,locA{k,1}(2)+5,locA{k,2}(1)-10,locA{k,2}(2)-10],'Edgecolor','g','LineWidth',2);
            end
        end
        for j=1:size(according,2)%j代表B中的序号
            if according(k,j)==1
                rectangle('Position',[locA{k,1}(1),locA{k,1}(2),locA{k,2}(1),locA{k,2}(2)]);
                line([locA{k,1}(1)+round(locA{k,2}(1)/2),locB{j,1}(1)+round(locB{j,2}(1)/2)+aw],[locA{k,1}(2)+round(locA{k,2}(2)/2),locB{j,1}(2)+round(locB{j,2}(2)/2)],'Color',co);
                %[~,idx] = sort(dist(peopleNumA2(imgA.no)+k,:),'ascend');%竹排排序
                rectangle('Position',[locB{j,1}(1)+aw,locB{j,1}(2),locB{j,2}(1),locB{j,2}(2)]);
                %rank=find(idx==(peopleNumB2(bID)+j));
                %text(round((locA(k,1)+locB(j,1)+aw)/2)-5,round((locA(k,2)+locB(j,2))/2)-5,[num2str(dist(peopleNumA2(aID)+k,peopleNumB2(bID)+j)),'(',num2str(rank),')'],'color','green');
            end
        end
    end
elseif option==2
    according=relationsPe{imgA.no,imgB.no};
    for k=1:size(according,1)%k代表A中的序号
        for j=1:size(according,2)%j代表B中的序号
            if according(k,j)==1
                figure(10);
                line([imgA.person.data{k,1}(1)+round(imgA.person.data{k,2}(1)/2),...
                    imgB.person.data{j,1}(1)+round(imgB.person.data{j,2}(1)/2)+aw],...
                    [imgA.person.data{k,1}(2)+round(imgA.person.data{k,2}(2)/2),...
                    imgB.person.data{j,1}(2)+round(imgB.person.data{j,2}(2)/2)]...
                    ,'Color','g');%前面两个是横坐标，后面两个是纵坐标
                
            end
        end
    end
end
if option==4
    for iForwBlank2=1:size(wBlank2,2)
        if wBlank2(iForwBlank2)==1
            tmp=imgA.wE(iForwBlank2);
            a=imgA.edgefeat(1,iForwBlank2);b=imgA.edgefeat(2,iForwBlank2);
            line([locA{a,1}(1)+round(locA{a,2}(1)/2),locA{b,1}(1)+round(locA{b,2}(1)/2)],[locA{a,1}(2)+round(locA{a,2}(2)/2),locA{b,1}(2)+round(locA{b,2}(2)/2)],'Color','k','LineWidth',2);
            text(round(((locA{a,1}(1)+locA{a,2}(1)/2)+(locA{b,1}(1)+locA{b,2}(1)/2))/2),...
                    round(((locA{a,1}(2)+locA{a,2}(2)/2)+(locA{b,1}(2)+locA{b,2}(2)/2))/2)-10,...
                    num2str(tmp));
        end
    end
end
scnsize = get(0,'MonitorPosition');
set(gcf, 'position', [0 0 scnsize(3)/2 scnsize(4)]);
end
function [A]=myLine(A,from,to,val)
%在矩阵A中画直线
k=(to(2)-from(2))/(to(1)-from(1));
b=from(2)-k*from(1);
if from(1)>to(1)
    for x=to(1):from(1)
        y=round(k*x+b);
        A(y,x)=val;
    end
else
    for x=from(1):to(1)
        y=round(k*x+b);
        A(y,x)=val;
    end
end
end
