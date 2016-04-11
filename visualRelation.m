function [ ] = visualRelation( relationsPa,relationsPe,distForGroup,peopleNumA,peopleNumB,imgSetA,imgSetB)
%UNTITLED Summary of this function goes here
%使用说明，把下方代码复制到命令窗口中运行，然后搞一下。
%load('../peopleNum.mat');load('../data/imgSet')visualRelation( relations,distForGroup,distForGroupT,peopleNumA,peopleNumB,imgSetA,imgSetB)
% global relations;
% global dist;
% load('peopleNum.mat');
% load('imgSet2.mat');
% imA包括pic,no,person,feature,edgefeat,wN,wE。person包括patchNum,data{n,2}，num，patch{n,pN,2},patchfeature{n,pN},feature{n}
%distForGroup{1,1}是patch分配元胞,{1,2}是patch相似度矩阵,{2,1}是人分配元胞,{2,2}是人相似度矩阵
groupNum=162;
option=1;
peopleNumA2=peopleNumA;peopleNumA2(1)=0;
peopleNumB2=peopleNumB;peopleNumB2(1)=0;
addpath('../');
pN=imgSetA{1,1}.person.patchNum;
hh=0;ww=0;
for i=2:size(peopleNumA,2)
    peopleNumA2(i)=peopleNumA2(i-1)+peopleNumA(i-1);%peopleNumA2中每个元素储存的是前面已经有了多少个人
    peopleNumB2(i)=peopleNumB2(i-1)+peopleNumB(i-1);
end
countera=1;counterb=1;
[value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
global fig;
global haxis;
fig=figure(10);
aID=1;
bID=1;%idx(counterb);
wNBlank=zeros(1,pN*imgSetA{1,aID}.person.num);
wEBlank=zeros(1,nchoosek(imgSetA{1,aID}.person.num,2)*pN^2);
draw2( aID,bID,distForGroup,distForGroup);
[hh,ww]=drawn(imgSetA{1,aID},imgSetB{1,bID},relationsPa,relationsPe,peopleNumA2,peopleNumB2,option,wNBlank,wEBlank);
set(fig, 'WindowButtonDownFcn',@click_ceshi);
set(fig,'WindowKeyPressFcn',@keypressfcn,'WindowKeyReleaseFcn',@keyreleasefcn);
    function keypressfcn(h,evt)
        fprintf('************press \n');
        evt
        fprintf('************ \n');
    end
    function keyreleasefcn(h,evt)
        fprintf('************release \n');
        if strcmp(evt.Key,'uparrow')
            countera=countera-1;
            if countera<1
                countera=1;
            end
            counterb=1;
            [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
            aID=countera;bID=idx(counterb);
            wNBlank=zeros(1,imgSetA{1,aID}.person.patchNum*imgSetA{1,aID}.person.num);
            wEBlank=zeros(1,nchoosek(imgSetA{1,aID}.person.num,2)*pN^2);
        elseif strcmp(evt.Key,'downarrow')
            countera=countera+1;
            counterb=1;
            [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
            aID=countera;bID=idx(counterb);
            wNBlank=zeros(1,imgSetA{1,aID}.person.patchNum*imgSetA{1,aID}.person.num);
            wEBlank=zeros(1,nchoosek(imgSetA{1,aID}.person.num,2)*pN^2);
        elseif strcmp(evt.Key,'rightarrow')
            counterb=counterb+1;
            bID=idx(counterb);
        elseif strcmp(evt.Key,'leftarrow')
            counterb=counterb-1;
            if counterb<1
                counterb=1;
            end
            bID=idx(counterb);
        elseif strcmp(evt.Key,'space')
            counterb=find(idx==aID);
            bID=idx(counterb);
        elseif strcmp(evt.Key,'1')
            counterb=1;
            bID=idx(counterb);
        elseif strcmp(evt.Key,'2')
            counterb=2;
            bID=idx(counterb);
        elseif strcmp(evt.Key,'3')
            counterb=3;
            bID=idx(counterb);
        elseif strcmp(evt.Key,'l')%按L，在显示线和显示权值之间切换
            option=1+option;
            if option==3   
                wNBlank=zeros(1,pN*imgSetA{1,aID}.person.num);
                wEBlank=zeros(1,nchoosek(imgSetA{1,aID}.person.num,2)*pN^2);
            end
            if option==4
                wNBlank=zeros(1,pN*imgSetA{1,aID}.person.num);
                wEBlank=zeros(1,nchoosek(imgSetA{1,aID}.person.num,2)*pN^2);
            end
            if option==5
                option=0;
            end
        elseif option==3
            if strcmp(evt.Key,'w')
                for kk=1:imgSetA{1,aID}.person.patchNum*imgSetA{1,aID}.person.num
                    if wNBlank(kk)==1
                        imgSetA{1,aID}.wN(kk,:)=imgSetA{1,aID}.wN(kk,:)+0.02;
                    end
                end
                for j=1:groupNum
                    [relationsPa{aID,j},relationsPe{aID,j},distForGroup(aID,j),~]...
                        =calcuMatch( imgSetA{1,aID}.feature,imgSetB{1,j}.feature,imgSetA{1,aID}.edgefeat,imgSetB{1,j}.edgefeat,imgSetA{1,aID}.wN,imgSetA{1,aID}.wE,1,'eig',imgSetA{1,1}.person.patchNum);
                end
                [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
            elseif strcmp(evt.Key,'s')
                for kk=1:imgSetA{1,aID}.person.patchNum*imgSetA{1,aID}.person.num
                    if wNBlank(kk)==1
                        imgSetA{1,aID}.wN(kk,:)=imgSetA{1,aID}.wN(kk,:)-0.02;
                    end
                end
                for j=1:groupNum
                    [relationsPa{aID,j},relationsPe{aID,j},distForGroup(aID,j),~]...
                        =calcuMatch( imgSetA{1,aID}.feature,imgSetB{1,j}.feature,imgSetA{1,aID}.edgefeat,imgSetB{1,j}.edgefeat,imgSetA{1,aID}.wN,imgSetA{1,aID}.wE,1,'eig',imgSetA{1,1}.person.patchNum);
                end
                [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
            elseif strcmp(evt.Key,'x')
                imgSetA{1,aID}.wN(:,:)=ones(size(imgSetA{1,aID}.wN,1),size(imgSetA{1,aID}.wN,2));
                for j=1:groupNum
                    [relationsPa{aID,j},relationsPe{aID,j},distForGroup(aID,j),~]...
                        =calcuMatch( imgSetA{1,aID}.feature,imgSetB{1,j}.feature,imgSetA{1,aID}.edgefeat,imgSetB{1,j}.edgefeat,imgSetA{1,aID}.wN,imgSetA{1,aID}.wE,1,'eig',imgSetA{1,1}.person.patchNum);
                end
                [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
                wNBlank=zeros(1,pN*imgSetA{1,aID}.person.num);
                wEBlank=zeros(1,nchoosek(imgSetA{1,aID}.person.num,2)*pN^2);
            end
        elseif option==4
            if strcmp(evt.Key,'w')
                for kk=1:nchoosek(imgSetA{1,aID}.person.num,2)*pN^2
                    if wEBlank(kk)==1
                        imgSetA{1,aID}.wE(kk,:)=imgSetA{1,aID}.wE(kk,:)+0.5;
                    end
                end
                for j=1:groupNum
                    [relationsPa{aID,j},relationsPe{aID,j},distForGroup(aID,j),~]...
                        =calcuMatch( imgSetA{1,aID}.feature,imgSetB{1,j}.feature,imgSetA{1,aID}.edgefeat,imgSetB{1,j}.edgefeat,imgSetA{1,aID}.wN,imgSetA{1,aID}.wE,1,'eig',imgSetA{1,1}.person.patchNum);
                end
                [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
            elseif strcmp(evt.Key,'s')
                for kk=1:nchoosek(imgSetA{1,aID}.person.num,2)*pN^2
                    if wEBlank(kk)==1
                        imgSetA{1,aID}.wE(kk,:)=imgSetA{1,aID}.wE(kk,:)-0.5;
                    end
                end
                for j=1:groupNum
                    [relationsPa{aID,j},relationsPe{aID,j},distForGroup(aID,j),~]...
                        =calcuMatch( imgSetA{1,aID}.feature,imgSetB{1,j}.feature,imgSetA{1,aID}.edgefeat,imgSetB{1,j}.edgefeat,imgSetA{1,aID}.wN,imgSetA{1,aID}.wE,1,'eig',imgSetA{1,1}.person.patchNum);
                end
                [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
            elseif strcmp(evt.Key,'x')
                imgSetA{1,aID}.wE(:,:)=ones(size(imgSetA{1,aID}.wE,1),size(imgSetA{1,aID}.wE,2));
                for j=1:groupNum
                    [relationsPa{aID,j},relationsPe{aID,j},distForGroup(aID,j),~]...
                        =calcuMatch( imgSetA{1,aID}.feature,imgSetB{1,j}.feature,imgSetA{1,aID}.edgefeat,imgSetB{1,j}.edgefeat,imgSetA{1,aID}.wN,imgSetA{1,aID}.wE,1,'eig',imgSetA{1,1}.person.patchNum);
                end
                [value,idx] = sort(distForGroup(countera,:),'ascend');%竹排排序
                wNBlank=zeros(1,pN*imgSetA{1,aID}.person.num);
                wEBlank=zeros(1,nchoosek(imgSetA{1,aID}.person.num,2)*pN^2);
            end
        end
        draw2( aID,bID,distForGroup,distForGroup);
        [hh,ww]=drawn(imgSetA{1,aID},imgSetB{1,bID},relationsPa,relationsPe,peopleNumA2,peopleNumB2,option,wNBlank,wEBlank);
        fprintf('************ \n');
    end
    function click_ceshi(src, event)
        xy = get(fig, 'CurrentPoint');
        %set(gca, 'Unit','Pixel');
        %sizeGca=get(gca, 'position');%gca的高永远与figure一样，
        sizeFig=get(fig, 'position');
        xy(2)=sizeFig(4)-xy(2);
        if (sizeFig(4)/hh)>(sizeFig(3)/ww)
            radio=sizeFig(3)/ww;
        else
            radio=sizeFig(4)/hh;
        end
        hh=hh*radio;ww=ww*radio;
        hMargin=(sizeFig(4)-hh)/2;wMargin=(sizeFig(3)-ww)/2;
        xy(1)=(xy(1)-wMargin)/radio;xy(2)=(xy(2)-hMargin)/radio;
        for ii=1:imgSetA{1,aID}.person.num
            for jj=1:imgSetA{1,aID}.person.patchNum
                tmpx=xy(1)-imgSetA{1,aID}.person.patch{ii,jj,1}(1);
                tmpy=xy(2)-imgSetA{1,aID}.person.patch{ii,jj,1}(2);
                if ((tmpx<imgSetA{1,aID}.person.patch{ii,jj,2}(1))&&(tmpx>0))&&( (tmpy<imgSetA{1,aID}.person.patch{ii,jj,2}(2))&&(tmpy>0) )
                    if option==3
                        if wNBlank((ii-1)*pN+jj)==1
                            wNBlank((ii-1)*pN+jj)=0;
                        else
                            wNBlank((ii-1)*pN+jj)=1;
                        end
                    elseif option==4
                        loc1=find(wNBlank==1);
                        if isempty(loc1) 
                            wNBlank((ii-1)*pN+jj)=1;
                        else
                            if wNBlank((ii-1)*pN+jj)==1;
                                wNBlank((ii-1)*pN+jj)=0;
                            else
                                wNBlank(loc1)=0;
                                [~,loc2]=find(imgSetA{1,aID}.edgefeat(1:2,:)==((ii-1)*pN+jj));
                                [~,loc3]=find(imgSetA{1,aID}.edgefeat(1:2,loc2)==loc1);
                                if wEBlank(loc2(loc3))==0
                                    wEBlank(loc2(loc3))=1;
                                else
                                    wEBlank(loc2(loc3))=0;
                                end
                            end
                        end
                    end
                    break;
                end
            end
        end
        draw2( aID,bID,distForGroup,distForGroup);
        [hh,ww]=drawn(imgSetA{1,aID},imgSetB{1,bID},relationsPa,relationsPe,peopleNumA2,peopleNumB2,option,wNBlank,wEBlank);
    end
end


