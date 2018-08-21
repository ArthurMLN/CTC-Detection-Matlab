clear all;close all;
tic; 
X2=0;
X2='5xy023';
nameD=['F:\Matlab\image\',X2,'c3.tif'];
nameF=['F:\Matlab\image\',X2,'c2.tif'];
nameT=['F:\Matlab\image\',X2,'c1.tif'];
nameO=['F:\Matlab\Cimage\5\',X2,'.tif'];
% D=DAPI original image
D=imread(nameD);
% F=FITC original image
F=imread(nameF);
% T=TRTC original image
T=imread(nameT);
% O=Overlay 
O=imread(nameO);
OD=O;
OF=O;
OT=O;
OC=O;
% Parameter
[m,n]=size(D);


% Clean
D1=zeros(size(D));
D2=zeros(size(D));
D1=uint8(D1);
D2=uint8(D2);
X13=find(D>30);
[p,q]=size(X13);
if (p*q)>(0.3*m*n)
    disp('no CTC1');
    break;
end
D2(D>200)=1;
L3 = bwlabel(D2,8);
Pro=0;
A=0;
Pro=regionprops(L3,'area');
p=size(Pro);
for i=1:p
    A(i)=Pro(i).Area;
end
X6=find(A<100);
[X7,X8]=size(X6);
for i=1:X8
    L3(L3==X6(i))=0;
end
X6=find(A>600);
[X7,X8]=size(X6);
for i=1:X8
    L3(L3==X6(i))=0;
end
D2=zeros(m,n);
D2=uint8(D2);
D2(D<200)=D(D<200);
D1=im2bw(D2,graythresh(D2));
D1(L3>0)=1;
D2(D1==1)=D(D1==1);



F1=zeros(size(F));
F2=F;
F2(find((F2-D2>50)))=0;
F2(D2<25)=0;
F2(F2>200)=0;
if  max(max(F2))<25
    disp('no CTC2');
    break;
end
F1=im2bw(F2,graythresh(F2));



T1=zeros(size(T));
T2=T;
T2(D2==0)=0;
T2(T2>200)=0;
if  max(max(T2))<25
    disp('no CTC3');
    break;
end
T1=im2bw(T2,graythresh(T2));



% % Delete meaningless DAPI
X=F1+T1;
D1(X==0)=0;
% 
% % Find CTC core,by DAPI and FITC positive
X2=F1+D1;
C1=zeros(m,n);

C1(X2>1)=1;

% % the fitc shall be stronger than tricy at lease 20%,later 50%
C1(F<T)=0;

if max(max(C1))==0;
    disp('no CTC found4');
    break;  
end
   
% Find WBC core by DAPI and TRIC positive
X3=T1+D1;
B1=zeros(m,n);
B1(X3>1)=1;
B1(F1==1)=0;
% caculate WBC FITC strength X4
X4=max(F(B1==1));
% caculatae CTC TRIC strength X5
X5=max(T(C1==1));



% Region label
L=0;
L = bwlabel(C1,8);
Pro=0;
Pro=regionprops(L,'area');
% p is the number of regions
[p,q]=size(Pro);
% because this is a nuclear,we clear area smaller than 250 px,and area
% bigger than 450
X6=0;
X7=0;
A=0;
X8=0;
for i=1:p
    A(i)=Pro(i).Area;
end
X6=find(A<100);
[X7,X8]=size(X6);
for i=1:X8
    L(L==X6(i))=-2;
end
C1(L==-2)=0;
X6=find(A>450);
[X7,X8]=size(X6);
for i=1:X8
    L(L==X6(i))=-2;
end
C1(L==-2)=0;
L=0;
L = bwlabel(C1,8);
if max(max(C1))==0;
    disp('no CTC found5');
    break;  
end

% region property
Pro=0;
Box=0;
Pro=regionprops(L,'area','Perimeter');
Box=regionprops(L,'BoundingBox');
[p,q]=size(Box);
A=0;
for i=1:p
    A(i)=Pro(i).Area;
end
for i=1:p
    P(i)=Pro(i).Perimeter;
end
for i=1:p
    Width(i,:)=Box(i).BoundingBox;
end
for i=1:p
WHRate(i)=abs((Width(i,3)/Width(i,4))-1);
end
L1=zeros(m,n);
% limit below WHRate
X6=find(WHRate>0.5);
[X7,X8]=size(X6);
for i=1:X8
    L1(L==X6(i))=-2;
end
% limit circularity
for i=1:p
     c(i)=((4*pi*A(i))/(P(i)^2));
     if ((4*pi*A(i))/(P(i)^2))<0.2;
     L1(L==i)=-2;
     end
end
C1(L1==-2)=0;
if max(max(C1))==0;
    disp('no CTC found6');
    break;  
end
% set L to the final area according to C1
L(C1<1)=0;
% The average strength limitation
for i=1:max(max(L))
    X10=F(L==i);
    X11=D(L==i);
    X12=T(L==i);
%     set area fitcy mean 20 to 100
    if mean(X10)<10
        C1(L==i)=0;
    end
     if mean(X10)>100
        C1(L==i)=0;
     end
%     set area fitcy mean 20% higher than tricy
    if (0.8*(mean(X10)))<mean(X12)
       C1(L==i)=0;
    end
%     set DAPI can't too low
    if mean(X11)<25
       C1(L==i)=0;
    end       
end
% limit DAPI can't too low than 20% of the average
for i=1:max(max(L))
      if mean(X11(i,:))<(0.2*mean(mean(X11)))
        C1(find(L==i))=0;
      end
end
if max(max(C1))==0;
    disp('no CTC found7');
    break;  
end


C1=imfill(C1,'holes');
C2=imfill(C1,'holes');
SE1=strel('disk',3);
SE2=strel('disk',1);
C1 = imdilate(C1,SE1);
C2 = imdilate(C2,SE2);
F4=zeros(m,n);
D4=zeros(m,n);
F4=F;
D4=D;
F4(C1<1)=0;
D4(C2<1)=0;
F4=im2bw(F4,graythresh(F4));
D4=im2bw(D4,graythresh(D4));
SE1=strel('disk',3);
SE2=strel('disk',1);
F4(F<T)=0;
F4=imerode(F4,SE1);
F4=imdilate(F4,SE1);
T4=mean(T(F4>0));
D4(F<T)=0;
D4=imdilate(D4,SE2);
D4=imerode(D4,SE2);
C1=F4;
C2=D4;
if max(max(C1))==0;
    disp('no CTC found8');
    break;  
end
% delet small area
L=0;
L2=0;
A=0;
Pro1=0;
Pro2=0;
L = bwlabel(F4,8);
L2 =bwlabel(D4,8);
for i=1:max(max(L))
    for j=1:max(max(L2))
         X1=find(L==i);
         X2=find(L2==j);
        if ~isempty(intersect(X1,X2));
             X3=find(L2==i);
             L2(L2==j)=i;
             L2(X3)=j;
         else
             L2(L2==j)=-1;             
        end
    end
end
L2(L2==-1)=0;
for i=1:max(L1)
    X4=find(L2==i);
    [p1,q1]=size(X4);
    if p1*q1==0
        L1(L1==i)=0;
    end
end
Pro1=regionprops(L,'area');
Pro2=regionprops(L2,'area');
[p,q]=size(Pro1);
for i=1:p
    if Pro1(i).Area<50
        F4(L==i)=0;
    end
end
[p,q]=size(Pro2);
for i=1:p
    if Pro2(i).Area<50
        D4(L==i)=0;
    end
end





% cell feature as whole analyze
L=0;
A=0;
Pro1=0;
Pro2=0;
L = bwlabel(F4,8);
L2 =bwlabel(D4,8);
for i=1:max(max(L))
    for j=1:max(max(L2))
         X1=find(L==i);
         X2=find(L2==j);
         if ~isempty(intersect(X1,X2));
             X3=find(L2==i);
             L2(L2==j)=i;
             L2(X3)=j;
         else
             L2(L2==j)=-1;             
        end
   end
end
L2(L2==-1)=0;
for i=1:max(L1)
    X4=find(L2==i);
    [p1,q1]=size(X4);
    if p1*q1==0
        L1(L1==i)=0;
    end
end
     
Pro1=regionprops(L,'area','Perimeter');
Pro2=regionprops(L2,'area','Perimeter','BoundingBox');
p2=0;
q2=0;
[p2,q2]=size(Pro1);
if size(Pro2)~=0
       for i=1:p2
    Width(i,:)=Pro2(i).BoundingBox;
end
end 

for i=1:p
%     index
    A(i,1)=i;
    A(i,2)=Pro1(i).Area;  
    A(i,3)=Pro1(i).Perimeter;
    X14=F(L==i);
    X15=T(L==i);
    A(i,4)=mean(X14); 
    A(i,5)=std2(X14);
    A(i,10)=mean(X15);
    A(i,6)=Pro2(i).Area;
    A(i,7)=A(i,6)/A(i,2);  
    A(i,8)=Pro2(i).Perimeter;
    A(i,9)=((4*pi*A(i,6))/(A(i,8)^2));
    A(i,11)=abs((Width(i,3)/Width(i,4))-1);
    
%    if A(i,2)<265
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
%     if A(i,3)<50
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
%     if 0.8*A(i,4)<A(i,10)
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
%     if A(i,5)<4
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
%      if A(i,6)>350
%        F4(L==i)=0;
%        D4(L2==i)=0;
%      end
%       if A(i,6)<190
%        F4(L==i)=0;
%        D4(L2==i)=0;
%       end
%       if A(i,7)<0.55
%        F4(L==i)=0;
%        D4(L2==i)=0;
%       end
%        if A(i,8)<45
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
%     if A(i,9)<1.03
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
%      if A(i,11)>0.35
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
%      if A(i,11)<1.05
%        F4(L==i)=0;
%        D4(L2==i)=0;
%     end
 end
C1=F4;
C2=D4;
if max(max(C1))==0;
    disp('no CTC found9');
    break;  
end


% Recorder;
L=0;
A=0;
Pro1=0;
Pro2=0;
L = bwlabel(F4,8);
L2 =bwlabel(D4,8);
for i=1:max(max(L))
    for j=1:max(max(L2))
         X1=find(L==i);
         X2=find(L2==j);
         if ~isempty(intersect(X1,X2));
             L2(L2==j)=i;
         end
   end
end
Pro1=regionprops(L,'area','Perimeter');
Pro2=regionprops(L2,'area','Perimeter');
for i=1:j
    S(i).Intensity=F(L==i);
end

[p,q]=size(Pro1);
for i=1:p
%     index
    A(i,1)=i*1000;
    A(i,2)=Pro1(i).Area;
    A(i,3)=Pro1(i).Perimeter;
    X14=F(L==i);
    X15=T(L==i);
    A(i,4)=sum(X14); 
    A(i,5)=std2(X14);
    A(i,6)=Pro2(i).Area;
    A(i,7)=A(i,6)/A(i,2);  
    A(i,8)=Pro2(i).Perimeter;
    A(i,9)=((4*pi*A(i,6))/(A(i,8)^2));
    A(i,10)=mean(X15);
end
% A=sortrows(A,-2);







% split channel display
D1e=edge(D1);
OD(D1e>0)=255;
OD(find(D1e>0)+m*n)=0;
OD(find(D1e>0)+2*m*n)=0;

F1e=edge(F1);
OF(F1e>0)=255;
OF(find(F1e>0)+m*n)=0;
OF(find(F1e>0)+2*m*n)=0;

T1e=edge(T1);
OT(T1e>0)=255;
OT(find(T1e>0)+m*n)=0;
OT(find(T1e>0)+2*m*n)=0;


% Overlay image display
C1e=zeros(m,n);
C2e=zeros(m,n);
C1e=edge(C1);
C2e=edge(C2);
B1e=B1;
OC(C1e>0)=255;
OC(find(C2e>0)+m*n)=255;
OC(find(C1e>0)+2*m*n)=0;




% Figure show
figure(1);subplot('Position',[0.01 0.7 0.3 0.30]);imshow(D);
figure(1);subplot('Position',[0.01 0.35 0.3 0.30]);imshow(D1);
figure(1);subplot('Position',[0.01 0.01 0.3 0.30]);imshow(OD);
figure(1);subplot('Position',[0.35  0.7 0.3 0.30]);imshow(F);
figure(1);subplot('Position',[0.35 0.35 0.3 0.30]);imshow(F1);
figure(1);subplot('Position',[0.35 0.01 0.3 0.30]);imshow(OF);
figure(1);subplot('Position',[0.7 0.7 0.3 0.30]);imshow(T);
figure(1);subplot('Position',[0.7 0.35 0.3 0.30]);imshow(T1);
figure(1);subplot('Position',[0.7 0.01 0.3 0.30]);imshow(OT);


figure(2);imshow(OC);
Time=toc;
disp(Time);

% figure(1),subplot(3,3,4);imshow(D1);
% figure(1),subplot(3,3,2);imshow(F);
% figure(1),subplot(3,3,5);imshow(F1);
% figure(1),subplot(3,3,3);imshow(T);
% figure(1),subplot(3,3,6);imshow(T1);
% figure(1),subplot(3,3,7);imshow(OD);
% figure(1),subplot(3,3,8);imshow(OF);
% figure(1),subplot(3,3,9);imshow(OT);
