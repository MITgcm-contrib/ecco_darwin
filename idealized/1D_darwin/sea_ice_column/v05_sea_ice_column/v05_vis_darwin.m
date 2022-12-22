clc;clear;
cd run/

tid=24:24:17520;
k=1;
for i=tid
no3(:,k)=squeeze(rdmds('PTRACER02',i));
k=k+1;
end

k=1;
for i=tid
po4(:,k)=squeeze(rdmds('PTRACER05',i));
k=k+1;
end

k=1;
for i=tid
sio(:,k)=squeeze(rdmds('PTRACER07',i));
k=k+1;
end

k=1;
for i=tid
doc(:,k)=squeeze(rdmds('PTRACER08',i));
k=k+1;
end

k=1;
for i=tid
poc(:,k)=squeeze(rdmds('PTRACER12',i));
k=k+1;
end

for chlid=27:31
    k=1;
    for i=tid
        chl(:,k,chlid-26)=squeeze(rdmds(['PTRACER' num2str(chlid,'%02d')],i));
        k=k+1;
    end
end

dz=[5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00,...
  5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 10.00,...
  10.00, 10.00, 10.00, 383.26];
z=cumtrapz(dz);
t=tid*3600/24/60/60;

subplot(3,2,1)
contourf(t,-z,no3,100,'linecolor','none')
ylim([-100 0])
hold on
cb=colorbar('east','color','w');
ylabel(cb,'NO_3 (\muM)')
ylabel('Depth (m)')
caxis([2 15])

subplot(3,2,3)
contourf(t,-z,sio,100,'linecolor','none')
ylim([-100 0])
hold on
cb=colorbar('east','color','w');
ylabel(cb,'SiO_2 (\muM)')
ylabel('Depth (m)')
caxis([1 20])

subplot(3,2,5)
contourf(t,-z,po4,100,'linecolor','none')
ylim([-100 0])
hold on
cb=colorbar('east','color','w');
ylabel(cb,'PO_4 (\muM)')
ylabel('Depth (m)')
caxis([0.5 1.2])
xlabel('Days')

subplot(3,2,2)
contourf(t,-z,poc,100,'linecolor','none')
ylim([-100 0])
hold on
cb=colorbar('east','color','w');
ylabel(cb,'POC (\muM)')

subplot(3,2,4)
contourf(t,-z,doc,100,'linecolor','none')
ylim([-100 0])
hold on
cb=colorbar('east','color','w');
ylabel(cb,'DOC (\muM)')

subplot(3,2,6)
contourf(t,-z,sum(chl,3),100,'linecolor','none')
ylim([-100 0])
hold on
cb=colorbar('east','color','w');
ylabel(cb,'Chl (\muM)')
xlabel('Days')
% set(gcf,'position',[9.8000  122.6000  887.2000  597.6000])

cd ..
print(gcf,'1d_darwin_v5.png','-dpng','-r500')
