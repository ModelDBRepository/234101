function dynamicDriveChange(alpha,varnames,from,to,tstep2)
% runs simulation for tstep ms and the changes variables varnames from
% value from to value to. If tstep2 is provided the value will be set back
% to the original one at tstep2. If varnamse is 'alpha' alpha will be
% changed.
% This script allows recration of the simulation in Figure 8. 
% here are some examples:
% Figure 8A dynamicDriveChange(0.02,{'alpha'},0.02,0.4)
% Figure 8B dynamicDriveChange(0.02,{'alpha'},0.85,0.6)
% Figure 8C dynamicDriveChange(0.02,{'alpha'},0.02,0.9,10500)
% Figure 8G dynamicDriveChange(0.5,{'V0Vlrdrive_off','V0Vfhdrive_off'},0,2)

executable = './icpg';
filename = './models/Danner-etal-eLife.txt';

tstep=10000;
if iscell(varnames)
    varname = varnames{1};
else
    varname = varnames;
end
varnamestr=varname;

if ~exist('dontrun','var')
    call=[executable, ' -f ', strrep(filename,' ','\ '), ' -a ', num2str(alpha), ' -U ' , varname,  ' ', num2str(tstep), ' ', num2str(from), ' ', num2str(to)];
    if exist('tstep2','var')
        call=[call, ' ', num2str(tstep2)];
        
    end
    if iscell(varnames)&&length(varnames)>1
        for i=2:length(varnames)
            call=[call, ' -V ', varnames{i}];
            varnamestr=[varnamestr, ' ', varnames{i}];
        end
    end
    system(call);
end
xn=h5read('./results/example.cdf','/data');
x=xn';
Vmin=-50;
Vmax=0;
Vs=x(:,2:end);
Vs=(Vs-Vmin)/(Vmax-Vmin);
Vs(Vs<=0)=0;
Vs(Vs>=1)=1;


PLOTNEURONS={'RGF_NaP_R_front','RGF_NaP_L_front','RGF_NaP_R_hind','RGF_NaP_L_hind'};
DURPLOTCLICK=3000;
xx=tstep-1000;

offi=1:4;
colors2=['bgrkyyyymmmmcccc'];


fid = fopen(filename);
indices=ones(length(PLOTNEURONS),1)*-1;
tline = fgets(fid);
neurons={};
nlineneuron=-1;
while ischar(tline)
    strs_temp = strsplit(tline);
    strs={};
    for st = 1:length(strs_temp)
        if ~isempty(strs_temp{st})
            strs{end+1}=strs_temp{st};
        end
    end
    if isempty(strs)
        tline = fgets(fid);
        continue
    end
    if strcmp(strs{1}, 'neuron')
        nlineneuron=nlineneuron+1;
        neurons{end+1,1}=nlineneuron;
        if length(strs)==3
            neurons{end,2}=strs{3};
        elseif length(strs)==2
            neurons{end,2}=strs{2};
        end
        for i=1:length(PLOTNEURONS)
            if length(strs)==3
                if strcmp(strs{3},PLOTNEURONS{i})
                    indices(i)=nlineneuron;
                    
                end
            elseif length(strs)==2
                if strcmp(strs{2},PLOTNEURONS{i})
                    indices(i)=nlineneuron;
                    
                end
            end
        end
        
    end
    tline = fgets(fid);
end
fclose(fid);
indices=indices+1;
for i = 1:length(PLOTNEURONS)
    PLOTNEURONS{i} = strrep(PLOTNEURONS{i}, '_', ' ');
    PLOTNEURONS{i} = strrep(PLOTNEURONS{i}, 'NaP', '');
    PLOTNEURONS{i} = strrep(PLOTNEURONS{i}, 'hind', 'h');
    PLOTNEURONS{i} = strrep(PLOTNEURONS{i}, 'front', 'f');
end

index = indices;
maxI=length(index);

fig=figure();
set(fig,'RendererMode','manual') ;
set(fig,'Renderer','painters');
excess=DURPLOTCLICK*0.5;

newx=x(x(:,1)>=(xx-excess)&x(:,1)<=(xx+DURPLOTCLICK+excess),:);
newVs=Vs(x(:,1)>=(xx-excess)&x(:,1)<=(xx+DURPLOTCLICK+excess),:);
subplot(11,1,1:5);

hold off;
for ij = 1:maxI
    plot(newx(:,1),newx(:,index(ij)+1)-(offi(ij)-2)*75-25,[colors2(ij),'-.']);hold on;
    plot(newx(:,1),(newVs(:,index(ij))*(Vmax-Vmin))-(offi(ij)-1)*75,colors2(ij));hold on;
    xlim([newx(1,1)+excess, newx(1,1)+DURPLOTCLICK+excess]);
    
    set(gca,'YTickLabel',PLOTNEURONS(length(PLOTNEURONS):-1:1),'YTick',-75*(length(PLOTNEURONS)-1):75:0);
end
for ij = DURPLOTCLICK/10:DURPLOTCLICK/10:DURPLOTCLICK*1.1+newx(1,1)+excess
    line([ij,ij],[-75*(length(PLOTNEURONS)-1) 50],'LineWidth',.5,'LineStyle',':','Color','k'); hold on;
end
line([tstep,tstep],[-75*(length(PLOTNEURONS)-1)-50 50],'LineWidth',1,'LineStyle','-','Color','k'); hold on;

ylim([-75*(length(PLOTNEURONS)-1)-50 50]);
title(['at alpha = ' , num2str(alpha),' ',varnamestr, ' from ', num2str(from), ' to ', num2str(to)]);
subplot(11,1,6:8);
hold off;
for iii =1:4
    [ons,~,~,offs]=findBursts(newx(:,1),newVs(:,index(iii)),.1);
    for jj=1:min([length(ons)-1,length(offs)])
        rectangle('Position',[offs(jj),-iii-.25,ons(jj+1)-offs(jj),.5],'EdgeColor',colors2(iii),'FaceColor',colors2(iii))
        set(gca,'YTickLabel',{'lHl','rHl','lFl','rFl'},'YTick',[ -4 -3 -2 -1]);
    end
end
for ij = (DURPLOTCLICK/10:DURPLOTCLICK/10:DURPLOTCLICK*1.1)+newx(1,1)+excess
    line([ij,ij],[-4.5 -0.5],'LineWidth',.5,'LineStyle',':','Color','k'); hold on;
end
line([tstep,tstep],[-4.5 -0.5],'LineWidth',1.,'LineStyle','-','Color','k'); hold on;

xlim([newx(1,1)+excess, newx(1,1)+DURPLOTCLICK+excess]);
ylim([-4.5 -0.5]);
title('extension phases')

subplot(11,1,9:10);
for iii=1:4
    phase=getPhases(index(4),index(iii),newx,newVs);
    plot(phase(:,2),phase(:,1),'Marker','.','MarkerSize',10,'LineWidth',1,'Color',colors2(iii));hold on;
end
phase=getPhases(index(2),index(1),newx,newVs);
plot(phase(:,2),phase(:,1),'Marker','.','MarkerSize',10,'LineWidth',1,'Color','m');hold on;
xlim([newx(1,1)+excess, newx(1,1)+DURPLOTCLICK+excess]);
ylim([-0.25,1.25]);
title('phase differences')

subplot(11,1,11);
[~,~,pdur,on]=findBursts(x(:,1),Vs(:,index(4)),.1);
if length(on)>length(pdur)
    on=on(1:length(pdur));
elseif length(pdur)>length(on)
    pdur=pdur(1:length(on));
end

plot(on+pdur/2,1000./pdur,'Marker','.','MarkerSize',10,'LineWidth',1,'Color','k');
xlim([newx(1,1)+excess, newx(1,1)+DURPLOTCLICK+excess]);
ylim([0 12])
title('frequency')

    function phases =  getPhases(i1,i2,x,Vs)
        [~,~,pdur1,on1]=findBursts(x(:,1),Vs(:,i1),.1);
        [~,~,~,on2]=findBursts(x(:,1),Vs(:,i2),.1);
        [~,II]=sort([on1,on2]);
        ONS=[on1,on2;zeros(1,length(on1)),ones(1,length(on2));1:length(on1),1:length(on2)];
        ONS=ONS(:,II);
        inons=strfind(ONS(2,:),[0,1]);
        phases=ones(length(inons)-1,2)*-1;
        for ii=1:length(inons)-1
            phases(ii,1)=(ONS(1,inons(ii)+1)-ONS(1,inons(ii)))/pdur1(ONS(3,inons(ii)));
            phases(ii,2)=ONS(1,inons(ii))+pdur1(ONS(3,inons(ii)))/2;
        end
        for iji = 2:length(phases)-1
            if phases(iji-1,1) > 0.4 && phases(iji,1) < -0.4 && phases(iji+1,1) > 0.4
                phases(iji,1)=phases(iji,1)+1;
            end
        end
        
    end
    function [on,bdur,pdur,neg] = findBursts(t,V,thr)
        dV = (V>=thr)';
        inon=strfind(dV,[0,1,1]);
        on=zeros(1,length(inon));
        for j=1:length(inon)
            off=V(inon(j));
            k=(V(inon(j)+1)-V(inon(j)))/(t(inon(j)+1)-t(inon(j)));
            on(j)=(thr-off)/k+t(inon(j));
        end
        inneg=strfind(dV,[1,0,0]);
        neg=zeros(1,length(inneg));
        for j=1:length(inneg)
            off=V(inneg(j));
            k=(V(inneg(j)+1)-V(inneg(j)))/(t(inneg(j)+1)-t(inneg(j)));
            neg(j)=(thr-off)/k+t(inneg(j));
        end
        
        if neg(1)<on(1)
            neg=neg(2:end);
        end
        ln=min([length(on),length(neg)]);
        neg=neg(1:ln);
        on=on(1:ln);
        bdur=neg-on;
        pdur=on(2:ln)-on(1:ln-1);
        
    end
end