function [  ] = create_noise_diags(sim,alpha, dontrun)
% Calculates distribution of phase differences when strong noise is applied
% to the model. Phase differences are then categorized them into three 
% equally sized bins (cf. Figure 7).
% Tse argument sim to specify which model shall be simulated
%       1: intact
%       2: all V0V ablated
%       3: diagonal V0V ablated
%       4: all V0V and V0D ablated
%       5: descending LPNs ablated
% alpha specifies the alpha value.
% 
% Examples
% Figure 7C:
% create_noise_diags(1,0.3)
% create_noise_diags(5,0.3)
% Figure 7D:
% create_noise_diags(1,0.6)
% create_noise_diags(5,0.6)
% Figure 7E:
% create_noise_diags(1,0.75)
% create_noise_diags(5,0.75)
% 
executable = './icpg'; % this needs to point to the executable
filename = './models/Danner-etal-eLife-noise.txt';

if ~exist('sim','var')
    sim = 1; % if not otherwise specified simulate the intact model
end
switch sim
    case 1
        % Simulate intact model
        argument = '';
    case 2
        % Simulate ablation of all V0V
        argument = ' -u V0VtoRGFdiagfh 0.0 -u V0VtoRGFdiaghf 0.0 -u inV0VtoRGF 0.0';
    case 3
        % Simulate ablation of only diagonal V0V
        argument = ' -u V0VtoRGFdiagfh 0.0 -u V0VtoRGFdiaghf 0.0';
    case 4
        % Simulate ablation of all V0V and V0D neurons
        argument = ' -u V0VtoRGFdiagfh 0.0 -u V0VtoRGFdiaghf 0.0 -u inV0VtoRGF 0.0 -u V0DtoRGFdiagfh 0.0 -u V0DtoRGF 0.0';
    case 5
        % Simulate ablation of all descending LPNs
        argument = ' -u V0VtoRGFdiagfh 0.0 -u V0DtoRGFdiagfh 0.0 -u inFH 0.0 -u V2aHomfh 0.0';
    otherwise
        argument = '';
end
if ~exist('dontrun','var')
    system([executable, ' -f ', strrep(filename,' ','\ '), ' -a ', num2str(alpha), ' ', argument]);
end

xn=h5read('./results/example.cdf','/data');
xn=xn';

x1=xn(xn(:,1)<=(max(xn(:,1))/2),:);
x2=xn(xn(:,1)>=(max(xn(:,1))/2),:);
x2(:,1)=x2(:,1)-min(x2(:,1));

Vmin=-50;
Vmax=0;

fid = fopen(filename);
tline = fgets(fid);
neurons={};
alphamax=0.93;
alphamin=0.0;
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
    if strcmp(strs{1}, 'stepwise')
        nSteps = sscanf(strs{2},'%d');
        stepDur = sscanf(strs{3},'%d');
        
        if length(strs)>=5 && ~isempty(sscanf(strs{5},'%f'))
            alphamin = sscanf(strs{4},'%f');
            alphamax = sscanf(strs{5},'%f');
        else
            alphamax = sscanf(strs{4},'%f');
        end
    end
    if strcmp(strs{1}, 'neuron')
        nlineneuron=nlineneuron+1;
        neurons{end+1,1}=nlineneuron;
        if length(strs)==3
            neurons{end,2}=strs{3};
        elseif length(strs)==2
            neurons{end,2}=strs{2};
        end
    end
    tline = fgets(fid);
end
fclose(fid);

N=length(neurons);
Vs=x1(:,2:N+1);

%Vs=1./(1+exp(-((Vs+30)./3)));
Vs=(Vs-Vmin)/(Vmax-Vmin);
Vs(Vs<=0)=0;
Vs(Vs>=1)=1;

Vs2=x2(:,2:N+1);

%Vs=1./(1+exp(-((Vs+30)./3)));
Vs2=(Vs2-Vmin)/(Vmax-Vmin);
Vs2(Vs2<=0)=0;
Vs2(Vs2>=1)=1;

if 1==1
    figure();
    plr=getPhases(1,2,x1,Vs,1);
    plrf=getPhases(5,6,x1,Vs,1);
    phom=getPhases(1,5,x1,Vs,1);
    pdiag=getPhases(1,6,x1,Vs,1);
    plr2=getPhases(1,2,x2,Vs2,2);
    plrf2=getPhases(5,6,x2,Vs2,2);
    phom2=getPhases(1,5,x2,Vs2,2);
    pdiag2=getPhases(1,6,x2,Vs2,2);
    plr=[plr;plr2];
    plrf=[plrf;plrf2];
    phom=[phom;phom2];
    pdiag=[pdiag;pdiag2];
    edges=0:0.125/8:1;
    plrh=histcounts(plr(:,1),edges);
    plrfh=histcounts(plrf(:,1),edges);
    phomh=histcounts(phom(:,1),edges);
    pdiagh=histcounts(pdiag(:,1),edges);
    mp=max([plrh,plrfh,phomh,pdiagh]);
    subplot(3,2,1);
    histogram(plr(:,1),edges);
    %polar(edges(1:end-1)*pi*2+pi/2,plrh)
    ylim([0 mp*1.1]);
    xlim([0 1]);
    title('left-right hind');
    subplot(3,2,2);
    histogram(plrf(:,1),edges);
    %polar(edges(1:end-1)*pi*2+pi/2,plrfh)
    ylim([0 mp*1.1]);
    xlim([0 1]);
    title('left-right fore');
    subplot(3,2,3);
    histogram(phom(:,1),edges);
    %polar(edges(1:end-1)*pi*2+pi/2,phomh)
    ylim([0 mp*1.1]);
    xlim([0 1]);
    title('fore-hind');
    subplot(3,2,4);
    histogram(pdiag(:,1),edges);
    title('diagonal');
    %polar(edges(1:end-1)*pi*2+pi/2,pdiagh)
    ylim([0 mp*1.1]);
    xlim([0 1]);
    
    subplot(3,2,[5,6]);
    barval=[sum(plr(:,1)<=0.167|plr(:,1)>=0.833),...
        sum(plr(:,1)>0.167&plr(:,1)<0.333)+...
        sum(plr(:,1)>0.667&plr(:,1)<0.833),...
        sum(plr(:,1)<=0.667&plr(:,1)>=0.333);...
        sum(plrf(:,1)<=0.167|plrf(:,1)>=0.833),...
        sum(plrf(:,1)>0.167&plrf(:,1)<0.333)+...
        sum(plrf(:,1)>0.667&plrf(:,1)<0.833),...
        sum(plrf(:,1)<=0.667&plrf(:,1)>=0.333);...
        sum(phom(:,1)<=0.167|phom(:,1)>=0.833),...
        sum(phom(:,1)>0.167&phom(:,1)<0.333)+...
        sum(phom(:,1)>0.667&phom(:,1)<0.833),...
        sum(phom(:,1)<=0.667&phom(:,1)>=0.333);
        sum(pdiag(:,1)<=0.167|pdiag(:,1)>=0.833),...
        sum(pdiag(:,1)>0.167&pdiag(:,1)<0.333)+...
        sum(pdiag(:,1)>0.667&pdiag(:,1)<0.833),...
        sum(pdiag(:,1)<=0.667&pdiag(:,1)>=0.333)
        ];
    bar1=bar(bsxfun(@rdivide,barval,sum(barval,2)));
    set(gca, 'XTick', 1:4, 'XTickLabel', {'left-right hind','left-right fore','homolateral','diagonal'});
    set(bar1(1),'DisplayName','Synchronized');
    set(bar1(2),'DisplayName','in-between');
    set(bar1(3),'DisplayName','Alternating');
    legend(gca,'show');
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

    function phases =  getPhases(i1,i2,x,Vs,Nrun)
        [on1,~,pdur1]=findBursts(x(:,1),Vs(:,i1),.1);
        [on2,~,~]=findBursts(x(:,1),Vs(:,i2),.1);
        if Nrun
            [~,II]=sort([on1,on2]);
        else
            [~,II]=sort([on1,on2],'descend');
        end
        ONS=[on1,on2;zeros(1,length(on1)),ones(1,length(on2));1:length(on1),1:length(on2)];
        ONS=ONS(:,II);
        inons=strfind(ONS(2,:),[0,1,0]);
        phases=ones(length(inons),2)*-1;
        for ii=1:length(inons)
            phases(ii,1)=(ONS(1,inons(ii)+1)-ONS(1,inons(ii)))/pdur1(ONS(3,inons(ii)));
            phases(ii,2)=ONS(1,inons(ii));
        end
        if (Nrun==2)
            phases(:,2)=-phases(:,2)+max(x(:,1));
        end
    end
end

