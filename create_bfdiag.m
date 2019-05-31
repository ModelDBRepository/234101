function [  ] = create_bfdiag(sim, dontrun)
% runs network simulation and plots bifurcation diagrams for all phase-
% differences, phase durations and frequency.
% use argument sim to specify which model shall be simulated
%       1: intact
%       2: all V0V ablated
%       3: diagonal V0V ablated
%       4: all V0V and V0D ablated
%       5: descending LPNs ablated
%
% Note: Simulations are by default run alpha being increased and then
% decreased in 50 equally spaced steps. This is done to keep simulation
% time low, but does not correspond to the figures in the paper. To create
% the same simulation as in the paper, the number 50 after stepwise on line
% 5 in Danner-etal-eLife.txt needs to be replaced by 500.

executable = './icpg';
filename = './models/Danner-etal-eLife.txt';

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
    system([executable, ' -f ', strrep(filename,' ','\ '), argument]);
end

xn=h5read('./results/example.cdf','/data');
xn=xn';
PLOTPHASELINES=1;
SPLITLINESPRECISION=0.025;

x1=xn(xn(:,1)<=(max(xn(:,1))/2),:);
x2=xn(xn(:,1)>=(max(xn(:,1))/2),:);
x2(:,1)=x2(:,1)-min(x2(:,1));

Vmin=-50;
Vmax=0;

fid = fopen(filename);
tline = fgets(fid);
neurons={};
nSteps=0;
stepDur=1;
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
if alphamin==alphamax
    alphamin = 0;
    alphamax = 1;
end
fclose(fid);

N=length(neurons);
createPlot(x1,1);
createPlot(x2,2);

    function createPlot(x,Nrun)
        Vs=x(:,2:N+1);
        Vs=(Vs-Vmin)/(Vmax-Vmin);
        Vs(Vs<=0)=0;
        Vs(Vs>=1)=1;
        
        subplot(6,1,1);
        [on,bdur,pdur,~]=findBursts(x(:,1),Vs(:,1),.1);
        if (Nrun==2)
            on = -on+max(x(:,1));
        end
        cats = floor(on/stepDur);
        cats=cats(1:min([length(pdur),length(bdur),length(cats)]));
        
        phasefq=zeros(4,max(cats)+1);
        for kk=0:max(cats)
            
            os=on(cats==kk);
            bd=bdur(cats==kk);
            phd=pdur(cats==kk);
            if length(os)>=2
                phasefq(1,kk+1)=abs(1000./(os(end)-os(end-1)));
                phasefq(2,kk+1)=bd(end-1);
                phasefq(3,kk+1)=phd(end-1)-bd(end-1);
                phasefq(4,kk+1)=kk*(alphamax-alphamin)/nSteps+alphamin;
            end
        end
        phasefq=phasefq(:,phasefq(1,:)~=0);
        plot(phasefq(4,:),phasefq(2,:),'r');hold on;
        plot(phasefq(4,:),phasefq(3,:),'b');
        xlim([alphamin alphamax])
        title(['burst duratation:']);
        
        subplot(612)
        plotPhase(1,2,x,Vs,Nrun,1);
        title('left-right hind phase (lHl-rHl)');
        
        subplot(613)
        plotPhase(5,6,x,Vs,Nrun,1);
        title('left-right fore phase (lFl-rFl)');
        
        subplot(614)
        
        plotPhase(1,5,x,Vs,Nrun,0);
        plotPhase(2,6,x,Vs,Nrun,0);
        title('homolateral phase (lHl-lFl)');
        
        subplot(615)
        plotPhase(1,6,x,Vs,Nrun,0);
        plotPhase(2,5,x,Vs,Nrun,0);
        title('diagonal phase (lHl-rFl)');
        
        subplot(616);
        [on1,~,~]=findBursts(x(:,1),Vs(:,1),.1);
        if (Nrun==2)
            on1 = -on1+max(x(:,1));
        end
        plot(on1(1:end-1),abs(1000./diff(on1))); hold on;
        xlim([0 max(x(:,1))]);
        title('frequency')
        
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
    function plotPhase(i1,i2,x,Vs,Nrun,dosym)
        colors={'b','r'};
        phase2=getPhases(i1,i2,x,Vs,Nrun);
        
        cats2 = floor(phase2(:,2)/stepDur)';
        lines={[]};
        ilines=1;
        first=1;
        for ijk=0:nSteps-1
            catPh=phase2(cats2==ijk,1);
            if length(catPh)>=2
                nplot=10;
                if length(catPh)<=nplot
                    nplot=length(catPh)-2;
                end
                catPh=catPh(end-(nplot):end-1);
                [mph,~]=phaseStat(catPh);
                if PLOTPHASELINES==1
                    if first==1
                        lines{ilines}=[lines{ilines};ijk.*(alphamax-alphamin)/nSteps+alphamin,mph];
                        first=2;
                    else
                        
                        ln = lines{ilines}(end,2);
                        if abs(ln-mph)>SPLITLINESPRECISION && (abs(ln-mph)<(1-SPLITLINESPRECISION) || dosym~=1)
                            ilines=ilines+1;
                            lines{ilines}=[ijk.*(alphamax-alphamin)/nSteps+alphamin,mph];
                        else
                            if abs(ln-mph)>=(1-SPLITLINESPRECISION)
                                lines{ilines}=[lines{ilines};ijk.*(alphamax-alphamin)/nSteps+alphamin,1-mph];
                            else
                                lines{ilines}=[lines{ilines};ijk.*(alphamax-alphamin)/nSteps+alphamin,mph];
                            end
                        end
                    end
                else
                    plot(ijk.*(alphamax-alphamin)/nSteps+alphamin,mph,['.',colors{Nrun}]);hold on;
                    if(dosym)
                        plot(ijk.*(alphamax-alphamin)/nSteps+alphamin,1-mph,['.',colors{Nrun}]);hold on;
                    end
                end
            end
            
        end
        if PLOTPHASELINES==1
            for lns=1:length(lines)
                if size(lines{lns},1)==1
                    options=['.',colors{Nrun}];
                else
                    options=colors{Nrun};
                end
                plot(lines{lns}(:,1),lines{lns}(:,2),options);hold on;
                if(dosym)
                    plot(lines{lns}(:,1),1-lines{lns}(:,2),options);hold on;
                end
            end
        end
        
        %xlim([alphamin alphamax]);
        xlim([0 alphamax]);
        ylim([0 1.0]);
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
    function [mphases,r]=phaseStat(phases)
        phases=phases*2*pi;
        mid=mean(exp(1i*phases));
        mphases=wrapTo2Pi(atan2(imag(mid),real(mid)))./(2.*pi);
        r=sqrt(imag(mid)^2+real(mid)^2);
    end
end

