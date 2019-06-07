% =========================================================================
% Lorentzian line fitting test
% =========================================================================

%% Inputs
% Constant
linewidthTK = 2;
linewidthTN = 0.75;
fontsz = 11.5;

% -------------------------------------------------------------------------
% Inputs
flag_debug = 1;
% dpath = 'C:\Lab_Desktop\Code\LorentzianFitting';
dnamem = 'iopamidol20mM_1stPointRemoved.mat';
rnamem = 'iopamidol20mM_1stPointRemoved_results.mat';
pHsel = [6.05,6.47,7.97];

nlline = 5;
llinestr = {'water','iop0.8','iop1.8','iop4.2','iop5.6'};
c0shift = [0,0.8,1.8,4.2,5.6];
cswaterub = 1;      % fitting upper bound for water
options = optimset('Display','notify','MaxFunEvals',100000,...
        'MaxIter',100000,'TolFun',1e-8,'TolX',1e-8,...
        'Algorithm','trust-region-reflective');
% -------------------------------------------------------------------------

%% Initialization
% load([dpath,'\',dnamem]);
load(dnamem);
if ~flag_debug
%     save([dpath,'\',rnamem],'nlline','llinestr','c0shift','cswaterub');
    save(rnamem,'nlline','llinestr','c0shift','cswaterub');
end

npHsel = length(pHsel);
ipHsel = zeros(npHsel,1);
for i = 1:npHsel
    ipHsel(i) = find(pH == pHsel(i),1);
end

Mzsel = squeeze(Mz(ipHsel,:));
Mzselnorm = zeros(size(Mzsel));
for i = 1:npHsel
    Mzselnorm(i,:) = Mzsel(i,:)./Mzsel(i,1);
end

%% Lorentzian fitting
for i = 1:npHsel
    c0 = zeros(1,3*nlline);
    for j = 1:nlline
        c0(3*j) = c0shift(j);
    end
    
    zspec = squeeze(Mzselnorm(i,:));
    [zspecf,lline,res,cf,c0,lb,ub] = lorentzianfit(zspec(:),cs(:),...
        nlline,c0,cswaterub,'llinestr',llinestr,'dfbr',0.5,'shiftcorr',1,...
        'updatewidth',0,'showfig',0,'unit','ppm');
    zspecf = 1-zspecf;
    
    % results display and save
    figure;
    hold all;
    plot(cs,zspec,'bo','LineStyle','none');
    plot(cs,zspecf,'r-','LineWidth',linewidthTK);
    plot(cs,res,'k-','LineWidth',linewidthTK);
    for j = 1:nlline
        plot(cs,1-squeeze(lline(j,:)),'-');
    end
    hold off;
    legend([{'Measured'},{'Fitted'},{'Residual'},llinestr(:)'],'Location','east');
    xlabel('Saturation frequency (ppm)');
    set(gca,'XDir','reverse');
    ylabel('Signal (au)');
    title(['Lorentzian fitting (',num2str(nlline),' lines,',...
        ' iop pH ',num2str(pHsel(i)),')']);
end

%%
% =========================================================================
% 20190228 SZ: 1st version. 
% 20190228 SZ: remove the 1st point in the spectrum.
% 20190606 SZ: use lorentzianfit function.
% =========================================================================