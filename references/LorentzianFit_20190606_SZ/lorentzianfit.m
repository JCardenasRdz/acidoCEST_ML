function [yf,ylline,res,cf,c0,lb,ub] = lorentzianfit(zspec,off,nlline,c0,offwaterub,varargin)
% =========================================================================
% Lorentzian fitting for CEST spectrum
% 1. The first pool must be water pool.
% =========================================================================

% Code path
addpath(genpath('C:\Lab_Desktop\Code\General'));    % For func findc1D
addpath('C:\Lab_Desktop\Code\CPE');                 % For func FWHMfinder2 

% Constants
linewidthTK = 2;
fconv = 128;

% Default values
llinestrdef = cell(1,nlline);
llinestrdef{1} = 'water';
for i = 2:nlline
    llinestrdef{i} = ['pool ',num2str(i-1)];
end
Albrdef = 0.5;          % lower bound ratio of the initial value
Aubrdef = 2;            % upper bound ratio of the initial value 
T2lbrdef = 0.5;
T2ubrdef = 2;
dfbrdef = 12.8;         % 0.1 ppm
fitoptionsdef = optimset('Display','notify','MaxFunEvals',100000,...
    'MaxIter',100000,'TolFun',1e-8,'TolX',1e-8,...
    'Algorithm','trust-region-reflective');
flag_shiftcorr_def = 1; % Use the min Z-spec value to replace the initial shift
flag_updateT2_def = 0;  % Update the T2 initial values for CEST pools
                        % or use fitted water T2 for all CEST pools 
flag_fig_def = 0;       % whether to show fitting result figure
unitdef = 'Hz';                     

% Inputs
p = inputParser;
addRequired(p,'zspec');
addRequired(p,'off');
addRequired(p,'nlline');
addRequired(p,'c0');
addRequired(p,'offwaterub');
addParameter(p,'llinestr',llinestrdef);
addParameter(p,'Albr',Albrdef);
addParameter(p,'Aubr',Aubrdef);
addParameter(p,'T2lbr',T2lbrdef);
addParameter(p,'T2ubr',T2ubrdef);
addParameter(p,'dfbr',dfbrdef);
addParameter(p,'fitoptions',fitoptionsdef);
addParameter(p,'shiftcorr',flag_shiftcorr_def);
addParameter(p,'updatewidth',flag_updateT2_def);
addParameter(p,'showfig',flag_fig_def);
addParameter(p,'unit',unitdef);

parse(p,zspec,off,nlline,c0,offwaterub,varargin{:});
llinestr = p.Results.llinestr;
Albr = p.Results.Albr;
Aubr = p.Results.Aubr;
T2lbr = p.Results.T2lbr;
T2ubr = p.Results.T2ubr;
dfbr = p.Results.dfbr;
fitoptions = p.Results.fitoptions;
flag_shiftcorr = p.Results.shiftcorr;
flag_updateT2 = p.Results.updatewidth;
flag_fig = p.Results.showfig;
unit = p.Results.unit;

% -------------------------------------------------------------------------
% Inputs
flag_debug = 0;
% -------------------------------------------------------------------------

% Initialization
noff = length(off);
if strcmp(unit,'ppm') && (dfbr == dfbrdef)
    dfbr = dfbr/fconv;
end
lb = zeros(1,3*nlline);
ub = zeros(1,3*nlline);
ylline = zeros(nlline,noff);

% update initial shift and water line fitting boundary
if flag_shiftcorr
    [~,ind2] = min(zspec);
    shift = off(ind2);
    for i = 1:nlline
        c0(3*i) = c0(3*i)+shift;
    end
    offwaterub = offwaterub+shift;
end
indwaterub = findc1D(off,offwaterub);

% update lower and upper bounds for shift
for i = 1:nlline
    lb(3*i) = c0(3*i)-dfbr;
    ub(3*i) = c0(3*i)+dfbr;
end

% fit the water line
% update initial amplitude
ind = findc1D(off,c0(3));
c0(1) = 1-zspec(ind);

% update lower and upper bounds for amplitude
lb(1) = Albr*c0(1);
ub(1) = Aubr*c0(1);

% update initial width
[FWHM,ind1,ind2] = FWHMfinder2(findc1D(off,c0(3)),(max(zspec)-min(zspec))/2,...
    off,max(zspec)-zspec);
c0(2) = FWHM;

% update lower and upper bounds for width
lb(2) = T2lbr*c0(2);
ub(2) = T2ubr*c0(2);

if flag_debug
    figure;
    hold all;
    plot(off,zspec);
    plot(off([ind2,ind1]),zspec([ind2,ind1]),'*');
    hold off;
    xlabel('Saturation frequency');
    set(gca,'XDir','reverse');
    ylabel('Signal (au)');
    title(['Find initial width for water line fitting: ',num2str(c0(2))]);
end

[cfbs,resnorm,residual,exitflag,output] = lsqcurvefit(@zspecfit2,...
    c0(1:3),off(1:indwaterub),1-zspec(1:indwaterub),lb(1:3),ub(1:3),...
    fitoptions,1);
fitparamdisp(1,llinestr{1},c0(1:3),lb(1:3),ub(1:3),cfbs);
yfbs = zspecfit2(cfbs,off,1);
res = 1-yfbs-zspec;

if flag_debug
    figure;
    hold all;
    plot(off,zspec,'bo','LineStyle','none');
    plot(off(1:indwaterub),zspec(1:indwaterub),'ro','LineStyle','none');
    plot(off,1-yfbs,'-');
    plot(off,res,'k-');
    hold off;
    legend('Measured not used in fitting','Measured used in fitting',...
        'Fitted','Residual');
    xlabel('Saturation frequency');
    set(gca,'XDir','reverse');
    ylabel('Signal (au)');
    title('water line fitting');
end

% fit all lines
% update initial amplitude
for i = 2:nlline
    ind = findc1D(off,c0(3*i));
    c0(3*(i-1)+1) = res(ind);
end

% update lower and upper bounds for amplitude
for i = 2:nlline
    lb(3*(i-1)+1) = Albr*c0(3*(i-1)+1);
    ub(3*(i-1)+1) = Aubr*c0(3*(i-1)+1);
end

if flag_updateT2
    % update initial width
    ind1 = zeros(1,nlline-1);
    ind2 = zeros(1,nlline-1);
    for i = 2:nlline
        [FWHM,ind1(i-1),ind2(i-1)] = FWHMfinder2(findc1D(off,c0(3*i)),c0(3*(i-1)+1)/2,off,res);
        c0(3*(i-1)+2) = FWHM;
    end
    
    if flag_debug
        figure;
        hold all;
        plot(off,res);
        for i = 1:(nlline-1)
            plot(off([ind1(i),ind2(i)]),res([ind1(i),ind2(i)]),'*');
        end
        hold off;
        xlabel('Saturation frequency');
        set(gca,'XDir','reverse');
        ylabel('Residual (au)');
        title('Find initial width for CEST line fitting');
    end
else
    % use water fitted width for all CEST pools
    for i = 2:nlline
        c0(3*(i-1)+2) = c0(2);
    end
end

% update lower and upper bounds for width
for i = 2:nlline
    lb(3*(i-1)+2) = T2lbr*c0(3*(i-1)+2);
    ub(3*(i-1)+2) = T2ubr*c0(3*(i-1)+2);
end

[cf,resnorm,residual,exitflag,output] = lsqcurvefit(@zspecfit2,...
    c0,off,1-zspec,lb,ub,fitoptions,nlline);
for i = 1:nlline
    fitparamdisp(i,llinestr{i},c0((3*(i-1)+1):(3*i)),...
        lb((3*(i-1)+1):(3*i)),ub((3*(i-1)+1):(3*i)),cf((3*(i-1)+1):(3*i)));
end
yf = zspecfit2(cf,off,nlline);
for i = 1:nlline
    ylline(i,:) = zspecfit2(cf((3*(i-1)+1):(3*i)),off,1);
end
res = 1-yf-zspec;

if flag_fig
    figure;
    hold all;
    plot(off,zspec,'bo','LineStyle','none');
    plot(off,1-yf,'r-','LineWidth',linewidthTK);
    plot(off,res,'k-','LineWidth',linewidthTK);
    for i = 1:nlline
        plot(off,1-squeeze(ylline(i,:)),'-');
    end
    hold off;
    legend([{'Measured'},{'Fitted'},{'Residual'},llinestr(:)']);
    xlabel(['Saturation frequency (',unit,')']);
    set(gca,'XDir','reverse');
    ylabel('Signal (au)');
    title(['Lorentzian fitting (',num2str(nlline),' lines)']);
end
end

% =========================================================================
% 20190228 SZ: 1st version.
% 20190423 SZ: make the main Lorentzian Fitting code a function.
% 20190423 SZ: update initial width for all pools.
% 20190523 SZ: make whether to update initial width for CEST pools optional
%              if not, use fitted water width for all CEST pools
% 20190523 SZ: make the shift boundary an absolute value instead of a ratio
% 20190523 SZ: add saturation frequency unit.
% 20190524 SZ: shift water line fitting boundary when shift corr is on.
% =========================================================================