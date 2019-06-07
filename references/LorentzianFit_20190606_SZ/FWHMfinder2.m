function [FWHM,ind1,ind2] = FWHMfinder2(indpk,hpkv,x,y)
% =========================================================================
% Find the FWHM
% =========================================================================

% -------------------------------------------------------------------------
% Inputs
flag_debug = 0;
% -------------------------------------------------------------------------

np = length(y);             % number of points

diff = y-hpkv;

indzc = NaN.*ones(1,np);
cnt = 0;
for i = 2:np
    if diff(i-1)*diff(i) <= 0
        cnt = cnt+1;
        [~,indmin] = min([diff(i-1),diff(i)]);
        indzc(cnt) = i-2+indmin;
    end
end

[~,indzc2] = sort(abs(indzc-indpk));
ind1 = indzc(indzc2(1));
ind2 = indzc(indzc2(2));
FWHM = abs(x(ind1)-x(ind2));

if flag_debug
    figure;
    hold all;
    hplot = zeros(4,1);
    hplot(1) = plot(x,y,'k.-');
    hplot(2) = plot(x(indpk),y(indpk),'ko');
    hplot(3) = plot(x(ind1),y(ind1),'r*');
    hplot(4) = plot(x(ind2),y(ind2),'b*');
    hlinehpkv = line(get(gca,'xlim'),[hpkv hpkv]);
    set(hlinehpkv,'Color','k','LineStyle','--','LineWidth',0.75);
    xlabel('x');
    xlim([min(x) max(x)]);
    ylabel('y');
    legend(hplot(:),'data','peak','FWHM (1)','FWHM (2)');
end

end
% =========================================================================
% 20181126 SZ: 1st version.
% 20190426 SZ: add the case if the input pk is not the true pk or the line
%              has multiple pks.
% =========================================================================