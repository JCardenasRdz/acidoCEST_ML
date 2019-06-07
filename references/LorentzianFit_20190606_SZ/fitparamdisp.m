function fitparamdisp(illine,llinelab,initval,lb,ub,fitval)
% =========================================================================
% Display the fitting params
% =========================================================================

fprintf('\nLorentzian lineshape #%d - %s:\n',illine,llinelab);
fprintf('Initial:\tA %3.2f,\twidth %3.2f,\tshift %3.2f\n',...
    initval(1),initval(2),initval(3));
fprintf('LBound:\t\tA %3.2f,\twidth %3.2f,\tshift %3.2f\n',...
    lb(1),lb(2),lb(3));
fprintf('UBound:\t\tA %3.2f,\twidth %3.2f,\tshift %3.2f\n',...
    ub(1),ub(2),ub(3));
fprintf('Fitted:\t\tA %3.2f,\twidth %3.2f,\tshift %3.2f\n',...
    fitval(1),fitval(2),fitval(3));

end
% =========================================================================
% 20190301 SZ: 1st version.
% =========================================================================