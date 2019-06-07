function yf = zspecfit2(c,off,n)
% =========================================================================
% Lorentzian fitting
% =========================================================================

% Check validity
if length(c) ~= 3*n
    error('ERROR: Number of coefficients does not match number of Lorentzian lines.');
end
% End check validity

yf = 0;
for i = 1:n
    eval(['A', num2str(i),' = c(',num2str(3*(i-1)+1),');']);
    eval(['T2',num2str(i),' = c(',num2str(3*(i-1)+2),');']);
    eval(['df',num2str(i),' = c(',num2str(3*i),');']);
    yf = eval(['yf+lorentzian2(A',num2str(i),',T2',num2str(i),...
        ',df',num2str(i),',off);']);
end

end
% =========================================================================
% 20190228 SZ: 1st version. Take from zspecfit() - 2 fit params per line.
% =========================================================================