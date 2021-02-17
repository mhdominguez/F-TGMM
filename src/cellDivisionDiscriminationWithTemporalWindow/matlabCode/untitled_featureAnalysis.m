nb = 5;

figure;
for ff = 1:48

% % qqDS = xTestDS(:,impIdx(ff));
% % qqZC = xTestZC(:,impIdx(ff));
% % qqZS = xTestZS(:,impIdx(ff));


qqDS = xTestDS(yTestDS>0.5,impIdx(ff));
qqZC = xTestZC(yTestZC>0.5,impIdx(ff));
qqZS = xTestZS(yTestZS>0.5,impIdx(ff));

subplot(8,6,ff);
[p,u] = hist(qqDS, nb);
plot(u,p / sum(p));
[p,u] = hist(qqZC, nb);
hold on;
plot(u,p / sum(p),'r');
hold off;
[p,u] = hist(qqZS, nb);
hold on;
plot(u,p / sum(p),'g');
hold off;
%title(['ff = ' num2str(ff)]);
title_ = [num2str(ff) '.-' featureName{impIdx(ff)}]; 
pp = strfind(title_,'_');
title_(pp) = ' ';
title_ = [title_(1:25) '\newline' title_(26:end)];

title(title_);
xlim([min([qqDS;qqZC;qqZS]) max([qqDS;qqZC;qqZS])]);

end