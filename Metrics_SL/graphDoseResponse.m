function [ID, HillC_nfkb, HillC_ktr, halfA_conc_nfkb, halfA_conc_ktr, equ_nfkb, equ_ktr, doses] = graphDoseResponse(IDs, varargin)

%20201209 graphs % responders for NFkB and KTR, with Hill coeff
%todo add functionality to plot/determine Hill coeff on any metric of dose
%response

%for any number of experiment IDs
p = inputParser;
addRequired(p,'IDs'); %vector containing IDs to be plotted
%parameters to be passed to metric function
expectedFlags = {'on','off'};
valid_conv = @(x) assert(isnumeric(x)&&(x>=0)&&(length(x)==1),...
    'Parameter must be single integer >= 0'); %checks whether parameters below are single integers
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));%checks whether optional name-value argument matches on or off %checks if x matches expectedFlags
addParameter(p,'MinLifetime',107, @isnumeric); %allows adjustment of minimum lifetime (?)
addParameter(p,'MinSize',90, valid_conv); %allows adjustment of minimum size (?)
addParameter(p,'TrimFrame',107, @isnumeric);
addParameter (p, 'OnThreshNFkB', 3, @isnumeric); %sigma threshold for determining responders
addParameter (p, 'GraphLimitsNFkB',[-0.25 7],@isnumeric);
addParameter (p, 'OnThreshKTR', 3, @isnumeric); %sigma threshold for determining responders
addParameter (p, 'GraphLimitsKTR',[-0.02,0.35],@isnumeric);
addParameter(p, 'StimulationTimePoint', 13, @isnumeric)


parse(p,IDs, varargin{:})

n = numel(IDs);
ID(n).metrics = [];
ID(n).graph = [];
ID(n).info = [];


for i= 1:n
    [ID(i).metrics,~,ID(i).graph,ID(i).info,~] = nfkb_ktr_ratio_metrics(IDs(i), 'MinLifetime',p.Results.MinLifetime,...
                            'OnThreshNFkB',p.Results.OnThreshNFkB,'OnThreshKTR',p.Results.OnThreshKTR,...
                            'MinSize', p.Results.MinSize,'Verbose', ... 
                            p.Results.Verbose, 'TrimFrame', p.Results.TrimFrame, 'StimulationTimePoint', p.Results.StimulationTimePoint);
end



doses(n).dose = [];
data_for_resp_nfkb = nan(1,n);
data_for_resp_ktr = nan(1,n);
for i = 1:n
    data_for_resp_nfkb(i) = ID(i).metrics.responders_fraction_nfkb;
    data_for_resp_ktr(i) = ID(i).metrics.responders_fraction_ktr;
    doses(i).dose = str2double(ID(i).info.dose{:});
    %doses = [doses, str2num(ID(i).info.dose{:})];
end

f_nfkb = fit([doses.dose]', data_for_resp_nfkb', '(maxA*x.^HillC)./(halfA_conc.^HillC + x.^HillC) + intersect','Upper',[10,max([doses.dose]),0.9,1], 'Lower', [0.1, min([doses.dose]),0, 0] );
f_ktr = fit([doses.dose]', data_for_resp_ktr', '(maxA*x.^HillC)./(halfA_conc.^HillC + x.^HillC) + intersect','Upper',[10,max([doses.dose]),0.9,1], 'Lower', [0.1, min([doses.dose]),0, 0] );
x_for_fit = logspace(log10(min([doses.dose])), log10(max([doses.dose])), 100);


%xy plot NFkB, KTR responders

figure
s1 = scatter([doses.dose],data_for_resp_nfkb, 'o');
s1.MarkerEdgeColor = [0.9290 0.6940 0.1250];
s1.MarkerFaceColor = [0.9290 0.6940 0.1250];

hold on
s2 = scatter([doses.dose],data_for_resp_ktr, '*');
s2.MarkerEdgeColor = [0, 0.4470, 0.7410];
s2.MarkerFaceColor = [0, 0.4470, 0.7410];

p1= plot(x_for_fit,f_nfkb(x_for_fit), 'Color', [0.9290 0.6940 0.1250]);
p2 = plot(x_for_fit,f_ktr(x_for_fit), 'Color', [0, 0.4470, 0.7410]);

ylim([0 1])
xlim([doses(1).dose, doses(n).dose])
set(gca, 'XTick', [doses.dose],'XTickLabels', [doses.dose]);
set(gca, 'xscale', 'log')
title({'Responding cells[%]',num2str(IDs)})
legend({'NFkB responders', 'KTR responders', 'NFkB Hill curve', 'KTR Hill curve'}, 'Location', 'northeastoutside')
ylabel('Responder Fraction')
xlabel('Dose')

equ_nfkb = ['NFkB: (' num2str(f_nfkb.maxA) '*x.^' num2str(f_nfkb.HillC) ')./(' num2str(f_nfkb.halfA_conc) '.^' num2str(f_nfkb.HillC) '+ x.^' num2str(f_nfkb.HillC) ') +' num2str(f_nfkb.intersect)];
equ_ktr = ['KTR: (' num2str(f_ktr.maxA) '*x.^' num2str(f_ktr.HillC) ')./(' num2str(f_ktr.halfA_conc) '.^' num2str(f_ktr.HillC) '+ x.^' num2str(f_ktr.HillC) ') +' num2str(f_ktr.intersect)];
text(min([doses.dose]),0.09,equ_nfkb, 'FontSize', 8, 'Interpreter', 'none');
text(min([doses.dose]),0.03,equ_ktr , 'FontSize', 8, 'Interpreter', 'none');

hold off

HillC_nfkb = f_nfkb.HillC;
halfA_conc_nfkb = f_nfkb.halfA_conc;


HillC_ktr = f_ktr.HillC;
halfA_conc_ktr= f_ktr.halfA_conc;



