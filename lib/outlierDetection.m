function [ y , eta ] = outlierDetection ( x , t , Setup )
%
%         Outlier Detection algorithm
%           Every samples outside interval [Quartile1-C*1.5*IQR, Quartile3+C*1.5*IQR ] are removed,
%           It is computed using the previous "Ne" samples, and "eta" is depending on the moving
%           IQR of those "N" samples, and a constant "C".
%
%          The main difference with "A Robust Method for ECG-Based Estimation of the Respiratory
%                   Frequency During Stress Testing" by R. Bailï¿½n et al.
%           is that eta here is defined by the maximum and minimum limits of the boxplot (outliers), not the STD.
%
%
%   In:   x      	=	unevenly sampled signal that may contain outliers
%         t         =	time vector for "x" (seconds). must be same size as x
%
%  Setup: N         =	number of previous samples which are used in "median_y" and "eta_y"
%         C         =	constant which multiplies IQR to obtain "eta_y"
%         plotflag	=	if true -> function shows some figures [Default: false]
%
%   Out:  y     	=	unevenly sampled signal x, outliers replaced with NaN values
%         eta       =	Limits defined by: Quartile (1 or 3) +/- C*(1.5*IQR)
%                           eta(1) = Lower Limit
%                           eta(2) = Upper Limit
%
% Created by Pablo Armaï¿½ac in 2019
% Contact: <parmanac@unizar.es>
%
%
%
% AUN FALTA QUE LA BUSQUEDA DE MUESTRAS SEA EN SEGUNDOS Y NO EN MUESTRAS
% (EL TIEMPO 100 MUESTRAS ANTES DE LO ACTUAL NO ES LO MISMO QUE 100 SEGUNDOS ANTES pq son uneven!)
%

%% Parse Inputs

% Check Inputs
if nargin <= 2,   Setup = struct();                       end
if nargin <  2,   error('Not enough input arguments.');   end

if islogical(Setup)
    aux = Setup; clear Setup
    Setup.plotflag	=	aux;
end

% Default Values
if ~isfield(Setup,'C'),             Setup.C         =   1.5;      			end
if ~isfield(Setup,'N'),             Setup.N         =   20;			    	end
if ~isfield(Setup,'wdw_width'),     Setup.wdw_width	=	0;      			end
if ~isfield(Setup,'eta'),     		Setup.eta		=	NaN(length(x),2);	end
if ~isfield(Setup,'plotflag'),      Setup.plotflag	=	true;				end

% Get the variable names
data_names          = fieldnames(Setup);

% Assign and store the variables
for ii = 1:length(data_names)
    eval([ data_names{ii} ' = Setup.' data_names{ii} ';'])
end
clear Setup data_names ii


%% Check inputs
t           =	t	(:); %(~isnan(sig_in(:)));
x           =	x	(:); %(~isnan(sig_in(:)));


%% Estimate moving median
if wdw_width > 0,   moving_median	=	movmedian 	( x , [wdw_width wdw_width] , 'omitnan' );
else,               moving_median 	=	zeros 		( length(x) , 1 ); end

x_aux	=	x - moving_median;

%% Compute Threshold of outliers: if eta is not defined, compute eta
if sum(isnan(eta(:,1))) >= length(eta)

    if length(x_aux)>N/2

        % Initialize
        [ eta(1:N/2,1) 	   		, eta(1:N/2,2) 			]	=	computeEta( x_aux(1:N/2),		  C );
        [ eta(end-N/2:end,1) 	, eta(end-N/2:end,2) 	]	=	computeEta( x_aux((end-N/2:end)), C );

        % loop over signal
        for ii = N/2:length(x)-N/2

            % Select the last N samples
            ini = ceil(ii-N/2);  if ini<1,          ini=1;          end
            fin = ceil(ii+N/2);  if fin>length(x),	fin=length(x);	end
            idx = ini:fin;
            idx(isnan(x_aux(idx))) = [];

            if length(idx) < N/4
                continue

            else

                % Se comprueba si entra en rango
                if ( x_aux(ii)> eta(ii-1,2) ) || ( x_aux(ii) < eta(ii-1,1) )
                    x_aux(ii) = NaN;
                end

                % Si entra en rango se aÃ±ade para el computo del nuevo IQR:
                [ eta(ii,1) , eta(ii,2) ] = computeEta(x_aux(idx),C);

            end

        end

    else
        eta (:,1)	=	NaN;
        eta (:,2)	=	NaN;
    end

end


%% Delete outliers
y = x;

y ( x_aux >= eta(:,2) ) = NaN;
y ( x_aux <= eta(:,1) ) = NaN;


%% Figure:
if plotflag

    figure

    ax(1) = subplot(2,1,1); hold on; box on;
    plot(t, x, '*', 'Color' ,[0.8500    0.3250    0.0980], 'DisplayName', 'Outliers');
    plot(t, y, '*', 'Color' ,[     0    0.4470    0.7410],'HandleVisibility','off'); ylims = get(gca,'YLim');
    plot(t, moving_median, 'k','DisplayName','$\overline x$');
    plot(t, moving_median+eta(:,1), 'k--','DisplayName','$\eta$');plot(t, moving_median+eta(:,2), 'k--','HandleVisibility','off');
    set(gca,'ylim',ylims); ylabel('Output'); title('Outlier Detection, based on IQR'); lgd=legend; set(lgd,'Interpreter','latex');

    ax(2) = subplot(2,1,2); hold on; box on;
    plot(t, x - moving_median, '*', 'Color' ,[0.8500    0.3250    0.0980]);
    plot(t, y - moving_median, '*', 'Color' ,[     0    0.4470    0.7410]);
    plot(t, eta(:,1:2), 'k--');
    plot(t, eta(:,3:4), 'k:');
    xlabel('Time'); ylabel('w/o baseline');
    linkaxes(ax, 'x'); try xlim([t(1) t(end)]); catch,  end
    zoom on

end



end




% % Compute threshold Backward
% eta_inf     =	NaN( length(x) , 1 ); %eta_y will contain the threshold
% eta_sup     =	NaN( length(x) , 1 );
% % Initialization
% % [ eta_inf(end-N:end) , eta_sup(end-N:end) ] = computeEta(x_aux(end-N:end),C);
%
% for ii = length(x)-N/2:-1:N/2
%
%     % Select the last N samples
%     ini = ceil(ii-N/2);  if ini<1, ini=1; end
%     fin = ceil(ii+N/2);  if fin>length(x), fin=length(x); end
%     idx = ini:fin;
%     idx(isnan(x_aux(idx))) = [];
%
%     if length(idx) < N/2
%         continue
%
%     else
%
%         % Se comprueba si entra en rango
%         if ( x_aux(ii)> eta_sup(ii+1) ) || ( x_aux(ii) < eta_inf(ii+1) )
%             x_aux(ii) = NaN;
%         end
%
%         % Si entra en rango se aï¿½ade para el computo del nuevo IQR:
%         [ eta_inf(ii) , eta_sup(ii) ] = computeEta(x_aux(idx),C);
%
%     end
% end
%
% eta = [ eta eta_inf(:)  eta_sup(:) ];


%%
function [ eta_inf , eta_sup ] = computeEta (sig,C)

Quartiles	=	prctile	(sig,[25 50 75] );
IQR         =	Quartiles(3) - Quartiles(1);

%% Normal whiskers boxplot
eta_inf     =	Quartiles(1) - (1.5 * IQR) * C ;
eta_sup     =	Quartiles(3) + (1.5 * IQR) * C ;

% % MedCouple estimation: Adjusted whiskers for Skewed distributions
% MC          =	medcouple( sig(~isnan(sig)) );
%
% % Whiskers adjustment [Hubert&Vandervieren2008]
% if      MC >= 0
%     eta_inf = Quartiles(1) - ( 1.5*C * exp(-3.5*MC) * IQR ) ;
%     eta_sup = Quartiles(3) + ( 1.5*C * exp( 4*MC) * IQR ) ;
% elseif	MC <  0
%     eta_inf = Quartiles(1) - ( 1.5*C * exp(-4*MC) * IQR ) ;
%     eta_sup = Quartiles(3) + ( 1.5*C * exp( 3.5*MC) * IQR ) ;
% end

end



function mc = medcouple(x)
%
% 'medcouple' computes the medcouple measure, a robust measure of skewness
% for a skewed distribution. It takes into account cases where the
% observations are equal to the median of the series.
%
% Data in 'x' are organized so that columns are the time series and rows
% are the time intervals. All series contain the same number of
% observations.
%
% [mc] = medcouple(x) returns the following:
% mc    - vector with the medcouple measure of the data series
%
% Created by Francisco Augusto Alcaraz Garcia
%            alcaraz_garcia@yahoo.com
%
% References:
%
% 1) G. Brys; M. Hubert; A. Struyf (2004). A Robust Measure of Skewness.
% Journal of Computational and Graphical Statistics 13(4), 996-1017.


[n, c] = size(x);

[s_x,ix] = sort(x);
x_med = median(x);
z = s_x - repmat(x_med,n,1);

mc = zeros(1,c);
for w = 1:c
    [ip, jp] = find(z(:,w)>=0); % These are the positions in z of z+
    [im, jm] = find(z(:,w)<=0); % These are the positions in z of z-

    p = size(ip,1);
    q = size(im,1);

    [mi, mj] = ind2sub([p,q],1:p*q); % Positions of all combinations of z+
    % and z- as elements in a pxq matrix

    zp = z(ip,w); % z+ repeated to account for all cells in the matrix
    zm = z(im,w); % z- repeated to account for all cells in the matrix

    h = (zp(mi)+zm(mj))./(zp(mi)-zm(mj)); % same size as mi, mj

    [ipz,jpz]= find(zp==0); % row numbers of z+ = 0, i.e., x_{i} = median(x)
    [imz,jmz]= find(zm==0); % row numbers of z- = 0, i.e., x_{i} = median(x)
    piz = ismember(mi,ipz); % positions in mi where z+=0
    pjz = ismember(mj,imz); % positions in mi where z-=0
    zmember = piz+pjz; % same size as mi, mj
    pijz = find(zmember == 2); % positions where z+ = z- = 0, i.e., x_{i} =
    % = x_{j} = median(x)
    [indi,indj] = ind2sub([p,q],pijz); % pxq matrix position of the zero entries
    indi = indi - min(indi) + 1; % row position of the zero entries as if they
    % were in a separated matrix
    indj = indj - min(indj) + 1; % column position of the zero entries as if they
    % were in a separated matrix

    for i=1:size(pijz,2),
        if (indi(i) + indj(i) - 1) > size(find(z==0),1),
            h(pijz(i)) = 1;
        elseif (indi(i) + indj(i) - 1) < size(find(z==0),1),
            h(pijz(i)) = -1;
        else
            h(pijz(i)) = 0;
        end
    end

    mc(w) = median(h);
end


end


% function mc = medcouple(x)
% % 'medcouple' computes the medcouple measure, a robust measure of skewness for a skewed distribution.
% % It takes into account cases where the observations are equal to the median of the series.
% %
% %  MC Computed by the Naïve algorithm.
% % This algorithm could be optimized by exploiting the sorted nature of the medcouple matrix H
% %
% % mc    - vector with the medcouple measure of the data series
% %
% % References:
% %
% % [1] = G. Brys; M. Hubert; A. Struyf (November 2004). "A Robust Measure of
% %       Skewness". Journal of Computational and Graphical Statistics. 13
% %       (4): 996-1017.
% % [2] = M. Hubert; E. Vandervieren (2008). "An adjusted boxplot for skewed
% %       distributions". Computational Statistics and Data Analysis. 52 (12):
% %       5186-5201. doi:10.1016/j.csda.2007.11.008.
%
% s_x = sort(x);
% x_med = median(x);
% z = s_x - x_med;
%
% idx_plus = find( z>=0 ); % These are the positions in z of z+
% idx_minus = find( z<=0 ); % These are the positions in z of z-
%
% zp = z(idx_plus); % z+ repeated to account for all cells in the matrix
% zm = z(idx_minus); % z- repeated to account for all cells in the matrix
%
% p = length(idx_plus);
% q = length(idx_minus);
%
% [mi, mj] = ind2sub( [p,q] , 1:p*q ); % Positions of all combinations of z+ and z- as elements in a pxq matrix
%
% h = ( zp(mi)+zm(mj) ) ./ ( zp(mi)-zm(mj) ); % same size as mi, mj
%
% ipz = find(zp==0); % row numbers of z+ = 0, i.e., x_{i} = median(x)
% imz = find(zm==0); % row numbers of z- = 0, i.e., x_{i} = median(x)
%
% piz = ismember(mi,ipz); % positions in mi where z+ = 0
% pjz = ismember(mj,imz); % positions in mi where z- = 0
%
% zmember = piz+pjz; % same size as mi, mj
% pijz = find(zmember == 2); % positions where z+ = z- = 0, i.e., x_{i} = x_{j} = median(x)
%
% [indi,indj] = ind2sub([p,q],pijz); % pxq matrix position of the zero entries
%
% indi = indi - min(indi) + 1; % row position of the zero entries as if they were in a separated matrix
% indj = indj - min(indj) + 1; % column position of the zero entries as if they were in a separated matrix
%
% for ii = 1:length(pijz)
%     if     (indi(ii) + indj(ii) - 1) > size(find(z==0),1),  h(pijz(ii)) = 1;
%     elseif (indi(ii) + indj(ii) - 1) < size(find(z==0),1),  h(pijz(ii)) = -1;
%     else,                                                   h(pijz(ii)) = 0;
%     end
% end
%
% mc = median(h);
%
% end


% function MC = medcouple(data_in,plotflag)
% % FUNCTION MC = medcouple(data_in)
% %
% % In statistics, the medcouple is a non-parametric robust statistic that measures the skewness of a univariate distribution.[1]
% % It is defined as a scaled median difference of the left and right half of a distribution.
% % Its robustness makes it suitable for identifying outliers in adjusted boxplots.[2][3]
% % Ordinary box plots do not fare well with skew distributions, since they label the longer
% % unsymmetrical tails as outliers.
% % Using the medcouple, the whiskers of a boxplot can be adjusted for skew distributions and thus
% % have a more accurate identification of outliers for non-symmetrical distributions.
% %
% %  MC Computed by the Naï¿½ve algorithm. It could be optimized by
% %    exploiting the sorted nature of the medcouple matrix H
% %
% %
% %
% % data_in       : the data input
% % MC            : MedCouple value
% %
% %
% % calculations are based on formulas 2.1, 2.2, & 2.3 from [1].
% %
% % [1] = G. Brys; M. Hubert; A. Struyf (November 2004). "A Robust Measure of
% %       Skewness". Journal of Computational and Graphical Statistics. 13
% %       (4): 996ï¿½1017.
% % [2] = M. Hubert; E. Vandervieren (2008). "An adjusted boxplot for skewed
% %       distributions". Computational Statistics and Data Analysis. 52 (12):
% %       5186ï¿½5201. doi:10.1016/j.csda.2007.11.008.
% % [3] = Pearson, Ron (February 6, 2011). "Boxplots and Beyond ï¿½ Part II:
% %       Asymmetry". exploringdata.blogspot.ca. Retrieved April 6, 2015.
% %
%
% if nargin<2, plotflag	=	false; end
%
% med     =	nanmedian(data_in);
%
% % clean & sort data
% X       =	sort(data_in(~isnan(data_in)));
%
% Xli     =   X (X<=med); % X-: all values <= to median (Q2)
% Xgj     =   X (X>=med); % X+: all values >= to median (Q2)
%
% % compare all possible combinations
% Xi      =   repmat(Xli',length(Xgj),1); % make VERT matrix of X-' to match X+
% Xj      =   repmat(Xgj ,1,length(Xli)); % make HORZ matrix of X+  to match X-
% h       =   ( (Xj-med) - (med-Xi) ) ./ (Xj-Xi) ; % formula 2.2
% ties	=   find( isnan(h) );% only ties with median result in NaN values and are to be replaced
%
% if plotflag
%     figure(1234);clf;hold on;
%     title(sprintf('median=%0.2f; ties=%d \nNumEl=%d (%2.2f%%)',med(n),length(ties),numel(h),100*length(ties)/numel(h)))
%     colormap(jet); %blue/cyan <0; red/orange/yellow >0; green ==0;
%     h(ties)=0; % just for plotting (green); to make sure rotations are correct for kernel
%     tt=imagesc(double(h),[-1 1]); axis equal tight ij % proportions and orientation are critical
%     set(gca,'yticklabel',Xgi(get(gca,'ytick')));
%     set(gca,'xticklabel',Xlj(get(gca,'xtick')));
%     xlabel('values $\leq$ to median'); ylabel('values $\geq$ to median');
%     drawnow;
% end
%
% if any(ties)
%     % ties always make a square (in upper right of plot) so use flip_leftright             % -1 -1  0
%     % (triangle_lower - triangle_upper) to quickly create an appropriately sized kernel    % -1  0  1
%     k       =   ones( sqrt(length(ties)) );                                                %  0  1  1
%     kernel	=   fliplr( tril(k)-triu(k) );
%     h(ties)	=   kernel; % formula 2.3
%
%     if plotflag
%         set(tt,'CData',h)
%         drawnow;
%     end
%
% end
%
% MC = nanmedian( h(:) );% formula 2.1
%
% end