function gridcompare(DEM, varargin )
% 
%  GRIDCOMPARE, compares two grids (e.g. Elevation and Slope)
% 
%  Syntax:
%         gridcompare(GRIDobj)
%         gridcompare(GRIDobj,'name',grid)
% 
% 
% Description:
% 
%  GRIDCOMPARE, compares two grids (e.g. Elevation and Slope) to highlight
%  any transient signal that may exist within a drainage basin or GRIDobj.
%  The resulting plot is similar to Fig. 7 Gallen et al., (2011) or Fig. 6 
%  Morriss and Wegmann, (2017).

%  You can either provide a DEM and this function will calculate the
%  gradient or relief, or you can provide your own grids (up to 3).
% 
% Required Input:
% 
%   DEM (as a topotoolbox GRIDobj)
% 
% Optional Input:
% 
%     slope grid - optional, must match area for elevation grid, specify with
%         'slope' before grid input.
%     relief grid - optional, must match area for elevation grid, specify with
%         'relief' before grid input.
%         
% 
% Examples:
% 
% DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
% gridcompare(DEM)
% 
% DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
% G = gradient8(DEM,'degrees');
% gridcompare(DEM,'slope',G);
% 
% % Look at just the largest drainage basin in BT grid
% DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
% FD = FLOWobj(DEM);
% S = STREAMobj(FD,'minarea',1e3);
% S = klargestconncomps(S);
% [L out] = drainagebasins(FD,S);
% DEMc = crop(DEM,L==1,nan);
% gridcompare(DEMc);
% 
% 
% Author: Morriss, Matthew (Matthew.c.morriss[at]gmail.com)
% with significant input from Sean Gallen 
% Date: September 2, 2018

	% Parse Inputs
p = inputParser;
p.FunctionName = 'gridcompare';
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    
addParameter(p,'SLOPE',[],@(x) isa(x,'GRIDobj'));
addParameter(p,'RELIEF',[],@(x) isa(x,'GRIDobj'));

parse(p,DEM,varargin{:});
DEM=p.Results.DEM;
SLOPE = p.Results.SLOPE;
RELIEF = p.Results.RELIEF;

home;
if isempty(SLOPE)
    disp([ datestr(now) ' No slope provided, calculating  ']);
    SLOPE = gradient8(DEM,'degree');
    
    
end

if isempty(RELIEF)
    disp([ datestr(now) ' No relief provided, calculating  ']);
    RELIEF = localtopography(DEM,500,'N',4);
    
end

elev.Max = nanmax(DEM.Z(:));
elev.Min = nanmin(DEM.Z(:));
elev.binWidth = 50; % may change in future verions. (10 m)
 

binVector = elev.Min:elev.binWidth:elev.Max; 

%%% pre allocate vectors %%
binMpt = nan(1,length(binVector)-1); % bin midpoint
elFreq = nan(1, length(binVector)-1); % elevation frequency
meanS = nan(1, length(binVector)-1); % mean slope
stdS = nan(1,length(binVector)-1); %std deviation of the slope


for i = 2:length(binVector); %for loop to loop through and fill all vectors
    
    binMpt(i) = (binVector(i) + binVector(i-1))/2; %should be bin midpiont

    %elFreq(i) = length(demGrid(demGrid <= binVector(i) & demGrid > binVector(i-1))); %elevation freq.
    elFreq(i) = length(~isnan(DEM.Z(DEM.Z <= binVector(i) & DEM.Z > binVector(i-1))));

    meanS(i) = nanmean(SLOPE.Z(DEM.Z <= binVector(i) & DEM.Z > binVector(i-1))); %mean slop
    stdS(i) = nanstd(SLOPE.Z(DEM.Z <= binVector(i) & DEM.Z > binVector(i-1))); %std deviation slope
    
    meanR(i) = nanmean(RELIEF.Z(DEM.Z <= binVector(i) & DEM.Z > binVector(i-1))); %mean slop
    stdR(i) = nanstd(RELIEF.Z(DEM.Z <= binVector(i) & DEM.Z > binVector(i-1))); %std deviation slope
    
end

%%% PLOT THE RESULT
figure('units','normalized','outerposition',[0 0 1 1]) 
set(gca, 'Units','normalized');
set(gcf, 'color', [1 1 1]);
hold on;
    
    yyaxis left
    bar(binMpt,elFreq,'DisplayName','Hypsometry'); 
    
    title('Elevation v Slope'); 
    xlabel('Elevation (m)'); 
    ylabel('Elevation bins');
    set(gca,'FontSize',16);

    yyaxis right
    hold on;
        plot(binMpt, meanS,'DisplayName','Mean Slope')
        plot(binMpt, meanS-stdS,'b--','DisplayName','Std of slope');
        plot(binMpt, meanS+stdS,'b--');
        ylabel('Slope (°)');
        set(gca,'FontSize',16);

    hold off;


figure('units','normalized','outerposition',[0 0 1 1]) 
set(gca, 'Units','normalized');
set(gcf, 'color', [1 1 1]);
hold on;
    
    yyaxis left
    bar(binMpt,elFreq,'DisplayName','Hypsometry'); 
    title('Elevation v Relief'); 
    xlabel('Elevation (m)'); 
    ylabel('Elevation bins');
    set(gca,'FontSize',16);


    yyaxis right
    hold on;
        plot(binMpt, meanR,'DisplayName','Mean Relief')
        plot(binMpt, meanR-stdR,'b--','DisplayName','Std of Relief');
        plot(binMpt, meanR+stdR,'b--');
        ylabel('Relief (M)');
        set(gca,'FontSize',16);
        
        
    hold off;



end

