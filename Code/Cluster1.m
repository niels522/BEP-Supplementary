function [ClustSize, diam, Loc2particle,massCenter]=Cluster1(FileName, Bandwidth, MinPts, MaxDiam, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CLUSTER1: perform clustering of (X,Y,T) coordinates with Mean-Shift
%%%%%% algorithm using a SINGLE channel to identify clusters having
%%%%%% particle-like shape with user defined selection parameters.


%%%%%% REQUIRED INPUTS:
%%%%%% FileName: name txt file with XYT coords, in nm and frame number
%%%%%% bandwidth: parameter for clustering, in nm
%%%%%% MinPts: minimum number of localizations in a cluster
%%%%%% MaxDiam: maximum diameter (longest axis), in nm

%%%%%% OPTIONAL INPUTS:
%%%%%% Elong: max ellipse elongation allowed, default=2.0 
%%%%%% ScaleFactor: scale factor in ellipse fit, default=1.0
%%%%%% MinClustDist: minimum distance between clusters (in nm), default=300 
%%%%%% AggrDist: minimum distance between clusters to be non-aggregates (in
%%%%%% nm), default like MinClustDist
%%%%%% FracThreshold: fraction of cluster localizations within fitted
%%%%%% sphere

%%%%%% OUTPUTS:
%%%%%% Loc2particle: cell array with all info for each selected cluster
%%%%%% ClusSize: number of localizations in each selected cluster
%%%%%% diam: diameter of each selected cluster (in nm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% beginning of function...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Parse input data =====================================================

p = inputParser; %init parser object
validNum = @(x) isnumeric(x) && (x > 0); %define valid inputs: positive num
validChar = @(x) ischar(x);

%define defaults values for optional param:
defaultElong = 2.0; %elongation in ellipse fit
defaultScaleFactor= 1.0; %scale factor in ellipse fit
defaultMinClustDist = 300; % min distance between clusters
defaultAggrDist = defaultMinClustDist; % min distance between clusters to be non-aggregate, by default like MinClustDist
defaultFracThreshold = 0.9; %fraction of cluster localizations within fitted sphere

%define required and optional input parameters:
addRequired(p,'FileName',validChar);
addRequired(p,'Bandwidth',validNum); 
addRequired(p,'MinPts',validNum);
addRequired(p,'MaxDiam', validNum);

addParameter(p,'Elong', defaultElong, validNum);
addParameter(p,'ScaleFactor', defaultScaleFactor, validNum);
addParameter(p,'MinClustDist', defaultMinClustDist, validNum);
addParameter(p,'AggrDist', defaultAggrDist, validNum);
addParameter(p,'FracThreshold', defaultFracThreshold, validNum);

%read input values:
parse(p,FileName, Bandwidth, MinPts, MaxDiam,varargin{:});

%assign the parsed values:
Last = Inf; % Limit the number of points processed (set to Inf to process all)
MaxParticleElongation = p.Results.Elong; % max elongation allowed in ellipse fit
EllipseFitScaleFactor = p.Results.ScaleFactor; % scale factor in ellipse fit 
MinClustDst = p.Results.MinClustDist; % min distance of closest cluster to be considered isolated
distanceAggregates = p.Results.AggrDist; % min distance between clusters to be non-aggregate
FracThreshold = p.Results.FracThreshold; % fraction of cluster localizations within fitted sphere 

% % if not introduced, by default distance aggregate is set as MinClustDst:
% if (MinClustDst ~= defaultMinClustDist) && (distanceAggregates == defaultAggrDist)
%     distanceAggregates = defaultMinClustDist;
%     

%%% Read coordinates =====================================================

disp('Reading data...');

% read XYT coords in txt file:
Coordinates=importdata(FileName);
% %%%read txt with more complex file structure:
% % delimiterIn = ' ';
% % headerlinesIn = 1;
% % Coordinates = importdata(FileName,delimiterIn,headerlinesIn);

% assign X-Y-T coordinates:
XCoords647 = Coordinates(:,1); %X coords in first column
YCoords647 = Coordinates(:,2); %Y coords in second column
TCoords647 = Coordinates(:,3); %T coords in third column

% Plot raw coordinates data:
XCoords647 = XCoords647(1:min(Last,numel(XCoords647)));
YCoords647 = YCoords647(1:min(Last,numel(YCoords647)));
plot(XCoords647,YCoords647,'.r'); axis equal; hold on;

%%% Clustering =====================================================

% Clustering using mean-shift:
disp('Clustering...');
Pts647 = [XCoords647 YCoords647];
[clustCent647,~,clustMembsCell647] = MeanShiftCluster(Pts647.',Bandwidth); % .' transpose of the matrix
ClustSize647 = cellfun(@numel, clustMembsCell647);
NClust647 = numel(clustMembsCell647);
disp(strcat(['Found ' num2str(NClust647) ' clusters']));

% Flag isolated Nanoparticles
[D, ~] = pdist2(clustCent647',clustCent647','euclidean','Smallest',2); %calculate the 2 smallest distances (the smallest is 0, because comparing vector with itself)
if (size(D,1) > 1)
    isolated647 = (D(2,:) >= MinClustDst); % logic array selecting clusters isolated, needed for Ellipse-Fit
else
    isolated647 = true;
end

%%% Filter by localization number and elongation ==========================
 
disp('Filtering...');
 
IndexRightClust647=false(NClust647,1); % init logical array to flag valid clusters
RightClust647=1:1:NClust647; % intialization cluster indexes
massCenter=clustCent647; % initialization centers coords of clusters, to be filterd
r=0;
for i = 1 : NClust647 %loop to screen all clusters
    A = [XCoords647(cell2mat(clustMembsCell647(i)))'; YCoords647(cell2mat(clustMembsCell647(i)))' ]; %exctract XY coord for each cluster
    if(size(A,2) >= MinPts)   
        % Filter particle based on ellipse elongation & major axis length
        % Plot ellipse as reference
        ellipse_t = fit_ellipse_app(A(1,:),A(2,:),EllipseFitScaleFactor,MaxParticleElongation,MaxDiam, isolated647(i));           
        valid = ellipse_t.valid;
        if valid == 1
            IndexRightClust647(i)=true; % flag as true the valid cluster (by localiz number and elongation)
        end
    else
      disp(strcat( ['#' num2str(i) ' cluster has few localizations and has been excluded']));  
      r=r+1; %count discarded clusters
    end
end

%filter indexes and centers by localiz N and elongation
RightClust647=RightClust647(IndexRightClust647);
massCenter=massCenter(:,IndexRightClust647);
%ClustSize647=ClustSize647(IndexRightClust647); % size of valid clusters

disp(strcat( [num2str(NClust647-length(RightClust647)-r) ' clusters have been excluded due to unrealistic size or elongation']));

hold on
plot(massCenter(1,:), massCenter(2,:), 'xk', 'LineWidth',3, 'MarkerSize',10); 
% this is to visualize the label of the identified clusters

for i = 1 : size(massCenter,2)
txt1 = ['\leftarrow ' num2str(i)];
text(massCenter(1,i),massCenter(2,i),txt1)
end

%%% Filter particle too close /aggregates =================================
%%% Note: this filter works on pre-filtered clusters

distanceCenters=pdist2(massCenter',massCenter'); %euclidean dist between previously filtered cluster centers
[rows,~]=find((distanceCenters~=0) & (distanceCenters<distanceAggregates)); %index of clusters closer than threshold dist
AggregateMember= rows';
Aggregates=RightClust647(AggregateMember); %store the index of previously filtered clusters corresponding to aggregates
AggrMembr=clustMembsCell647(Aggregates); %for every aggregate which points are in it
RightClust647(AggregateMember)=[]; %filter index of clusters corresponding to aggregates 
NPMembs=clustMembsCell647(RightClust647); %for every selected cluster which points are in it 
NPSize647 = cellfun(@numel, NPMembs);  % for each selected cluster, counts the number of points in it

disp(strcat( [num2str(length(AggregateMember)) ' clusters have been excluded because were forming aggregates']));
disp(strcat(['Identified ' num2str(length(RightClust647)) ' nanoparticles candidates from ' num2str(NClust647) ' candidate clusters']));

% plot (overlay) the localizations of aggregated clusters
figure
hold on
for i=1: length(Aggregates)
    AggregateData=[XCoords647(cell2mat(AggrMembr(i))), YCoords647(cell2mat(AggrMembr(i)))];
    plot(AggregateData(:,1), AggregateData(:,2), 'xy')
end

%%% Size Check ============================================================
%this performs a circle-fitting of the clusters to calculate the radius
%comprising a fraction of the localization, NPs with unrealistic size
%are discarded

DataType='SolidSphere';
nclusters=length(RightClust647);
%SizeCheck=zeros(nclusters,1);
C=zeros(nclusters,2);
R=zeros(nclusters,1);
Rcheck=zeros(nclusters,1);
RCheckLow=5; % lower limit for radius value, set to 5nm
RCheckHigh=MaxDiam*2; % upper limit for radius value, set to 2times MaxDiam
discard=0;
IndexToRemove=true(nclusters,1); % logical array to remove index of cluster discarded

for i=1:nclusters %loop to screen all previously filtered clusters
    switch DataType
        case 'SolidSphere'
            ClusterData=[XCoords647(cell2mat(NPMembs(i))), YCoords647(cell2mat(NPMembs(i)))]; %extract XY coords of cluster
            %plot(ClusterData(:,1), ClusterData(:,2), 'xg')
            Cinitial=[mean(ClusterData(:,1)) mean(ClusterData(:,2))]; %mean coord as initial guess
            CovMat=cov(ClusterData); %calculate covariance
            Rinitial=1.5*mean([sqrt(CovMat(1,1)) sqrt(CovMat(2,2))]); %initial guess of radius (1.5x the stDev)
            
            % Now the location of the center of the smallest sphere 
            % encompassing X% of the datapoints is determined.
            CenterSampleSizeAz=10; CenterSampleSizeRad=4;
            CenterLocation=zeros(CenterSampleSizeAz+1,2,CenterSampleSizeRad);
            TrackLocNumber=zeros(CenterSampleSizeAz+1,CenterSampleSizeRad);
            
            for n=1:CenterSampleSizeRad %screen circles of increasing radius
                t2 = linspace(0,2*pi,CenterSampleSizeAz);
                XCurrentCircle=0.05*n*Rinitial*cos(t2)+Cinitial(1); %Xcoords of an approx circle with initial center and increasing radius
                YCurrentCircle=0.05*n*Rinitial*sin(t2)+Cinitial(2); %Ycoords of an approx circle with initial center and increasing radius
                for o1=1:CenterSampleSizeAz %screen angular directions 
                        CenterLocation(o1,:,n)=[XCurrentCircle(o1) YCurrentCircle(o1)];%coordinate of each (XY) point of the circle, different angular directions
                        LocDistCheck=find(((ClusterData(:,1)-XCurrentCircle(o1)).^2+(ClusterData(:,2)-YCurrentCircle(o1)).^2) < Rinitial^2); %for each point (radius and angular direct), see if better than intial radius
                        TrackLocNumber(o1,n)=length(LocDistCheck); % for each (radius and angular dir): how many local of cluster are closer if compared to initial rad
                end
            end
            
            RefLocCheck=find(((ClusterData(:,1)-Cinitial(1)).^2+(ClusterData(:,2)-Cinitial(2)).^2) < Rinitial^2);%find distances from initial guess center shorter than initial guess radius
            RefLocNumber=length(RefLocCheck); %a score of reference: how many many local of cluster within initial guess radius
            MaxLoc=max(max(TrackLocNumber)); %max score for each point (n-radius and o1-angular dir) previously screeened
            if MaxLoc > RefLocNumber %if better than initial guess
                [I,J]=find(TrackLocNumber==MaxLoc); %find the indexes corresponding to angular direction (o1) and radius (n) screened having max score 
                MaxLocMinRadTemp=[I J]; %store the indexes of the angular dir and radius
                MaxLocMinRadInd=(J==min(J)); % find(J==min(J)); %logical label the indexes (o1-angular dir and n-radius) with smallest radius (n)
                MaxLocMinRad=MaxLocMinRadTemp(MaxLocMinRadInd,:); %filter to select indexes (o1-angular dir and n-radius) with smallest radius (n)
                Cfinal=mean(CenterLocation(MaxLocMinRad(:,1),:,min(J)),1); %move the circle center in the new position
            else
                Cfinal=Cinitial;
            end
            
            RadiusVec=linspace(0,3*Rinitial,1500);
            RadiusNumLoc=zeros(length(RadiusVec),1);
            for p=1:length(RadiusVec) %more accuarte screen of radius
                LocDistCheckRad=find(((ClusterData(:,1)-Cfinal(1)).^2+(ClusterData(:,2)-Cfinal(2)).^2) < RadiusVec(p)^2); %indexes of cluster localization within radius
                RadiusNumLoc(p)=length(LocDistCheckRad);%for each screened radius  a score: number of localizations within it
            end
            
            RadTrack=1;
            while RadiusNumLoc(RadTrack) < FracThreshold * size(ClusterData,1) %this loop find which of the screened radius has FractThresh of cluster localiz
                RadTrack=RadTrack+1;
            end
            Rfinal=RadiusVec(RadTrack); %the final radius is set to the minimum comprisin FractThresh of cluster localiz 
            
            %fore each cluster, set new values of center and radius
            C(i,:)=Cfinal';
            Rcheck(i)=Rfinal;       
    end
    
    %a final check that radius falls within a realistic range
    if Rcheck(i) >= RCheckLow && Rcheck(i) <= RCheckHigh
        R(i)=Rfinal;
    else
        disp(['Cluster #' num2str(i) ' has been excluded from the analysis due to an unrealistic radius:' num2str(round(Rcheck(i))*2)])
        IndexToRemove(i) = false; % collect the indexes of the clusters to remove
        discard=discard+1;
    end
end
disp(strcat(['Identified ' num2str(nclusters-discard) ' valid nanoparticles from ' num2str(nclusters) ' candidate nanoparticles']));
NPselect= R~=0;
R=R(NPselect);
R=round(R);
C=[C(NPselect,1), C(NPselect,2)];
NPSize647 = NPSize647(IndexToRemove);  %remove the discarded cluster, after size check
NPMembs=NPMembs(IndexToRemove); %remove discarded clusters, after size check

%%% Generate output information ===========================================

disp('Creating output...');

figure
title({'black ellipse: circle fitting';'red dots: selected-NP localiz, black dots: other localiz'});grid on;axis equal;
hold on
plot(XCoords647,YCoords647,'.k'); axis equal; hold on;

% Build a Loc2particle cell array with selected X,Y,T coords for each particle:

Lgt=length(NPMembs);
Loc2particle = cell(Lgt,1); %initialize cell array, number of entries corresponding to numb of selected NPs

%fill the cell array: for each selected NP the X,Y,T coords of localizations 
for k = 1 : Lgt
    A = [XCoords647(cell2mat(NPMembs(k)))'; YCoords647(cell2mat(NPMembs(k)))'; TCoords647(cell2mat(NPMembs(k)))' ];
    Loc2particle{k}=A;
    plot(A(1,:), A(2,:), 'r.'); hold on; %plot selected localiz in red
end

%this is to plot the retrieved radius for selected NPs:
for m=1:length(R)
       t = linspace(0,2*pi,100);
    plot(R(m)*cos(t)+C(m,1),R(m)*sin(t)+C(m,2),'k','LineWidth', 1)
    hold on
%     axis image
    %axis([C(m,1)-500 C(m,1)+500 C(m,2)-500 C(m,2)+500])
end

% Statistics and export:
ClustSize=NPSize647; %how many MAIN Localizations are in each nanoparticle
% ClustSize = cellfun(@numel, Loc2particle); %how many MAIN Localizations are in each nanoparticle
% ClustSize=ClustSize/3; %IMPORTANT! Because count 1 localization as 3 (3 coords)
diam=R*2;

%histograms of local and diameter:
% figure(3);hist(ClustSize,30);grid on;title('Number of 647 localizations/nanoparticle');
% figure(4);hist(R*2,10);grid on;title('Nanoparticle diam');
% figure(5);scatter(diam,ClustSize,'filled');grid on;title('Diameter/localizations');xlabel('Diameter (nm)')