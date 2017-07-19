%% ---------------------- Preparation_MRCAT_MRair -------------------------
%
% Goal: Preparation of MRCAT with air as detected on CT and on MR.
% The procedure is described in the paper.
% Note that the matching of the body contour as described in the paper
% is not here presented.
% This method is provided withouht any warraty, and is presented for
% reproducibility purposes. The code cannot be used for medical purposed
% without any further check. No responsability can be related to the
% developers from any misuse of the code.
% Reproduction and modification of the code is possible upon reference to
% the original code and the published paper.
% The code has been developed by Matteo Maspero in 2017 at UMC Utrecht.
% For info feel free to contact:
% m.maspero@umcutrecht.nl/matteo.maspero.it@gmail.com

% Brief description of the used functions:
% view3Dgui     -> Visualization 3D with window levelling allowed
% Body_CT       -> calcualtion of the body mask (usable both CT and MRI)

%% General Initialization

clc; clear;
close all;
% closereq is needed in case view3dgui is used as a viewer. This is a bug
% of view3dgui (which can be found online at https://nl.mathworks.com/
% matlabcentral/fileexchange/31370-3d-image-slice-viewer?
% focused=5189889&tab=function)
closereq;closereq;closereq;closereq;closereq;closereq
closereq;closereq;closereq;closereq;closereq;closereq
clc;

% General Settings for the uniformation of the figures font/fontsize
Size_ftn=20; Name_ftn='Times New Roman';
set(0,'DefaultAxesFontSize',Size_ftn,'DefaultAxesFontName',Name_ftn,...
    'DefaultTextFontname','Times New Roman');

%------------------------ Path to include with utils-----------------------
addpath(genpath('./utils/'));

%% Settings Preliminary Params

VisualVol=1234;  % 1234 to Visualize the volumes

%Selection of the patients
Pti=[3;4;7;8;9];

%% Pseudo-HU utilized in the pseudo-CT
% Huonsfield Unit used in MRCAT v257 (Philips, Vantaa, Finland), check
% also Maspero, Matteo et al., Quantification of confounding factors in
% MRI-based dose calculations as applied to prostate IMRT, 2017, Phys. Med.
% Biol. 62(3): 948-965. DOI:10.1088/1361-6560/aa4fe7

Air=-968;           %(-inf;-210)
Adip_tis=-86;       %[-210;-20)
Soft_tis=42;        %[-20;120)
Marrow=198;         %[120;555)
Cort_Bon=949;       %[555;inf) -> if not marrow=> [200;inf)

% The following values have been found after iterative optimization as
% described in the paper
Air_int=-600;
Marrow_new=265;
Cort_Bon_new=1250;

% Threshold used to detect the Air and the Body Contour on CT
Air_Thre=-200;

%% open Info MR (previously extracted from Dicom)6592.0
load('./MAt/DataMR.mat')

tiic=tic;
% Loop over the patient
for ssssl=5%1:numel(Pti)
    
    tStart = tic;
    Pt=Pti(ssssl);
    disp('-------------------------------');
    disp(['Processing Pt',num2str(Pt,'%04u')])
    
    % Open the selected patient
    load(['MAt/P',num2str(Pt,'%i'),'_CT.mat'])
    load(['MAt/P',num2str(Pt,'%i'),'_MRCAT.mat'])
    load(['MAt/P',num2str(Pt,'%i'),'_MRCAT_F.mat'])
    load(['MAt/P',num2str(Pt,'%i'),'_MRCAT_W.mat'])
    load(['MAt/P',num2str(Pt,'%i'),'_MRCAT_IP.mat'])
        
    %view3dgui(CT,VoxCT)
    if VisualVol==1234
        view3dgui(CT,DataMR.VoxelSize,'transverse')
        view3dgui(MRcat,DataMR.VoxelSize,'transverse')
    end
    
    %% Find Body massk on CT and MR
    
    tic
    [ BodyCT,BodyErode_CT] = Body_CT(CT,Air_Thre);
    [ BodyMR,BodyErode_MR] = Body_CT(MRcat,Air_Thre);
    
    disp(['Body mask calculated in ',num2str(toc,'%4.1f'),' s'])
    CT_Body=CT;
    CT_Body(logical(~BodyCT))=Air;
    
    %% Source MR Image Normalisation
    
    % Generate normalised W/F fraction
    WaterFraction=(MRcatW./(MRcatW+MRcatF));
    FatFraction=(MRcatF./(MRcatW+MRcatF));
    % Normalisation on In-Phase image
    IP=MRcatIP/max(MRcatIP(:));
    
    %% Generation of a binary mask for the air as detected on the source MRI
    
    % First gross step: air is what has low intensity in water and fat
    AirMR=BodyErode_MR&(MRcatIP/(max(MRcatIP(:)))<0.1);
    AirMR(WaterFraction>0.75)=0;
    AirMR(:,:,[1:8 end-2:end])=0;
    
    %Dilated Bone Mask from MRCAT
    Bones=BodyErode_MR;
    Bones(MRcat<100)=0;
    BonesDil=zeros(size(Bones));
    for ii=1:size(Bones,3)
        I=squeeze(Bones(:,:,ii));
        BonesDil(:,:,ii)= imdilate(I,ones(9,9));
    end
    for ii=1:size(Bones,2)
        I=squeeze(BonesDil(:,ii,:));
        BonesDil(:,ii,:)= imdilate(I,ones(5,5));
    end
    AirMR(BonesDil==1)=0;
    
    % Remove high water content
    AirMR(WaterFraction>0.85)=0;
    
    %Remove small components with volume < 75 ml
    AirMR1=AirMR;
    ml=75;
    for ii=1:size(AirMR,3)
        I=bwareaopen(squeeze(AirMR(:,:,ii)),...
            ceil(ml/prod(DataMR.VoxelSize)));
        AirMR1(:,:,ii)= I;
    end
    
    %% Assignment of new pseudo-HU in the bones    \
    
    MRcat_M265_C1250=MRcat;
    MRcat_M265_C1250(MRcat==Marrow)=Marrow_new;
    MRcat_M265_C1250(MRcat==Cort_Bon)=Cort_Bon_new;
    
    % MRCAT with air and corrected inside
    MRcatAir=MRcat_M265_C1250;
    MRcatAir(AirMR1&BodyMR)=Air_int;
    
    disp(['Total volume of Air detected on MR = ',...
        '                                    ',...
        num2str(numel(MRcat(AirMR1==1&BodyMR))*...
        prod(DataMR.VoxelSize)/1000,'%5.1f'),' ml'])
    disp(['Total volume of Air detected on MR within',...
        ' intersection of body contours= ',...
        num2str(numel(AirMR1(AirMR1==1&BodyMR&BodyCT))*...
        prod(DataMR.VoxelSize)/1000,'%5.1f'),' ml'])
    
    
    if VisualVol==1234
        view3dgui(MRcatAir,DataMR.VoxelSize,'transverse')
    end
    
    %% Detect Air in CT
    
    AirinCT=(CT<Air_Thre)&(BodyCT);
    AirinCT1=zeros(size(AirinCT));
    
    for ii=1:size(AirinCT,3)
        I=bwareaopen(squeeze(AirinCT(:,:,ii)),ceil(ml/prod(DataMR.VoxelSize)));
        AirinCT1(:,:,ii)= I;
    end
    
    MRcat_AirCT=MRcat_M265_C1250;
    MRcat_AirCT(AirinCT1&BodyMR)=Air_int;
    
    disp(['Total volume of Air detected on CT = ',...
        '                                     ',...
        num2str(numel(MRcat(AirinCT1&BodyMR==1))*...
        prod(DataMR.VoxelSize)/1000,'%5.1f'),' ml'])
    disp(['Total volume of Air detected on CT within',...
        ' intersection of body contours = ',...
        num2str(numel(CT(AirinCT1&BodyMR&BodyCT==1))*...
        prod(DataMR.VoxelSize)/1000,'%5.1f'),' ml'])
    
    if VisualVol==1234
        view3dgui(MRcat_AirCT,DataMR.VoxelSize,'transverse')
    end
    
    %% Difference MRCAT_MRair-CT
    
    if VisualVol==1234
        view3dgui(MRcatAir-CT,DataMR.VoxelSize,'transverse')
        % Tips: use a colormap (press the button 'c' and 'b') and set the
        % center to 0 and the width to 2000
    end
    
end

toc(tiic)