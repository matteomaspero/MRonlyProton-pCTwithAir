function [ Body_Mask_only,BodyErode] = Body_CT( data,thre )
%Body_Mask Body mask extracted using a hard thresholding and a filling, 
% after removal of small components. An additional mask eroded of 5 voxels
% is also produced.
%-----------------------------Nota Bene -----------------------------------
% This method is provided withouht any warraty, and is presented for
% reproducibility purposes. The code cannot be used for medical purposed
% without any further check. No responsability can be related to the
% developers from any misuse of the code.
% Reproduction and modification of the code is possible upon reference to
% the original code and the published paper.
% The code has been developed by Matteo Maspero in 2017 at UMC Utrecht.
% For info feel free to contact:
% m.maspero@umcutrecht.nl/matteo.maspero.it@gmail.com

dim=size(data);
Body_Mask_only=zeros(dim);
BodyErode=Body_Mask_only;
Body_Mask=data>thre;
 for ii=1:dim(3)
     I=squeeze(Body_Mask(:,:,ii));
     II=logical(bwareaopen(I,ceil(dim(1)*dim(2)*0.01)));
     III=imfill(II,'holes');
     Body_Mask_only(:,:,ii) = III;
     BodyErode(:,:,ii)= bwmorph(III,'erode',5);
end

