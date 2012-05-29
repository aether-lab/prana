function [MAD_ratio,MAD_ratio_hdr]=PTVval_meanandmedian(tracks,valprops)
%
% [MAD_ratio,MAD_ratio_hdr]=PTVval_meanandmedian(tracks,valprops);
%
% PROGRAM DESCRIPTION
% This function validates each particle track by statistically comparing it
% to neighboring tracks within a user defined search radius.  Once deemed
% an incorrect track, an position estimate is provided based on the motion
% of the surrounding particles
% 
% INPUTS
%   tracks - main  array of the matched particle pairs:
%       [X1 X2 Y1 Y2 Z1 Z2 d1 d2 I1 I2 p#1 p#2 match_probability]
%   valprops - control parameters for the validation of tracks
%       valprops.method - 'mean' or 'median' ; controls if the validation
%           is performed with the mean or median statistical properties
%       valprops.s_radius - pixel value to search for neighboring tracks
%       valprops.MAD_U - allowable ratio of mean/median absolute deviation
%           of each track vs. the avg. MAD of the neighboring tracks
%       valprops.MAD_V - same as above, but for the V-component of vel.
%
% OUTPUTS
%   MAD_ratio - estimated location of the ID'ed bad vectors for next
%       tracking step (2-D array)
%           [MAD_ratio_U; MAD_ratio_V; bad_U(0|1); bad_V(0|1); p#1; 
%           est_locx; est_locy]
%
%(v1) N. Cardwell - 10.31.2009

%get info and preallocate arrays
num_tracks=size(tracks,1);
MAD_ratio=zeros(num_tracks,7);

%main program loop for all tracks
for i=1:num_tracks

    %remove the track currently being evaluated (replace w/ NaN)
    trk_i=tracks(i,:);   tracks_i=tracks;   tracks_i(i,:)=NaN;
    
    %locate neighboring particles within user specified s_radius
    dist_i=sqrt( (tracks_i(:,1)-trk_i(1)).^2 + (tracks_i(:,3)-trk_i(3)).^2 + (tracks_i(:,5)-trk_i(5)).^2 );
    [r_i,c_i]=find(dist_i <= valprops.s_radius);

    %decision block - which method to use to ID and recalc "bad" tracks?
    switch lower(valprops.method)
        %use MEAN absoulte deviation
        case {'mean'}
            MAD_ratio_hdr={'MeanADratio Ui' 'MeanADratio Vi' 'bad vec U' 'bad vec V'...
                'part index' 'new est_locx' 'new est_locy'};

            %compute the mean of the neighboring vectors
            Ui=tracks_i(r_i,2) - tracks_i(r_i,1);
            Vi=tracks_i(r_i,4) - tracks_i(r_i,3);
            Mean_Ui=mean( Ui );   Mean_Vi=mean( Vi );

            %compute the mean absolute deviation of all neighboring vectors
            MeanAD_Ui=mad(Ui,0);   MeanAD_Vi=mad(Vi,0);

            %compute the mean deviation for trk_i
            trk_i_U = trk_i(2)-trk_i(1);   trk_i_V = trk_i(4)-trk_i(3);
            MeanAD_trk_iU = abs(trk_i_U - Mean_Ui);
            MeanAD_trk_iV = abs(trk_i_V - Mean_Vi);

            %compute the ratio of the MeanAD of the neighboring tracks
            %with the MeanAD of trk_i
            MAD_ratio(i,1) = MeanAD_trk_iU/MeanAD_Ui;
            MAD_ratio(i,2) = MeanAD_trk_iV/MeanAD_Vi;

            %determine if the MeanAD_ratio (either U or V) exceeds the user
            %specified limits (MAD_U and MAD_V)
            MAD_ratio(i,3) = MAD_ratio(i,1) > valprops.MAD_U;
            MAD_ratio(i,4) = MAD_ratio(i,2) > valprops.MAD_V;

            if ( MAD_ratio(i,3)==1 || MAD_ratio(i,4)==1 )==1
                %grab the particle number (related to SIZE1)
                MAD_ratio(i,5) = trk_i(11);
                
                %use the mean U and V velocity to get the new est_loc X & Y
                MAD_ratio(i,6) = trk_i(1) + Mean_Ui;
                MAD_ratio(i,7) = trk_i(3) + Mean_Vi;
            end

            %use MEDIAN absolute deviation
        case {'median'}
            MAD_ratio_hdr={'MedianADratio Ui' 'MedianADratio Vi' 'bad vec U' 'bad vec V'...
                'part index' 'new est_locx' 'new est_locy'};
            
            %compute the median of the neighboring vectors
            Ui=tracks_i(r_i,2) - tracks_i(r_i,1);
            Vi=tracks_i(r_i,4) - tracks_i(r_i,3);
            Median_Ui=median( Ui );   Median_Vi=median( Vi );

            %compute the median absolute deviation of all neighboring vectors
            MedianAD_Ui=mad(Ui,1);   MedianAD_Vi=mad(Vi,1);

            %compute the median deviation for trk_i
            trk_i_U = trk_i(2)-trk_i(1);   trk_i_V = trk_i(4)-trk_i(3);
            MedianAD_trk_iU = abs(trk_i_U - Median_Ui);
            MedianAD_trk_iV = abs(trk_i_V - Median_Vi);

            %compute the ratio of the MeanAD of the neighboring tracks
            %with the MeanAD of trk_i
            MAD_ratio(i,1) = MedianAD_trk_iU/MedianAD_Ui;
            MAD_ratio(i,2) = MedianAD_trk_iV/MedianAD_Vi;

            %determine if the MeanAD_ratio (either U or V) exceeds the user
            %specified limits (MAD_U and MAD_V)
            MAD_ratio(i,3) = MAD_ratio(i,1) > valprops.MAD_U;
            MAD_ratio(i,4) = MAD_ratio(i,2) > valprops.MAD_V;

            if ( MAD_ratio(i,3)==1 || MAD_ratio(i,4)==1 )==1
                %grab the particle number (related to SIZE1)
                MAD_ratio(i,5) = trk_i(11);
                
                %use the mean U and V velocity to get the new est_loc X & Y
                MAD_ratio(i,6) = trk_i(1) + Median_Ui;
                MAD_ratio(i,7) = trk_i(3) + Median_Vi;
            end

    end
end

% %Plotting code for all tracks and validated tracks
% im_bounds=[1 64 1 64 -1 1];
% figure;  scatter(tracks(:,1),tracks(:,3),'.b');  hold on
% scatter(tracks(:,2),tracks(:,4),'+r');
% for j=1:size(tracks,1)
%     line([tracks(j,1);tracks(j,2)],[tracks(j,3);tracks(j,4)]);
% end
% scatter(MAD_ratio(MAD_ratio(:,5)~=0,6),MAD_ratio(MAD_ratio(:,5)~=0,7),'xg');
% new_vecs=MAD_ratio(MAD_ratio(:,5)~=0,5:7);
% for j=1:nnz(MAD_ratio(:,5))
%     line([tracks(tracks(:,11)==new_vecs(j,1),1) ; new_vecs(j,2) ],...
%         [tracks(tracks(:,11)==new_vecs(j,1),3) ; new_vecs(j,3)],'Color','g');
% end
% legend({'image1','image2','im1-im2','image2:val','im1-im2:val'},...
%     'Location','NorthOutside','Orientation','horizontal')
% set(gca,'DataAspectRatio',[1 1 1]);
% set(gca,'YDir','reverse');
% axis(im_bounds(1:4))

%compress output array
MAD_ratio=MAD_ratio(MAD_ratio(:,5)~=0,:);

end