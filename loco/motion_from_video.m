function motion_yn = motion_from_video(avi_full)   
    %determine where the mouse was moving or not by clustering the difference
    %btw pixels in the neighboring frames.
    
    % read video parameters
    v = VideoReader(avi_full);
    imw=v.Width;
    imh=v.Height;
    alldimsum = imw*imh;
    
    nrames =ceil(v.FrameRate*(v.Duration+1));    
    im=readFrame(v);
    motion_btw_frames = zeros(nrames-1,1);
    i=1;
    while hasFrame(v)
        im1=readFrame(v);
        dim = im1(:,:,1)-im(:,:,1);
        motion_btw_frames(i)=sum(dim(:))/alldimsum;
        im=im1;
        i=i+1;
    end
    nfr=i-1;
    motion_btw_frames=motion_btw_frames(1:nfr);
    
    %k-means kluster into two:
    [motion_inds, motc] = kmeans(motion_btw_frames,2);
    %set to zero indices with lower centroid values
    motion_yn=zeros(size(motion_inds));
    [~,mot_yes_ind]=max(motc);
    motion_yn(motion_inds==mot_yes_ind)=1;
    %smooth the results
    motion_yn = movmedian(motion_yn,31);
    motion_yn(motion_yn>0)=1;
     
%     %%create a new video with the idicator whether the mouse is moving 
%     avi_new = 'motion_yn.avi';
%     v = VideoReader(avi_full);    
%     vout = VideoWriter(avi_new);
%     open(vout);
%     i=1;
%     while i<=nfr
%         im=readFrame(v);
%         
%         if motion_yn(i)==1
%             im(50:75,50:75,2)=255;
%             im(50:75,50:75,1)=0;
%             im(50:75,50:75,3)=0;
%         else
%             im(50:75,50:75,:)=0;
%         end
%         
%         writeVideo(vout,im);
%         i=i+1;
%     end
%     close(vout);    
   
end