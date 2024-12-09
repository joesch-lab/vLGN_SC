function [bar_loc, barwidth_pixels] = vertical_random_bar_optointervals_reconstruct(logfile)
    %reconstruct the location of the random vertical bar using params in the
    %logfile
    load(logfile);

    barwidth_pixels=round(barwidth_deg/240*608);

    tex_width=608;
    barwidth_pixels_2=barwidth_pixels/2;

    first_rand_index=ceil(hstart*tex_width);
    last_rand_index=min(tex_width, hend*tex_width)-barwidth_pixels;

    rng(rseed);
    nframes=frame_count;
    bar_loc = zeros(1,nframes);
    frame_count=0;
    while frame_count < nframes      
        %get a new bar location
        if mod(frame_count,pattnfr)==0        
            %get a new random source rect
            randleft=randi([first_rand_index,last_rand_index])-1;           
        end
        bar_loc(frame_count+1) = round(tex_width - randleft-barwidth_pixels_2);
        % increment the frame counter for the next loop            
        frame_count=frame_count+1;      
    end
end
