function aniid_all = get_animal_id_from_recordingname(recname_list)
%get id of the animal from the filename         

if ~iscell(recname_list) % only one recording is given
    recname_list_in=recname_list;
    recname_list={};
    recname_list{1}=recname_list_in;
end

aniid_all =[];
for i=1:length(recname_list)
      if iscell(recname_list(i))
          recname=cell2mat(recname_list(i));
      else
        recname=recname_list(i);
    end

    %old recording id is 1-2 digits after 'fiber' and 'ipsi'
    if contains(recname, 'fiber','IgnoreCase', true) && contains(recname, 'ips','IgnoreCase', true)
        m1=strfind(lower(recname), 'fiber')+length('fiber');
        m2=strfind(lower(recname), 'ips')-1;        
        aniid=str2num(recname(m1:m2));
    elseif contains(recname, '_T'+ digitsPattern(1,5)) %when animal is encoded between '_T' and '_'
        m1=strfind(recname, '_T')+length('_T');
        m2=strfind(recname(m1:end), '_'); m2=m1+(m2(1)-1)-1;       
        aniid=str2num(recname(m1:m2));
    else %animal id is 3-4 numbers surrounded by -/_/MJ
        recnamesplit=strsplit(recname,{'-','_', 'MJ'});
        celllen=cellfun(@length, recnamesplit);
        celldigits=cellfun(@str2double, recnamesplit);
        lensel=find(celllen>=3 & celllen<=4 & ~isnan(celldigits));
        aniid=str2num(recnamesplit{lensel});
    end
    aniid_all=[aniid_all, aniid];   
end
end


