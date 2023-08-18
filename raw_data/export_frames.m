files = dir ('*.etk');
k = length(files);

for c = 1:k
    file = files(c).name;
    load('-mat',file);
    filename = split(file,'.');
    directory = strcat(filename{1},"_frames");
    mkdir(directory);
    for i = 1:length(Tracked.Frames)
    C = Tracked.Frames(i).Cells;
    if length(C) == 1
        fulltab = struct2table(C,AsArray=true);
    else    
        fulltab = struct2table(C);
    end    
    tab = fulltab(:,["pos","progenitor","descendants"]);
    writetable(tab,strcat(directory,"/frame",int2str(i),".txt"));
    end    
end

