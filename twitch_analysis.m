function [amp,t95]=twitch_analysis(twitch,t)

    amp=max(twitch)-min(twitch);

    [pks,idx]=findpeaks(-abs(twitch-min(twitch)-amp*0.05));
    [~,i]=sort(pks,'descend');
    idx=idx(i);
    t95=abs(t(idx(2))-t(idx(1)));
    


end