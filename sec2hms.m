function hmsStr = sec2hms( tSec )
%%sec2hms  Convert time in seconds to HH:MM:SS format.
%  [George-Gate @2017-06-08]
%  [Usage]
%       hmsStr=sec2hms( timeInSeconds )
%   
    if tSec<0
        error('Time must be a positive double number.');
    end
    if tSec<60
        hmsStr=[num2str(tSec,'%.3f'),'s'];
    elseif tSec<3600
        tSec=round(tSec);
        MM=floor(tSec/60);
        SS=tSec-MM*60;
        hmsStr=[num2str(MM,'%02d'),':',num2str(SS,'%02d')];
    else
        tSec=round(tSec);
        HH=floor(tSec/3600);
        tSec=tSec-HH*3600;
        MM=floor(tSec/60);
        SS=tSec-MM*60;
        hmsStr=[num2str(HH,'%02d'),':',num2str(MM,'%02d'),':',num2str(SS,'%02d')];
    end
end

