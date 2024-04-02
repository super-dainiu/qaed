function V = qtfm_version
% Find QTFM version number (of current installation or Sourceforge latest).
% Returns the version as a string, e.g. '2.3'.
%
% Copyright © 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(0, 0); nargoutchk(0, 2);

T = ver('qtfm'); % Get the version of the currently running version. This
V = T.Version;   % comes from the file Contents.m.

% Read the latest version number from Sourceforge.

try
    data = webread('http://sourceforge.net/projects/qtfm/best_release.json');
    F = data.release.filename; % Should have format: '/qtfm/2.3/qtfm_2_3.zip'
    if strcmp(F(2:5), 'qtfm')
        SV = F(7:9); % Sourceforge version
        if ~strcmp(V, SV)
            if str2num(SV) > str2num(V)
                disp(['***** A more recent version (' SV ...
                      ') of QFTM is available at Sourceforge *****'])
            else
                disp(['***** The Sourceforge version (' SV ...
                    ') of QFTM is different to the installed version (' ...
                                 V ') *****'])
            end
        end
    else
        warning('Data retrieved from Sourceforge has unexpected format')
    end
catch
    warning('Unable to access Sourceforge to check the latest version')
end

end

% $Id: qtfm_version.m 1168 2022-11-20 19:07:25Z sangwine $
