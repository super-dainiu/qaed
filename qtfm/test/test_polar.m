function test_polar
% Test code for the quaternion/octonion polar function.

% Copyright © 2018 Stephen J. Sangwine and Nicolas Le Bihan
% See the file : Copyright.m for further details.

disp('Testing polar decomposition ...');

T = 1e-10;

% Test biquaternions.

b = complex(randq(10), randq(10));

% The function called below is not the MATLAB polar function, therefore we
% suppress the Code Compatibility warning about using the polarplot
% function instead.

[t, h, n] = polar(b, 'R'); %#ok<*POLAR>

compare(b, t .* h .* n, T, 'quaternion/polar failed test 1.');

[t, h, n] = polar(b, 'L');

compare(b, h .* t .* n, T, 'quaternion/polar failed test 2.');

b = unit(b);

[t, h] = polar(b, 'R');

compare(b, t .* h, T, 'quaternion/polar failed test 3.');

[t, h] = polar(b, 'L');

compare(b, h .* t, T, 'quaternion/polar failed test 4.');

% Now repeat for octonions.

T = 1e-10;

b = complex(rando(10), rando(10));

[t, h, n] = polar(b, 'R');

compare(b, t .* h .* n, T, 'octonion/polar failed test 1.');

[t, h, n] = polar(b, 'L');

compare(b, h .* t .* n, T, 'octonion/polar failed test 2.');

b = unit(b);

[t, h] = polar(b, 'R');

compare(b, t .* h, T, 'octonion/polar failed test 3.');

[t, h] = polar(b, 'L');

compare(b, h .* t, T, 'octonion/polar failed test 4.');


disp('Passed');

end

% $Id: test_polar.m 1113 2021-02-01 18:41:09Z sangwine $
