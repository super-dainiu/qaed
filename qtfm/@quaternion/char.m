function str = char(q)
% CHAR Create character array (string).
% (Quaternion overloading of standard Matlab function.)

% Note: the Matlab char function converts arrays of numeric values into
% character strings. This is not what this function does, but the Matlab
% guidance on user-defined classes suggests writing a char function and
% a disp/display function. This advice has been followed.

% TODO Merge this code into DISP. The advice mentioned above is now
% obsolete, in particular the overloading of DISPLAY is not now
% recommended. There is no need for a CHAR function. CAUTION test_qfft.m
% calls char to format a string for error output, so this function is
% needed, and removing it would require the test code to be modifiied.

% Copyright © 2005, 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

% There are three cases to be handled. The argument is one of: empty, a
% pure quaternion, a full quaternion.

if isempty(q)
    str = '[] quaternion'; % This case must be handled first because an
    return;                % empty quaternion is not scalar, and the
end                        % next check would fail.

if ~isscalar(q)
    error('char cannot handle a vector or a matrix of quaternions.')
end

% TODO The format here should be dependent on the current setting of the
% format (using the format command). This can be obtained using f = get(0,
% 'format') but unfortunately it returns a string like 'Short' rather than
% a C-style format code. There doesn't seem to be a built-in way to map one
% to the other, so we would have to do it here, with a switch statement.

f = '%.4g'; % Control over the format of each numeric value.

% Lines of greater than about 5950 characters seem to cause text overwrite
% in the MATLAB command window. To prevent this from happening, we insert a
% newline after component of the quaternion if the length exceeds 5000.
% This doesn't fix a problem if one component has a length greater than
% 5000, but that is done in the MATLAB char function. (Long strings are
% only likely in the symbolic case.)

last_newline = 0; % Used below to control output of newlines as
                  % needed (only when the line lengths exceed 5000 chars).

if isempty(q.w)
  % The scalar part is empty, so we begin with the x component of q.

  str = addnewline([nullminus(q.x) comp2str(q.x) ' * I']);
else
  % There is a scalar part, so we start with that, and then add in the x
  % component of q.
  % TODO What if the scalar part is zero, should we elide it?
  str = addnewline([     nullminus(q.w)     comp2str(q.w) ' ']);
  str = addnewline([str, plusminus(q.x) ' ' comp2str(q.x) ' * I']);
end

str = addnewline([str ' ' plusminus(q.y) ' ' comp2str(q.y) ' * J']);
str = addnewline([str ' ' plusminus(q.z) ' ' comp2str(q.z) ' * K']);

    function S = addnewline(X)
        % Add a newline if the length added since the last newline is
        % greater than 5000.
        if length(X) - last_newline > 5000
            S = [X, newline];
            last_newline = length(S);
        else
            S = X;
        end
    end

    function S = comp2str(X)
        % Create a string representation of one component of the quaternion,
        % which may be numeric (real or complex) or a symbolic, or logical.
 
        if isnumeric(X)
            if isreal(X)
                S = num2str(abs(X), f);
            else
                S = ['(' num2str(X, f) ')']; % Complex numeric.
            end
        elseif islogical(X)
            % Logical values cannot be complex (try it, any way you try to
            % make a complex with logical elements, they get converted to
            % double). Hence we can display logicals using the same code as
            % for numerics. We do it separately in case in the future we
            % decide to display logical values as T/F or true/false.
            S = num2str(abs(X), f);
        elseif isa(X, 'sym')
            % We form a string here with surrounding parentheses in some
            % cases and without in others. If the string contains an
            % operator like + ^ or / we put the parentheses in place. If
            % the string is just a variable name, we don't need to,
            % although we don't test for that explicitly at present. This
            % may need review.
            S = char(X);
            % We have just called MATLAB's char function on a single
            % component of a quaternion. If the output is longer than 5000
            % characters or so, it is going to be garbled on the command
            % window, so we truncate it and add an ellipsis and newline.
            % Not ideal, but it's a workaround for now (September 2022).
            if length(S) > 5000
                S = [S(1:5000), '...', newline];
            end
            if ~isempty(regexp(S, '[\^/+-]', 'once'))
                S = ['(' S ')'];
            end
        else
            error(['Cannot convert quaternion with ', class(q.x), ...
                   'elements to a character representation'])
        end
    end
end

function S = plusminus(X)
% Extracts the sign of X and returns the character '-', '+'. The sign
% function doesn't exist for non-numeric values (e.g. logical) which is why
% we check whether X is numeric.

if isnumeric(X) && sign(X) == -1
    S = '-';
else
    S = '+';
end
end

function S = nullminus(X)
% Returns a space or minus according to the sign of X.
if isnumeric(X) && sign(X) == -1
    S = '-';
else
    S = ' ';
end

end

% $Id: char.m 1154 2022-10-13 09:27:42Z sangwine $
