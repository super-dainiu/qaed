function str = char(o)
% CHAR Create character array (string).
% (Octonion overloading of standard Matlab function.)

% Note: the Matlab char function converts arrays of numeric values into
% character strings. This is not what this function does, but the Matlab
% guidance on user-defined classes suggests writing a char function and
% a disp/display function. This advice has been followed.

% TODO Merge this code into DISP. The advice mentioned above is now
% obsolete, in particular the overloading of DISPLAY is not now
% recommended. There is no need for a CHAR function.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1) 

% There are three cases to be handled. The argument is one of: empty, a
% pure octonion, a full octonion.

if isempty(o)
    str = '[] octonion'; % This case must be handled first because an
    return;                % empty octonion is not scalar, and the
end                        % next check would fail.

if ~isscalar(o)
    error('char cannot handle a vector or a matrix.')
end

f = '%.4g'; % Control over the format of each numeric value.

% First deal with the scalar part (if it exists) and the X or I component
% because the latter needs to be handled differently according as to
% whether the scalar parts exists or is empty.

t = x(o.a);
if isempty(s(o.a))
    % No scalar part exists, so the X component is the first component to
    % be output.

    str = [nullminus(t) comp2str(t) ' * I'];
else
    % There is a scalar part, so we start with that, and then add in the X
    % component.

    str = [nullminus(s(o.a)) comp2str(s(o.a)) ' ' ...
           plusminus(t)  ' ' comp2str(t) ' * I'];
end

% Now output the Y/Z/A/B/C/D components, with J/K/L/M/N and suitable
% signs and spaces.

t = y(o.a); str = [str ' ' plusminus(t) ' ' comp2str(t) ' * J'];
t = z(o.a); str = [str ' ' plusminus(t) ' ' comp2str(t) ' * K'];
t = s(o.b); str = [str ' ' plusminus(t) ' ' comp2str(t) ' * L'];
t = x(o.b); str = [str ' ' plusminus(t) ' ' comp2str(t) ' * M'];
t = y(o.b); str = [str ' ' plusminus(t) ' ' comp2str(t) ' * N'];
t = z(o.b); str = [str ' ' plusminus(t) ' ' comp2str(t) ' * O'];

    % TODO See the quaternion/char function, addnewline function. The same
    % idea needs to be implemented here, to avoid garbled command window
    % output with long symbolic expressions.

    function S = comp2str(X)
        % Create a string representation of one component of the octonion,
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
            if ~isempty(regexp(S, '[\^/+-]', 'once'))
                S = ['(' S ')'];
            end
        else
            error(['Cannot convert octonion with ', class(z(o.a)), ...
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

% $Id: char.m 1151 2022-09-28 16:30:10Z sangwine $
