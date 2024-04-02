function r = subsasgn(a, ss, b)
% SUBSASGN Subscripted assignment.
% (Octonion overloading of standard Matlab function.)

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

switch ss.type
case '()'
    if length(ss) ~= 1
        error(['Multiple levels of subscripting are not supported ', ...
               'for octonions.'])
    end
    
    if ~isa(a, 'octonion')
        % See the quaternion version of this function for a discussion as
        % to why we don't attempt to handle this case.
        error(['Left side must be an octonion in subscripted ', ...
               'assignment if right side is an octonion.'])    
    end

    if ~isa(b, 'octonion')

        % Argument a (the left-hand side of the assignment) is an
        % octonion, but b (the right-hand side) is not (it could be
        % double, for example, zero; or it could be empty). To handle
        % some of these cases, we need to convert b to an octonion.
           
        if isempty(b)
            % The right hand side of the assignment is empty. See the
            % quaternion version of this function for a discussion as to
            % why we must pass an explicit empty, and not b.a or b.b here.
            
            r = a; % Copy the input array to create an octonion result.
            
            r.a = subsasgn(a.a, ss, []); % These calls are to the
            r.b = subsasgn(a.b, ss, []); % quaternion subsagn function.
            return
        end
        
        if isnumeric(b) % but not empty because we dealt with that above
            b = octonion(b); % Convert b to an octonion (with zero vector
                             % part), and carry on, the code below will
                             % handle it.
        else
            error(['Cannot handle right-hand argument of class ', ...
                    class(b)]);
        end

    end
    
    % Both parameters, a and b, were, or are now, octonions. To implement
    % indexed assignment, we operate separately on the quaternion
    % components, replace the components with the modified components, then
    % copy to the output parameter. The a component of either octonion
    % could be pure or full according as to whether the octonion itself is
    % pure or full. All the complications are handled or detected by the
    % quaternion subsasgn function, so we don't need to deal with it here.
    
    % To perform subscripted assignment, we split the octonion into two
    % quaternions and operate on the two quaternions separately. Finally we
    % re-assemble the two quaternion results to make the octonion result.
    
    r = a; % Copy the input array to create an octonion result.
    r.a = subsasgn(a.a, ss, b.a); % These calls are to the quaternion
    r.b = subsasgn(a.b, ss, b.b); % subsagn function, of course.
case '{}'
    error('Cell array indexing is not valid for octonions.')
case '.'
    error('Structure indexing is not implemented for octonions.')
otherwise
    error('Octonion subsasgn received an invalid subscripting type.')
end

end

% $Id: subsasgn.m 1113 2021-02-01 18:41:09Z sangwine $
