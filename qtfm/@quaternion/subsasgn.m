function r = subsasgn(a, ss, b)
% SUBSASGN Subscripted assignment.
% (Quaternion overloading of standard Matlab function.)

% Copyright © 2005, 2010, 2015, 2020
%               Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

switch ss.type
case '()'
    if length(ss) ~= 1
        error(['Multiple levels of subscripting are not supported ', ...
               'for quaternions.'])
    end
    
    if ~isa(a, 'quaternion')
        % The left hand side of the assignment is NOT a quaternion. Matlab
        % is able to handle an assignment such as r(2) = randn, even if r
        % does not exist (an array is created before the assignment is
        % performed). We choose not to handle this here, because an
        % assignment to an existing numeric array such as r(2) = randq
        % does not work because Matlab attempts to cast the quaternion to a
        % numeric, and handling that seems a lot of trouble for what it's
        % worth.

        error(['Left side must be a quaternion in subscripted ', ...
               'assignment if right side is a quaternion.'])    
    end

    if ~isa(b, 'quaternion')

           % Argument a (the left-hand side of the assignment) is a
           % quaternion, but b (the right-hand side) is not (it could be
           % double, for example, zero; or it could be empty). To handle
           % some of these cases, we need to convert b to a quaternion.
           
           if isempty(b)
               % The right parameter is empty. This occurs when the
               % subscripted assignment is attempting to delete some part
               % of the array a, for example a(2, :) = [] to delete the
               % second row. There is a subtle issue here. We have to
               % explicitly assign [], rather than attempting to assign the
               % empty that is passed to us in b, otherwise an error
               % occurs. The class of the empty array does not matter,
               % since it is not present in the result.

               if ~isempty(a.w)  
                   sa = a.w; sa(ss.subs{:}) = []; a.w = sa;
               end
               
               xa = a.x; xa(ss.subs{:}) = []; a.x = xa;
               ya = a.y; ya(ss.subs{:}) = []; a.y = ya;
               za = a.z; za(ss.subs{:}) = []; a.z = za;

               r = a;
               return
           end
           
           if isnumeric(b) % but not empty because we dealt with that above
               b = quaternion(b); % Convert b to a quaternion (with zero
                                  % vector part), and carry on, the code
                                  % below will handle it.
           else
               error(['Cannot handle right-hand argument of class ', ...
                       class(b)]);
           end
           
    end
    
    % Both parameters, a and b, were, or are now, quaternions. To implement
    % indexed assignment, we operate separately on the components, replace
    % the components with the modified components, then copy to the output
    % parameter. However, we need to deal with the complication that the
    % two sides of the assignment may differ in whether they are pure or
    % full quaternions.

    pa = isempty(a.w);
    pb = isempty(b.w);
        
    if pa && ~pb
        
       % The left side is pure but the right side is full. To handle this
       % we need to convert the left side to a full quaternion.
       
        sa = zeros(size(a), class(a.x)); sa(ss.subs{:}) = b.w; a.w = sa;
        xa = a.x; xa(ss.subs{:}) = b.x;                        a.x = xa;
        ya = a.y; ya(ss.subs{:}) = b.y;                        a.y = ya;
        za = a.z; za(ss.subs{:}) = b.z;                        a.z = za;
        r = a;
        return
    end
    
    if ~pa && pb

        % The left side is full, the right side is pure.  We can cope with
        % this by providing explicit zeros for the scalar part of the right
        % side.
        
        sa = a.w; sa(ss.subs{:}) = zeros(size(b), class(b.x)); a.w = sa;
        xa = a.x; xa(ss.subs{:}) = b.x;                        a.x = xa;
        ya = a.y; ya(ss.subs{:}) = b.y;                        a.y = ya;
        za = a.z; za(ss.subs{:}) = b.z;                        a.z = za;
        r = a;
        return
    end
    
    % a and b are both pure, or both full. These cases can be dealt with
    % together provided we treat the scalar part separately,

    if ~pa % && ~pb, since they must be either both pure, or both full.
        sa = a.w; sa(ss.subs{:}) = b.w; a.w = sa;
    end

    xa = a.x; xa(ss.subs{:}) = b.x; a.x = xa;
    ya = a.y; ya(ss.subs{:}) = b.y; a.y = ya;
    za = a.z; za(ss.subs{:}) = b.z; a.z = za;
    r = a;
    return
case '{}'
    error('Cell array indexing is not valid for quaternions.')
case '.'
    error('Structure indexing is not implemented for quaternions.')
    %
    % Possible use of structure indexing is to implement the following
    % sorts of assignment:
    %
    % q.x = blah
    %
    % However, there are some issues to be considered before implementing
    % such a scheme. Would it work, for example, if q doesn't exist? What
    % should happen if q is pure and q.s = blah occurs? Should q become a
    % full quaternion, or should an error be raised? What about mixed use
    % of array indexing and structure indexing, e.g. q.x(:,1)? Would it
    % work for q.x = q.x .* 2 where the structure indexing occurs on the
    % right hand side as well as on the left. Guess: q.x on the right would
    % be handled by subsref, not subsassgn.
    %
    % (Notes added after discussion with Sebastian Miron, 12 July 2005.)
    %
otherwise
    error('Quaternion subsasgn received an invalid subscripting type.')
end

end

% $Id: subsasgn.m 1113 2021-02-01 18:41:09Z sangwine $
