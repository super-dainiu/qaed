% Quaternion Toolbox for Matlab (QTFM)
%
% Quaternion class definition and constructor method/function.

% Copyright © 2005, 2009, 2010, 2011, 2016, 2020, 2021
%               Stephen J. Sangwine and Nicolas Le Bihan.
%
% See the file : Copyright.m for further details.

classdef (InferiorClasses = {?sym}) quaternion % Note 1.
    properties (SetAccess = 'private', Hidden = false)
        w = [];
        x = [];
        y = [];
        z = [];
    end
    methods
        function q = quaternion(a0, a1, a2, a3)
            % QUATERNION   Construct quaternions from components. Accepts
            % the following possible arguments, which may be scalars,
            % vectors or matrices; numeric, logical or symbolic. In
            % addition character arrays may be passed to create symbolic
            % components with the given names (these must not already
            % exist, otherwise they should be passed by name).
            %
            % No argument     - returns an empty quaternion scalar, vector
            %                   or matrix.
            % One argument    - A quaternion argument returns the argument
            %                   unmodified. A non-quaternion argument
            %                   returns the argument in the scalar part and
            %                   supplies a zero vector part with elements
            %                   of the same type as the scalar part.
            % Two arguments   - returns a quaternion, provided the first
            %                   argument is numeric/symbolic/char and the
            %                   second is a compatible pure quaternion.
            % Three arguments - returns a pure quaternion scalar, vector or
            %                   matrix, with an empty scalar part.
            % Four arguments  - returns a full quaternion scalar, vector or
            %                   matrix.
            
            nargoutchk(0, 1) % The number of input arguments is checked
                             % as part of the switch statement below.
            switch nargin
                case 0
                    
                    % Construct an empty quaternion.
                    
                    q.w = []; q.x = []; q.y = []; q.z =[];
                    
                case 1
                    
                    if isa(a0, 'quaternion')
                        q = a0; % a0 is already a quaternion, so return it.
                    else
                        % a0 is not a quaternion ...
                        if isnumeric(a0) || isa(a0, 'sym')
                            % ... but it is numeric or symbolic. We can
                            % handle this, but we have to supply the vector
                            % components of the quaternion, and they must
                            % have the same type as a0 (e.g. double, int8,
                            % sym), for which we use the class option to
                            % the zeros function.
                            
                            t = zeros(size(a0), 'like', a0);
                            q.w = a0; q.x = t; q.y = t; q.z = t;
                        elseif islogical(a0)
                            % ... but it is logical, which we can handle,
                            % with different code, because zeros won't work
                            % for a logical type (sadly).
                            f = false(size(a0));
                            q.w = a0; q.x = f; q.y = f; q.z = f;
                        elseif isa(a0, 'char')
                            % ... but it is a character array,
                            % and a potential name for a symbolic scalar
                            % part. This character array must represent ONE
                            % name (therefore it must be a row vector). We
                            % cannot handle creation of a quaternion with
                            % size greater than [1,1] by this mechanism.
                            
                            if ~isrow(a0)
                                error(['char parameter must be a ', ...
                                       'row vector representing one ', ...
                                       'symbolic variable name.'])
                            end
                            
                            % What follows is tricky code, because we want
                            % to test for or create a symbolic variable in
                            % the caller's workspace, not in this
                            % function's workspace. However, where a0
                            % occurs below, it is evaluated here to yield
                            % the string of characters that was passed in
                            % the first parameter.
                           
                            if evalin('caller', ...
                                     ['exist(''', a0, ''',', '''var'')'])
                                error([a0, ' already exists in the workspace, ', ...
                                           'try using the variable name directly.'])
                                % The above means of course the workspace from
                                % which this function was called.
                            else
                                evalin('caller', ['syms ', a0]); % Create variable..
                                % ... and now reference it and assign to the scalar
                                % part of q.
                                q.w = evalin('caller', ['sym(''', a0, ''')']);
                                q.x = zeros(1, 'like', q.w);
                                q.y = q.x;
                                q.z = q.x;
                            end
                        elseif isstring(a0)
                            % TODO - this case unlike the char case above
                            % could be an ARRAY of strings, requiring
                            % iteration to create the variables. Maybe too
                            % complex to make it worthwhile.
                            error('String parameters are not yet supported.');
                        else
                            error(['Cannot construct a quaternion with a ', ...
                                    class(a0), ' in the scalar part.']);
                        end
                    end
                    
                case 2
                    
                    % In this case, the first argument must be the scalar
                    % part, and thus numeric (or logical), or symbolic, and
                    % the second must be the vector part, and thus a pure
                    % quaternion (with components of the same class and
                    % size as the first argument).
                    
                    % TODO If the first parameter is a symbolic zero, do we
                    % make a pure quaternion and omit the scalar part? What
                    % if the scalar part is symbolic and simplifies to
                    % zero? See the simplify function, which deals with
                    % this.
                                                            
                    if (isnumeric(a0) || islogical(a0)) || isa(a0, 'sym') ...
                            && isa(a1, 'quaternion')

                        if any(size(a0) ~= size(a1))
                            error('Arguments must have the same dimensions.')
                        end
                    
                        c0 = class(a0); c1 = class(a1.x);
                        if isempty(a1.w)
                            if strcmp(c0, c1)
                                % a1 is already a quaternion, so copy it to
                                % the output and insert a0 as the w
                                % component.
                                q = a1; q.w = a0;
                            else
                                error(['Classes of given scalar and vector ', ...
                                       'parts differ: ', c0, ' ', c1])
                            end
                        else
                            error('The second argument must be a pure quaternion.');
                        end
                        
                    elseif isa(a0, 'char') && isa(a1, 'quaternion') && isempty(a1.w)
                        
                        % This needs special treatment because the first
                        % parameter need not have the same size as the
                        % second since it will be an identifier. See case 1
                        % above for commentary on the other aspects here.
                        
                        if ~isa(a1.x, 'sym')
                            error(['Second parameter must have ',  ...
                                   'symbolic components. Found: ', ...
                                   class(a1.x)])
                        end
                        
                        if ~isrow(a0)
                            error(['char parameter must be a ', ...
                                   'row vector representing one ', ...
                                   'symbolic variable name.'])
                        end
                        
                        if any(size(a1) ~= [1,1])
                            error('Second parameter must be of size [1,1].')
                        end
                        
                        q = a1; % a1 is already a (pure) quaternion so use
                                % it for the output value, and then below
                                % construct the scalar part.
                        
                        if evalin('caller', ['exist(''', a0, ''',', '''var'')'])
                            error([a0, ' already exists in the workspace, ', ...
                                       'try using the variable name directly.'])
                        else
                            evalin('caller', ['syms ', a0]);
                            q.w = evalin('caller', ['sym(''', a0, ''')']);
                        end

                    elseif isstring(a0)
                        % TODO - this case unlike the char case above
                        % could be an ARRAY of strings, requiring
                        % iteration to create the variables.
                        error('String parameters are not yet supported.');
                    else
                        error(['First argument must be numeric, ',    ...
                               'logical, symbolic or char; and the ', ...
                               'second must be a pure quaternion.']);
                    end
                    
                case 3
                    
                    % Construct a pure quaternion, since we are given three
                    % arguments and this is the only possibility. All three
                    % arguments must be numeric, logical or symbolic and of
                    % the same class, or be char arrays which will form the
                    % names of symbolic variables.
                    
                    if isnumeric(a0) || islogical(a0) || isa(a0, 'sym')
                        
                        s0 = size(a0);
                        if any(s0 ~= size(a1)) || any(s0 ~= size(a2))
                            error('All three arguments must have the same dimensions.')
                        end
                        
                        c0 = class(a0); c1 = class(a1); c2 = class(a2);
                        
                        if ~(strcmp(c0, c1) && strcmp(c0, c2))
                            error(['All three arguments must be of the', ...
                                   ' same class. Given: ', ...
                                   c0, ' ', c1, ' ', c2]);
                        end
                        
                        q.w = cast([], 'like', a0);
                        q.x = a0; q.y = a1; q.z = a2;
                        
                    elseif isa(a0, 'char')
                        
                        if ~isa(a1, 'char')
                            error(['First parameter is char, but second', ...
                                   ' is not, found ', class(a1)])
                        end
                        if ~isa(a2, 'char')
                            error(['First parameter is char, but third', ...
                                   ' is not, found ', class(a2)])
                        end
                        
                        par = {a0, a1, a2};
                        field = cell(1,3); % Preallocate for sym variables.
                        
                        for j = 1:3
                            a = par{j};
                            if ~isrow(a)
                                error(['char parameters must each be a ', ...
                                       'row vector representing one ', ...
                                       'symbolic variable name.' ...
                                       ' Problem with parameter ', num2str(j)])
                            end
                            if evalin('caller', ['exist(''', a, ''',', '''var'')'])
                                error([a, ' already exists in the workspace, ', ...
                                    'try using the variable name directly.'])
                            else
                                evalin('caller', ['syms ', a]);
                                field{j} = evalin('caller', ['sym(''', a, ''')']);
                            end
                        end
                        q.w = sym.empty;
                        [q.x, q.y, q.z] = field{:}; % RHS is equiv: deal(field)
                    else
                        error(['Arguments must be numeric, logical, ', ...
                               'symbolic or char. Cannot handle ', class(a0)]);    
                    end
                    
                case 4 % Return a full quaternion.
                    
                    % TODO What if the scalar part supplied is zero
                    % (numeric or symbolic)? Should we not make a pure
                    % quaternion? See the simplify function, which deals
                    % with this.
                    
                    if isnumeric(a0) || islogical(a0) || isa(a0, 'sym')
                        
                        s0 = size(a0);
                        if any(s0 ~= size(a1)) || ...
                           any(s0 ~= size(a2)) || ...
                           any(s0 ~= size(a3))
                            error('All four arguments must have the same dimensions.')
                        end
                        
                        c0 = class(a0); c1 = class(a1); ...
                        c2 = class(a2); c3 = class(a3);
                        
                        if ~(strcmp(c0, c1) && strcmp(c0, c2) && strcmp(c0, c3))
                            error(['All four arguments must be of the', ...
                                   ' same class. Given: ', ...
                                   c0, ' ', c1, ' ', c2, ' ', c3]);
                        end
                        
                        q.w = a0; q.x = a1; q.y = a2; q.z = a3;
                        
                    elseif isa(a0, 'char')
                        
                        if ~isa(a1, 'char')
                            error(['First parameter is char, but second', ...
                                   ' is not, found ', class(a1)])
                        end
                        if ~isa(a2, 'char')
                            error(['First parameter is char, but third', ...
                                   ' is not, found ', class(a2)])
                        end
                        if ~isa(a3, 'char')
                            error(['First parameter is char, but fourth', ...
                                   ' is not, found ', class(a3)])
                        end                        

                        % This needs special treatment because the four
                        % parameters could well have different sizes - they
                        % are identifiers. See case 1 above for commentary
                        % on the other aspects here.
                        
                        par = {a0, a1, a2, a3};
                        field = cell(1,4); % Preallocate.
                        
                        for j = 1:4
                            a = par{j};
                            if ~isrow(a)
                                error(['char parameters must each be a ', ...
                                       'row vector representing one ', ...
                                       'symbolic variable name.' ...
                                       ' Problem with parameter ', num2str(j)])
                            end

                            if evalin('caller', ['exist(''', a, ''',', '''var'')'])
                                error([a, ' already exists in the workspace, ', ...
                                    'try using the variable name directly.'])
                            else
                                evalin('caller', ['syms ', a]);
                                field{j} = evalin('caller', ['sym(''', a, ''')']);
                            end
                        end
                        [q.w, q.x, q.y, q.z] = field{:}; % RHS is equiv: deal(field)
                    else
                        error(['Arguments must be numeric, logical, ', ...
                               'symbolic or char. Cannot handle ', class(a0)]);    
                    end
                    
                    % If more than 4 arguments are passed, Matlab catches
                    % it and outputs the error message 'Too many input
                    % arguments'. Hence we do not need an 'otherwise' case,
                    % since it could not be reached.
            end
        end
        
        function n = numArgumentsFromSubscript(~,~,~)
            % Introduction of this function with Matlab R2015b and QTFM 2.4
            % permitted numel to revert to its obvious function of
            % providing the number of elements in an array.
            n = 1;
        end

        % TODO Consider how DUMP and DETAILS interact. Should we replace
        % DUMP with DETAILS or make DETAILS a synonym? DETAILS was
        % introduced in R2013b, so was not in existence when the toolbox
        % was first developed.

        function dump(q)
            % Output a diagnostic dump showing the content of the quaternion
            % q for debugging purposes (particularly in development when the
            % quaternion somehow contains inconsistent data in each field).
            s = inputname(1);
            if ~isempty(s)
                s = [s, '.']; % If q has a name, display 'name'.w etc.
            end
            disp([s, 'w: [', num2str(size(q.w)), '] ', class(q.w)])
            disp([s, 'x: [', num2str(size(q.x)), '] ', class(q.x)])
            disp([s, 'y: [', num2str(size(q.y)), '] ', class(q.y)])
            disp([s, 'z: [', num2str(size(q.z)), '] ', class(q.z)])
        end
    end
    methods (Static = true)
        function q = empty(varargin)
            % This function makes it possible to use the dotted notation
            % quaternion.empty or quaternion.empty(0,1) to create an empty
            % array.
            % TODO Make it possible to write quaternion.empty('double') and
            % specify the class of the empty components. The code for this
            % is present in the Clifford Multivector Toolbox for Matlab and
            % would require very little adaptation to work here, as it is
            % mostly concerned with parameter validation.
            d = double.empty(varargin{:});
            q = quaternion(d, d, d, d);
        end
    end
end

% Note 1. The InferiorClasses attribute ensures that in an expression
% involving a symbolic expression or variable, and a quaternion, the
% quaternion function will be called, and not the symbolic function. This
% is essential because the symbolic toolbox knows nothing about
% quaternions. Example: the expression a + quaternion(x, y, z) where all
% the variables are symbolic, must result in a call to @quaternion/plus. So
% must the expression quaternion(x, y, z) + a. Without the attribute, only
% the second expression will be passed to @quaternion/plus, the first one
% would result in a call of the symbolic plus function and therefore an
% error.

% $Id: quaternion.m 1144 2022-09-07 14:43:23Z sangwine $
