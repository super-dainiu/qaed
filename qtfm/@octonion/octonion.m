% Quaternion Toolbox for matlab (QTFM)
%
% Octonion class definition and constructor method/function. The octonion
% class is built on top of the quaternion class and unlike the quaternion
% class, needs a class definition to establish that it is superior to the
% quaternion class. This means that an expression such as q + o which
% attempts to add a quaternion to an octonion (not possible of course) will
% always call the octonion function plus.m. Therefore any error handling
% may be included only in the octonion code and we do not have to rewrite
% existing quaternion code. This also fits with the primary purpose of QTFM
% which is to be a quaternion toolbox - the octonion class is a bonus, but
% likely to remain somewhat experimental for some years to come.

% Copyright © 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

classdef (InferiorClasses = {?quaternion, ?sym}) octonion % Note 1.
    properties (SetAccess = 'private', Hidden = false)
        a = quaternion();
        b = quaternion();
    end
    methods
        function o = octonion(a0, a1, a2, a3, a4, a5, a6, a7)
            % OCTONION   Construct octonions from components. Accepts the
            % following possible arguments, which may be scalars, vectors
            % or matrices; numeric, logical or symbolic. In addition
            % character arrays may be passed to create symbolic components
            % with the given names (these must not already exist, otherwise
            % they should be passed by name).
            %
            % No argument     - returns an empty octonion scalar, vector or
            %                   matrix.
            % One argument    - An octonion argument returns the argument
            %                   unmodified. A non-octonion argument returns
            %                   the argument in the scalar part and
            %                   supplies a zero vector part with elements
            %                   of the same type as the scalar part.
            % Two arguments   - returns an octonion, two cases are
            %                   supported:
            %                   1. the first argument is numeric and the
            %                   second is a pure octonion, with the same
            %                   size, and the elements of the pure octonion
            %                   must have the same class as the scalar;
            %                   2. both arguments are quaternions of the
            %                   same size, with elements of the same class.
            % Seven arguments - returns a pure octonion scalar, vector or
            %                   matrix, with an empty scalar part.
            % Eight arguments - returns a full octonion scalar, vector or
            %                   matrix.
            
            nargoutchk(0, 1) % The number of input arguments is checked
                             % as part of the switch statement below.
            switch nargin
                case 0
                    
                    % Construct an empty octonion with two empty quaternion
                    % components.
                    
                    o.a = quaternion(); o.b = o.a;
                    
                case 1
                    
                    if isa(a0, 'octonion')
                        o = a0; % a0 is already an octonion, so return it.
                    else
                        % a0 is not a octonion ...
                        if isnumeric(a0)
                            % ... but it is numeric. We can handle this,
                            % but we have to supply the vector components
                            % of the octonion, and they must have the same
                            % type as a0 (e.g. double, int8).
                            
                            o.a = quaternion(a0);
                            o.b = zerosq(size(a0), class(a0));
                        elseif islogical(a0)
                            % ... but it is logical, which we can handle,
                            % with different code, because zeros won't work
                            % for a logical type (sadly).
                            f = false(size(a0));
                            o.a = quaternion(a0, f, f, f);
                            o.b = quaternion(f,  f, f, f);
                        elseif isa(a0, 'sym')
                            % ... but it is symbolic, and we can handle
                            % that. Actually it could be done with the same
                            % code as for numeric, but we keep it separate
                            % until we understand it better.
                            t = zeros(size(a0), 'like', a0);
                            o.a = quaternion(a0, t, t, t);
                            o.b = quaternion(t,  t, t, t);
                        elseif isa(a0, 'char')
                            % ... but it is a character array, and a
                            % potential name for a symbolic scalar part.
                            % What follows is tricky code, because we want
                            % to test for or create a symbolic variable in
                            % the caller's workspace, not in this
                            % function's workspace. However, where a0
                            % occurs below, it is evaluated here to yield
                            % the string of characters that was passed in
                            % the first parameter.
                            
                            if evalin('caller', ['exist(''', a0, ''',', '''var'')'])
                                error([a0, ' already exists in the workspace, ', ...
                                    'try using the variable name directly.'])
                                % The above means of course the workspace
                                % from which this function was called.
                            else
                                evalin('caller', ['syms ', a0]); % Create variable..
                                % ... and now reference it and assign to
                                % the scalar part of o.
                                t = evalin('caller', ['sym(''', a0, ''')']);
                                z = zeros(size(t), 'like', t);
                                o.a = quaternion(t, z, z, z);
                                o.b = quaternion(z, z, z, z);
                            end
                        else
                            % ... or it isn't numeric.
                            error(['Cannot construct an octonion with a ',...
                                class(a0), ' in the scalar part.']);
                        end
                    end
                    
                case 2
                    % In this case, the first argument must be the scalar
                    % part, and thus numeric (or logical), or symbolic, and
                    % the second must be the vector part, and thus a pure
                    % octonion (with components of the same class as the
                    % first argument).
                    
                    % TODO If the first parameter is a symbolic zero, do we
                    % make a pure octonion and omit the scalar part? What if
                    % the scalar part is symbolic and simplifies to zero? See
                    % the simplify function, which deals with this. See also
                    % the quaternion constructor for the same issue.
                    
                    if (isnumeric(a0) || islogical(a0) || isa(a0, 'sym')) ...
                            && isa(a1, 'octonion')
                        
                        if any(size(a0) ~= size(a1))
                            error('Arguments must have the same dimensions')
                        end
                        
                        c0 = class(a0); c1 = class(x(a1.a));
                        if isempty(s(a1.a))
                            if strcmp(c0, c1)
                                % a1 is already an octonion, so copy it to
                                % the output and insert a0 as the scalar
                                % component. TODO It would be better here
                                % if we could just assign a0 to the scalar
                                % part of o.a and thus avoid copying three
                                % components onto themselves.
                                o = a1; o.a = quaternion(a0, v(o.a));
                            else
                                error(['Classes of given scalar and vector parts', ...
                                    ' differ: ', c0, ' ', c1])
                            end
                        else
                            error('The second argument must be a pure octonion.');
                        end
                    elseif isa(a0, 'quaternion') && isa(a1, 'quaternion')
                        if ispure(a1)
                            error('Second quaternion argument cannot be pure.')
                        end
                        % But the first quaternion argument may be pure,
                        % and we make a pure octonion (because the second
                        % quaternion argument must be full, due to the
                        % preceding check).
                        c0 = class(x(a0));
                        c1 = class(x(a1));
                        if strcmp(c0, c1)
                            % All is well, make the octonion from the two
                            % quaternion components.
                            o.a = a0;
                            o.b = a1;
                        else
                            error(['Classes of quaternion elements must', ...
                                ' agree, found: ' c0, ' and ', c1])
                        end
                    elseif isa(a0, 'char') && ...
                           isa(a1, 'octonion') && isempty(s(a1.a))
                        
                        % We are dealing here with a character array
                        % representing a symbolic scalar part, so a1 must
                        % be a symbolic pure octonion.
                        
                        % This needs special treatment because the first
                        % parameter need not have the same size as the
                        % second since it will be an identifier.
                        
                        if ~isa(x(a1.a), 'sym')
                            error(['Second parameter must have ',  ...
                                   'symbolic components. Found: ', ...
                                   class(x(a1.a))])
                        end
                        
                        o = a1; % a1 is already a (pure) octonion so use it
                        % for the output value, and then below insert the
                        % scalar part.
                        % TODO A better idea would be to assign the scalar
                        % part to the w component of o.a, but we don't have
                        % a way to do this (yet). Hence the clunky method
                        % actually used.
                        
                        if evalin('caller', ['exist(''', a0, ''',', '''var'')'])
                            error([a0, ' already exists in the workspace, ', ...
                                'try using the variable name directly.'])
                        else
                            evalin('caller', ['syms ', a0]);
                            temp = evalin('caller', ['sym(''', a0, ''')']);
                            o.a = quaternion(temp, vector(o.a)); % TODO Clunky!
                        end
                    else
                        error(['Parameters must be numeric/logical/sym + '...
                            'pure octonion, or two quaternions.']);
                    end
                    
                case 7
                    
                    % Construct a pure octonion, since we are given seven
                    % arguments and this is the only possibility. All seven
                    % arguments must be numeric, logical, or symbolic and
                    % of the same class.
                    
                    if isnumeric(a0) || islogical(a0) || isa(a0, 'sym')
                        
                        s0 = size(a0);
                        if any(s0 ~= size(a1)) || any(s0 ~= size(a2)) || ...
                           any(s0 ~= size(a3)) || any(s0 ~= size(a4)) || ...
                           any(s0 ~= size(a5)) || any(s0 ~= size(a6))
                            error('All seven arguments must have the same dimensions')
                        end
                        
                        c0 = class(a0); c1 = class(a1); c2 = class(a2); ...
                        c3 = class(a3); c4 = class(a4); c5 = class(a5); ...
                        c6 = class(a6);

                        % We test only the first argument for numeric,
                        % logical, or symbolic, because if the classes
                        % match, the others must also be
                        % numeric/logical/symbolic.
                        
                        if ~(strcmp(c0, c1) && strcmp(c0, c2) && ...
                             strcmp(c0, c3) && strcmp(c0, c4) && ...
                             strcmp(c0, c5) && strcmp(c0, c6))
                            error(['All seven arguments must be numeric ', ...
                                   'and of the same class. Given: ', ...
                                   c0, ' ', c1, ' ', c2, ' ', ...
                                   c3, ' ', c4, ' ', c5, ' ', c6]);
                        end
                        o.a = quaternion(    a0, a1, a2);
                        o.b = quaternion(a3, a4, a5, a6);

                    elseif isa(a0, 'char')
                        
                        % TODO The other args must also be char so we need
                        % to check this. See TODO in quaternion
                        % constructor.
                        
                        % This needs special treatment because the seven
                        % parameters could well have different sizes - they
                        % are identifiers.
                        
                        par = {a0, a1, a2, a3, a4, a5, a6};
                        var = cell(1,7); % Preallocate.
                        
                        for j = 1:7
                            a = par{j};
                            if evalin('caller', ['exist(''', a, ''',', '''var'')'])
                                error([a, ' already exists in the workspace, ', ...
                                    'try using the variable name directly.'])
                            else
                                evalin('caller', ['syms ', a]);
                                var{j} = evalin('caller', ['sym(''', a, ''')']);
                            end
                        end
                        o.a = quaternion(var{1}, var{2}, var{3});
                        o.b = quaternion(var{4}, var{5}, var{6}, var{7});
                    end
                    
                case 8 % Return a full octonion. All eight arguments must be
                       % numeric, logical, or symbolic and of the same class.
                    
                    if isnumeric(a0) || islogical(a0) || isa(a0, 'sym')

                        s0 = size(a0);
                        if any(s0 ~= size(a1)) || any(s0 ~= size(a2)) || ...
                           any(s0 ~= size(a3)) || any(s0 ~= size(a4)) || ...
                           any(s0 ~= size(a5)) || any(s0 ~= size(a6)) || ...
                           any(s0 ~= size(a7))
                            error('Arguments must have the same dimensions')
                        end
                        
                        c0 = class(a0); c1 = class(a1); c2 = class(a2); c3 = class(a3);
                        c4 = class(a4); c5 = class(a5); c6 = class(a6); c7 = class(a7);
                        
                        % We test only the first argument for numeric, because if
                        % the classes match, the other three must also be numeric.
                        
                        if ~(strcmp(c0, c1) && strcmp(c0, c2) && ...
                             strcmp(c0, c3) && strcmp(c0, c4) && ...
                             strcmp(c0, c5) && strcmp(c0, c6))
                            error(['All eight arguments must be numeric ', ...
                                   'and of the same class. Given: ', ...
                                   c0, ' ', c1, ' ', c2, ' ', c3, ' ', ...
                                   c4, ' ', c5, ' ', c6, ' ', c7]);
                        end
                        o.a = quaternion(a0, a1, a2, a3);
                        o.b = quaternion(a4, a5, a6, a7);
                    elseif isa(a0, 'char')
                        
                        % This needs special treatment because the eight
                        % parameters could well have different sizes - they
                        % are identifiers.
                        
                        par = {a0, a1, a2, a3, a4, a5, a6, a7};
                        var = cell(1,8); % Preallocate.
                        
                        for j = 1:8
                            a = par{j};
                            if evalin('caller', ['exist(''', a, ''',', '''var'')'])
                                error([a, ' already exists in the workspace, ', ...
                                    'try using the variable name directly.'])
                            else
                                evalin('caller', ['syms ', a]);
                                var{j} = evalin('caller', ['sym(''', a, ''')']);
                            end
                        end
                        o.a = quaternion(var{1}, var{2}, var{3}, var{4});
                        o.b = quaternion(var{5}, var{6}, var{7}, var{8});
                    end
                otherwise
                    error('Octonion constructor takes 0, 1, 2, 7 or 8 arguments.');
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

        function dump(o)
            % Output a diagnostic dump showing the content of the octonion
            % o for debugging purposes (particularly in development when
            % the octonion somehow contains inconsistent data in each
            % field).
            s = inputname(1);
            if ~isempty(s)
                s = [s, '.']; % If o has a name, display 'name'.w etc.
            end
            q = o.a;
            disp([s, 'w: [', num2str(size(q.w)), '] ', class(q.w)])
            disp([s, 'x: [', num2str(size(q.x)), '] ', class(q.x)])
            disp([s, 'y: [', num2str(size(q.y)), '] ', class(q.y)])
            disp([s, 'z: [', num2str(size(q.z)), '] ', class(q.z)])
            q = o.b;
            disp([s, 'a: [', num2str(size(q.w)), '] ', class(q.w)])
            disp([s, 'b: [', num2str(size(q.x)), '] ', class(q.x)])
            disp([s, 'c: [', num2str(size(q.y)), '] ', class(q.y)])
            disp([s, 'd: [', num2str(size(q.z)), '] ', class(q.z)])
        end
    end
    methods (Static = true)
        function q = empty(varargin)
            % This function makes it possible to use the dotted notation
            % octonion.empty or octonion.empty(0,1) to create an empty
            % array.
            % TODO Make it possible to write octonion.empty('double') and
            % specify the class of the empty components.
            d = double.empty(varargin{:});
            q = octonion(d, d, d, d, d, d, d, d);
        end
    end
end

% Note 1. The InferiorClasses attribute ensures that in an expression
% involving a symbolic expression or variable, and an octonion, the
% octonion function will be called, and not the symbolic function. This is
% essential because the symbolic toolbox knows nothing about octonions.
% Example: the expression a + octonion(b, c, d, e, f, g, h) where all the
% variables are symbolic, must result in a call to @octonion/plus. So must
% the expression octonion(b, c, d, e, f, g, h) + a. Without the attribute,
% only the second expression will be passed to @octonion/plus, the first
% one would result in a call of the symbolic plus function and therefore an
% error. Similarly, but for different reasons, an expression involving an
% octonion and a quaternion must result in a call to the relevant octonion
% function. This is so that the octonion function can deal with the
% resulting error, leaving quaternion code to deal only with quaternion
% errors.

% $Id: octonion.m 1156 2022-10-13 16:12:30Z sangwine $
