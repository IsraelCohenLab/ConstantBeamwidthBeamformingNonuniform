function varargout =  patternCustomEditted(MagE, theta, phi, varargin)
%PATTERNCUSTOM  Plot 2D and 3D radiation pattern of the antenna in 
%               polar and rectangular coordinate system.
%
% patternCustom(MagE, theta, phi) plots the 2D and 3D radiation pattern
% over the specified phi and theta vectors.
%
% patternCustom(__, Name, Value) plots the 2D and 3D radiation pattern from
% the input data with additional options specified by the Name, Value
% pairs.
%
% hPlot = patternCustom(___) returns handles of the lines or surface in the 
% figure window, using any of the input arguments in the previous syntaxes.
%
% Input Arguments
%
% MagE: It can be either the Magnitude of the quantity plotted specified as
%       a vector or matrix, a phased.CustomAntennaElement object or a
%       phased.IsotropicAntennaElement object.
%       If the quantity is a vector, it should be of the same size
%       as phi and theta. If the quantity is a matrix, it should be of size 
%       phi x theta.
%
% theta: Theta angles in spherical coordinate system specified as a vector.
%        This needs to be provided only when MagE is a vector or matrix.
%
% phi:  Phi angles in the spherical coordinate system specified as a vector.
%       This needs to be provided only when MagE is a vector or matrix.
%       
%
% Below are the list of Name-Value pairs available in the patternCustom
% function:
%
% CoordinateSystem: Coordinate system to visualize the radiation pattern
% specified as polar | rectangular. The default is polar.
%
% Slice: Plane to visualize the 2D data specified as phi | theta. There are
% no defaults for this pair.
%
% SliceValue: Values for the Slice (phi or theta) specified as a vector.
% There are no defaults for this pair.

%   Copyright 2019 The MathWorks, Inc.
if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

% Parsing through the inputs
parseobj = inputParser;
parseobj.FunctionName = 'patternCustom';

% Checking whether the MagE input is a vector or a matrix
if ~isobject(MagE)
    if isvector(MagE)
        narginchk(3,9);
        % Making sure the size of theta and phi is same
        if numel(theta) == numel(phi)
            typeValidationMat = @(x) validateattributes(x,{'numeric'},      ...
                {'vector','numel',numel(theta),'nonempty','real'}, 'patternCustom');
        else
            error(message('shared_channel:patternCustom:DimensionMismatch'));
        end
    elseif ismatrix(MagE)
        narginchk(3,9);
        typeValidationMat = @(x) validateattributes(x,{'numeric'},          ...
            {'size',[numel(phi),numel(theta)],'nonempty','real'},'patternCustom');
    end
    
else
    if isa(MagE,'phased.CustomAntennaElement')||isa(MagE,'phased.IsotropicAntennaElement')
        narginchk(1,7);
        if isa(MagE, 'phased.CustomAntennaElement')
            freqpat = sum(MagE.FrequencyVector)/numel(MagE.FrequencyVector);
        else
            freqpat = sum(MagE.FrequencyRange)/2;
        end
        [MagE,az,el] = pattern(MagE,freqpat); % Dummy frequency
        MagE = MagE';
        phi = az';
        theta = 90-el;
        typeValidationMat = @(x) validateattributes(x,{'numeric'},          ...
            {'size',[numel(phi),numel(theta)],'nonempty','real'},'patternCustom');
    else
        error(message('shared_channel:patternCustom:UnsupportedObject'));
    end    
end

expectedcoord = {'polar','rectangular'};
expectedslice = {'phi','theta'};
addRequired(parseobj,'MagE',typeValidationMat);
typeValidationVec = @(x) validateattributes(x,{'numeric'},              ...
    {'vector','nonempty', 'real','finite', 'nonnan'},'patternCustom');
addRequired(parseobj,'theta', typeValidationVec);
addRequired(parseobj,'phi',typeValidationVec);
addParameter(parseobj,'CoordinateSystem','polar',                   ...
    @(x) any(validatestring(x,expectedcoord)));
addParameter(parseobj,'Slice','',                   ...
    @(x) any(validatestring(x,expectedslice)));
typeValidationSlice = @(x) validateattributes(x,{'numeric'},          ...
    {'vector','nonempty','real','finite', 'nonnan'}, 'patternCustom');
addParameter(parseobj,'SliceValue',[],typeValidationSlice);
parse(parseobj, MagE, theta, phi, varargin{:});

% Calling the function to get back the rearranged data
[MagEBlocks,theta1,phi1] = antennashared.internal.rearrangeData(MagE,theta,phi);

% Checking for the output argument
if nargout <= 1
    % This loop will be executed if 2D data is to be plotted
    if ((~isequal(parseobj.Results.Slice,'')) &&                        ...
            (~isequal(parseobj.Results.SliceValue,[])))
        
        % Calling the radpattern2Ddata function to get the data in the
        % appropriate format
        [MagStore,AngStore] = antennashared.internal.radpattern2Ddata(MagEBlocks,  ...
            phi1,theta1,1e9,parseobj.Results.Slice,                     ...
            parseobj.Results.SliceValue);
        
        % This loop is executed if the coordinatesystem specified is polar
        if strcmpi(parseobj.Results.CoordinateSystem,'polar')                              
            hPlot = polarpattern(AngStore(:,1),MagStore,                ...
                'DrawGridToOrigin', 1, 'LineWidth', 2, 'GridWidth', 1.5,...
                'AngleResolution', 30);
            if strcmpi(parseobj.Results.Slice,'theta')
                createLabels(hPlot, 'theta=%d#deg',                     ...
                    parseobj.Results.SliceValue(1:end));
            elseif strcmpi(parseobj.Results.Slice,'phi')
                createLabels(hPlot, 'phi=%d#deg',                       ...
                    parseobj.Results.SliceValue(1:end));
            end       
        % coordinatesystem specified is rectangular
        elseif strcmpi(parseobj.Results.CoordinateSystem,'rectangular')
            str = cell(length(parseobj.Results.SliceValue),1);            
            hPlot = plot(gca,AngStore,MagStore);
            set(hPlot,'LineWidth',2);
            grid on;
            ylabel('Magnitude');
            if strcmpi(parseobj.Results.Slice,'theta')
                xlabel('Phi (degree)');
                for m = 1:length(parseobj.Results.SliceValue)
                    str{m}= strcat('theta =', num2str(parseobj.Results.SliceValue(m)),' deg');
                end
            elseif strcmpi(parseobj.Results.Slice,'phi')
                xlabel('Theta (degree)');
                for m = 1:length(parseobj.Results.SliceValue)
                    str{m}= strcat('phi =', num2str(parseobj.Results.SliceValue(m)),' deg');
                end
            end
            legend(str);
            legend('Location','Best');
        end
        
        % This loop will be executed if 3D data is to be plotted
    elseif ((isequal(parseobj.Results.Slice,'')) &&                     ...
            (isequal(parseobj.Results.SliceValue,[])))
        
        % This loop is executed if the coordinatesystem specified is polar
        if strcmpi(parseobj.Results.CoordinateSystem,'polar')
            [hPlot] = radiationpattern3DEditted(MagEBlocks,theta1,phi1,'CurrentAxes', 1);
            axesHand = hPlot.Parent;
%             axis(axesHand,'normal') % I removed
            axesHand.DataAspectRatio = [1 1 1];
            
            % This loop is executed if the coordinatesystem specified is
            % rectangular
        elseif strcmpi(parseobj.Results.CoordinateSystem,'rectangular')
            [thetaBlocks,phiBlocks] = meshgrid(theta1, phi1);
            [hPlot,axesHand] = antennashared.internal.RectangularPlot3D(thetaBlocks, ...
                phiBlocks,MagEBlocks, 1);
            axis(axesHand,'normal')
            axesHand.DataAspectRatio = [1 1 1];
            %set(axesHand,'Position',[0.1300 0.1100 0.7750 0.8150]);
            xlabel('Theta (degree)');
            ylabel('Phi (degree)');
            zlabel('Magnitude');
        end
        
        % These loops are executed if the number of input arguments are not enough
    elseif ((isequal(parseobj.Results.Slice,'')) && (~isequal(parseobj.Results.SliceValue,[])))
        error(message('shared_channel:patternCustom:UnspecifiedCut'));
    elseif (~(isequal(parseobj.Results.Slice,'')) && (isequal(parseobj.Results.SliceValue,[])))
        error(message('shared_channel:patternCustom:UnspecifiedCutValue'));
    end
    
    if nargout == 1
        varargout{1} = hPlot;
    end
elseif nargout > 1
    error(message('shared_channel:patternCustom:IncorrectNumArguments','output','output','1'));
end
end