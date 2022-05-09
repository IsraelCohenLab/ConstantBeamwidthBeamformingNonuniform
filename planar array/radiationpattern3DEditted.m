function varargout =  radiationpattern3DEditted(MagE1,theta1,phi1,varargin) 
%RADIATIONPATTERN3D plots 3D radiation pattern for the antenna
%
%   The function plots the 3D radiation pattern of the antenna 
%
%   OBJ is the antenna object
%
%   FREQ is a single frequency point at which the radiation pattern is
%   plotted

%   Copyright 2015-2020 The MathWorks, Inc.

% Compute the max and min values of directivity
minval = min(min(MagE1));

parserObj          = inputParser();
addParameter(parserObj,'offset',minval);
addParameter(parserObj,'spherical',false);
addParameter(parserObj,'CurrentAxes',0);
addParameter(parserObj,'plottype','');

parse(parserObj,varargin{:});
temp = parserObj.Results.CurrentAxes;
minval = parserObj.Results.offset;
if(parserObj.Results.spherical)
    r = ones(length(phi1),length(theta1));
elseif isequal(unique(MagE1),minval)        % The isotropic case
    r = MagE1;
else
    r = MagE1 - minval;
end

[~, units] = em.FieldAnalysisWithFeed.getfieldlabels(parserObj.Results.plottype);

if ~temp
    % Create the axis for the radiation pattern
    if ~isempty(get(groot,'CurrentFigure'))
        clf(gcf);
    end
    haxRadPat = axes('Parent', gcf, 'Position', [0.28 0.24 0.7 0.7]);
    
    hfig = ancestor(haxRadPat,'figure'); %hfig = get(haxRadPat, 'Parent');
    try
        set(hfig,'toolbar','figure'); % Show standard toolbar %none
    catch
    end
%     cameratoolbar(hfig,'hide'); % Hide camera toolbar    
    tempval = 0; 
else
    haxRadPat = gca;
    hfig = ancestor(haxRadPat,'figure');
    tempval = 1;
end

hold on;
add_x_y_z_labels(haxRadPat);
add_az_el_labels();
draw_circle(90, 1:1:360, haxRadPat,'b'); % Circle in the x-y plane
draw_circle(1:1:360, 0, haxRadPat, 'g'); % Circle in the x-z plane
draw_circle(1:1:360,90, haxRadPat, 'r'); % Circle in the y-z plane
surfHdl = draw_3d_plot(MagE1,theta1,phi1, r, haxRadPat, tempval);

z = zoom;
z.setAxes3DPanAndZoomStyle(haxRadPat,'camera');

% This is for patternCustom function.
if nargout == 1
    varargout{1} = surfHdl;
end

% Add colorbar
% Till geck G1115830 is fixed

cbar = colorbar('peer',haxRadPat);
ylabel(cbar, units);
view(haxRadPat,135,20);
hold off;
try
    set(hfig, 'NextPlot', 'replace');
    if antennashared.internal.figureForwardState(hfig) 
        shg;
    end
catch
end

end% of radiationpattern3D

function surfHdl =  draw_3d_plot(MagE1,theta1,phi1,r1, axes1, val)

[theta,phi] = meshgrid(theta1, phi1);
MagE        = reshape(MagE1,length(phi1),length(theta1));
r           = reshape(r1,length(phi1),length(theta1));
[X, Y, Z]   =  antennashared.internal.sph2cart(phi, theta, r./max(max(r)));

surfHdl = surf(axes1,X,Y,Z,MagE, 'FaceColor','interp');
set(surfHdl,'LineStyle','none','FaceAlpha',1.0,'Tag','3D polar plot');% change to 0.9
if ~val   
    set(axes1, 'Position',[0.2 0.2 0.7 0.7], 'DataAspectRatio',[1 1 1]); %[0.28 0.24 0.7 0.7] %[0.05 0.05 0.9 .85]
    axis vis3d
    axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);
    axis off
else
    %set(axes1,'DataAspectRatio',[1 1 1]);
    axis vis3d    
    axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);
    axis off
    axis equal
end
colormap(jet(256));

end% of draw_3d_plot

function draw_circle (theta, phi, axes1, color)

[theta,phi] = meshgrid(theta, phi);
[X, Y, Z]   = antennashared.internal.sph2cart(phi, theta, 1.1);
plot3(axes1,X,Y,Z,color,'LineWidth',2);

end % of draw_circle

function add_x_y_z_labels(axes1)
% Create pseudo-axes and x/y/z mark ticks

r      = 1.2;
XPos   = r;
YPos   = r;
ZPos   = r;
plot3( axes1, [0,XPos],[0,0],[0,0],'r','LineWidth',1.5 );
text(1.2*XPos,0,0, '$x$','Interpreter','latex','FontName','Times New Roman','HorizontalAlignment','center');
plot3( axes1, [0,0],[0,YPos],[0,0],'g','LineWidth',1.5 );
text(0,1.2*YPos,0, '$y$','Interpreter','latex','FontName','Times New Roman','HorizontalAlignment','center');
plot3( axes1, [0,0],[0,0],[0,ZPos],'b','LineWidth',1.5 );
text(0,0,1.15*ZPos, '$z$','Interpreter','latex','FontName','Times New Roman','HorizontalAlignment','center');

% figures
global title_text
text(0,0,ZPos+0.5, title_text,'HorizontalAlignment','center','Interpreter','latex','FontName','Times New Roman');%,'FontSize',hcfontsize,'FontName','Times New Roman')

end% of add_x_y_z_labels

function add_az_el_labels()
% Display azimuth/elevation

% Create arrows to show azimuth and elevation variation
XPos = 1.15;
ZPos = 1.15;
% draw_arrow([XPos 0],[XPos 0.1],1.5);
% text(1.2,0.12,0.0, texlabel('az'));
% draw_arrow([XPos 0],[XPos 0.1],1.5,  0, 'xz');
% text(1.2,-0.025,0.15, texlabel('el'));
% draw_arrow([0 ZPos],[0.1 ZPos],1.5, 0, 'xz');
% text(0.16,0,1.2, texlabel('el'));

end% of add_az_el_labels

function draw_arrow(startpoint,endpoint,headsize, offset, plane)

if nargin == 3
    plane = 'xy';
    offset= 0;
end

v1 = headsize*(startpoint-endpoint)/2.5;

theta      = 22.5*pi/180;
theta1     = -1*22.5*pi/180;
rotMatrix  = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
rotMatrix1 = [cos(theta1) -sin(theta1) ; sin(theta1) cos(theta1)];

v2 = v1*rotMatrix;
v3 = v1*rotMatrix1;
x1 = endpoint;
x2 = x1 + v2;
x3 = x1 + v3;
if strcmpi(plane, 'xy')
    fill3([x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],[offset offset offset],'k');    
    plot3([startpoint(1) endpoint(1)],[startpoint(2) endpoint(2)],      ...
        [offset offset],'linewidth',1.5,'color','k');
elseif strcmpi(plane,'xz')    
    fill3([x1(1) x2(1) x3(1)],[offset offset offset],[x1(2) x2(2) x3(2)],'k');    
    plot3([startpoint(1) endpoint(1)],[offset offset],                  ...
        [startpoint(2) endpoint(2)],'linewidth',1.5,'color','k');
elseif strcmpi(plane,'yz')    
    fill3([offset offset offset],[x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],'k');    
    plot3([offset offset],[startpoint(1) endpoint(1)],                  ...
        [startpoint(2) endpoint(2)],'linewidth',1.5,'color','k');
end
end% of draw arrow