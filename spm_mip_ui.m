function varargout = spm_mip_ui(varargin)
% GUI for displaying MIPs with interactive pointers
% FORMAT hMIPax = spm_mip_ui(Z,XYZ,M,DIM,F,units)
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% M       - voxels - > mm matrix
% DIM     - image dimensions {voxels}
% F       - Figure (or axes) to work in [Defaults to gcf]
% hMIPax  - handle of MIP axes
% units   - units of space
%
% FORMAT xyz = spm_mip_ui('GetCoords',h)
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
% xyz     - Current Talairach coordinates of cursor
%
% FORMAT [xyz,d] = spm_mip_ui('SetCoords',xyz,h,hC)
% xyz     - (Input) {3 x 1} vector of desired Talairach coordinates
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
% hC      - Handle of calling object, if used as a callback. [Default 0]
% xyz     - (Output) {3 x 1} vector of voxel centre nearest desired xyz co-ords
% d       - Euclidian distance from desired co-ords & nearest voxel
%__________________________________________________________________________
%
% spm_mip_ui displays a maximum intensity projection (using spm_mip)
% with draggable cursors.
%
% See spm_mip.m for details of MIP construction, display, and the brain
% outlines used.
%                           ----------------
%
% The cursor can be dragged to new locations in three ways:
%
% (1) Point & drop: Using the primary "select" mouse button, click on a
%     cursor and drag the crosshair which appears to the desired location.
%     On dropping, the cursors jump to the voxel centre nearest the drop
%     site.
%
% (2) Dynamic drag & drop: Using the middle "extend" mouse button, click on
%     a cursor and drag it about. The cursor follows the mouse, jumping to
%     the voxel centre nearest the pointer. A dynamically updating
%     information line appears above the MIP and gives the current
%     co-ordinates. If the current voxel centre is in the XYZ pointlist,
%     then the corresponding image value is also printed.
%
% (3) Magnetic drag & drop: As with "Dynamic drag & drop", except the cursors
%     jump to the voxel centre in the pointlist nearest to the cursor. Use
%     the right "alt" mouse button for "magnetic drag & drop".
%
% In addition a ContextMenu is provided, giving the option to jump the
% cursors to the nearest suprathreshold voxel, the nearest local
% maximum, or to the global maximum. (Right click on the MIP to bring up
% the ContextMenu.) A message in the MATLAB command window describes the
% jump.
%
%                           ----------------
%
% The current cursor position (constrained to lie on a voxel) can be
% obtained by xyz=spm_mip_ui('GetCoords',hMIPax), and set with
% xyz=spm_mip_ui('SetCoords',xyz,hMIPax), where hMIPax is the handle of
% the MIP axes, or of the figure containing a single MIP [default gcf].
% The latter rounds xyz to the nearest voxel center, returning the
% result.
%
% spm_mip_ui handles all the callbacks required for moving the cursors, and
% is "registry" enabled (See spm_XYZreg.m). Programmers help is below in the
% main body of the function.
%
%__________________________________________________________________________
% Copyright (C) 1996-2015 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_mip_ui.m 6445 2015-05-21 17:38:59Z guillaume $


%==========================================================================
% - FORMAT specifications for embedded CallBack functions
%==========================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take.          )
%
% FORMAT hMIPax = spm_mip_ui(Z,XYZ,V,F,units)
% [ShortCut] Defaults to hMIPax=spm_mip_ui('Display',Z,XYZ,V,F,units)
%
% FORMAT hMIPax = spm_mip_ui('Display',Z,XYZ,M,DIM,F,units)
% Displays the MIP and sets up cursors
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% M       - voxels - > mm matrix
% DIM     - image dimensions {voxels}
% F       - Handle of figure (or axes) to work in [Defaults to gcf]
% hMIPax  - handle of MIP axes
% units   - units of space
%
% FORMAT xyz = spm_mip_ui('GetCoords',h)
% Returns coordinates of current cursor position
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
% xyz     - Current Talairach coordinates of cursor
%
% FORMAT [xyz,d] = spm_mip_ui('SetCoords',xyz,h,hC)
% Sets cursor position
% xyz     - (Input) {3 x 1} vector of desired Talairach coordinates
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
% hC      - Handle of calling object, if used as a callback. [Default 0]
% xyz     - (Output) {3 x 1} vector of voxel centre nearest desired xyz co-ords
% d       - Euclidian distance from desired co-ords & nearest voxel
%
% FORMAT spm_mip_ui('PosnMarkerPoints',xyz,h)
% Utility routine: Positions cursor markers
% xyz     - {3 x 1} vector of Talairach coordinates for cursor
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
%
% FORMAT [xyz,d] = spm_mip_ui('Jump',h,loc)
% Utility routine (CallBack of UIcontextMenu) to jump cursor
% h       - Handle of MIP axes [defaults to spm_mip_ui('FindMIPax')]
% loc     - String defining jump: 'dntmv' - don't move
%                                 'nrvox' - nearest suprathreshold voxel
%                                 'nrmax' - nearest local maximum
%                                 'glmax' - global maximum
% xyz     - co-ordinates of voxel centre jumped to {3 x 1} vector
% d       - (square) Euclidian distance jumped
%
% FORMAT hMIPax = spm_mip_ui('FindMIPax',h)
% Looks for / checks MIP axes 'Tag'ged 'hMIPax'... errors if no valid MIP axes
% h       - Handle of MIP axes, or figure containing MIP axis [default gcf]
% hMIPax  - Handle of valid MIP axis found (errors if multiple found)
%
% FORMAT spm_mip_ui('MoveStart')
% Utility routine: CallBack for starting cursor dragging
% This is the ButtonDownFcn for the cursor markers
%
% FORMAT spm_mip_ui('Move',DragType)
% Utility routine: Initiate cursor move
% DragType - 0 = Point & drop (no dragging)
%            1 = Drag'n'drop with dynamic coordinate/value updating
%            2 = Magnetic drag'n'drop with dynamic coordinate updating
%
% FORMAT spm_mip_ui('MoveEnd')
% Utility routine: End cursor move
%__________________________________________________________________________

%-Condition arguments
%==========================================================================
if nargin==0
    error('Insufficient arguments')
elseif ~ischar(varargin{1})
    varargout={spm_mip_ui('Display',varargin{1:end})}; return
end


%-Axis offsets for 3d MIPs:
%==========================================================================
%-MIP pane dimensions, origin offsets and #pixels per mm
%-See spm_project.c for derivation
mipmat = char(spm_get_defaults('stats.results.mipmat'));
load(mipmat, 'DXYZ', 'CXYZ', 'scale');
% DMIP = [DXYZ(2)+DXYZ(1), DXYZ(1)+DXYZ(3)];
%-Coordinates of Talairach origin in multipane MIP image (Axes are 'ij' + rot90)
% Transverse: [Po(1), Po(2)]
% Sagittal  : [Po(1), Po(3)]
% Coronal   : [Po(4), Po(3)]
% 4 voxel offsets in Y since using character '<' as a pointer.
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) +2;
Po(3)  = DXYZ(3)         -CXYZ(3) +2;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;



%==========================================================================
switch lower(varargin{1}), case 'display'
    %==========================================================================
    % hMIPax = spm_mip_ui('Display',Z,XYZ,M,DIM,F,units)
    if nargin<5
        F      = gcf;
        hMIPax = [];
    else
        F = varargin{6};
        if ischar(F), F=spm_figure('FindWin',F); end
        if ~ishandle(F), error('Invalid handle'), end
        switch get(F,'Type'), case 'figure'
            hMIPax = [];
            case 'axes'
                hMIPax = F;
                F      = get(hMIPax,'Parent');
            otherwise
                error('F not a figure or axis handle')
        end
    end
    if nargin<4, error('Insufficient arguments'), end
    Z       = varargin{2};
    XYZ     = varargin{3};
    M       = varargin{4};
    DIM     = varargin{5};
    try
        units = varargin{7};
    catch
        units = {'mm' 'mm' 'mm'};
    end
    
    xyz = spm_XYZreg('RoundCoords',[0;0;0],M,DIM);
    
    %-Display (MIP) transformation matrix
    %----------------------------------------------------------------------
    % Md maps various non-spatial data
    % Ms maps spatial (mm) to pixel coordinates, it must not change
    % mapping of non-spatial coordinates
    % Ms is left-multiplied to Md
    if isequal(units,{'mm' 'mm' 'mm'})
        Md      = eye(4);
        Ms      = diag([scale(1:3) 1]);
    elseif isequal(units,{'mm' 'mm' ''})
        Md      = eye(4);
        Md(3,4) = -100;  % out of field of view (MNI)
        Ms      = diag([scale(1:2) 1 1]);
    elseif isequal(units,{'mm' 'mm' 'ms'}) || isequal(units,{'mm' 'mm' 'Hz'})
        Md      = eye(4);
        Md(3,3) = 80 / (M(3,3)*DIM(3));
        Md(3,4) = -80 * M(3,4) / (M(3,3)*DIM(3));
        Ms      = diag([scale(1:2) 1 1]);
    elseif isequal(units,{'Hz' 'ms' ''})
        Md      = eye(4);
        Md(1,1) = -136 / (M(1,1)*DIM(1));
        Md(1,4) = M(1,4)*136 / (M(1,1)*DIM(1)) + 68;
        Md(2,2) = 172 / (M(2,2)*DIM(2));
        Md(2,4) = -M(2,4)*172 / (M(2,2)*DIM(2)) - 100;
        Md(3,4) = -100;
        Ms      = eye(4);
    elseif isequal(units,{'Hz' 'Hz' ''})
        Md      = eye(4);
        Md(1,1) = 136 / (M(1,1)*DIM(1));
        Md(1,4) = -M(1,4)*136 / (M(1,1)*DIM(1)) - 68;
        Md(2,2) = 172 / (M(2,2)*DIM(2));
        Md(2,4) = -M(2,4)*172 / (M(2,2)*DIM(2)) - 100;
        Md(3,4) = -100;
        Ms      = eye(4);
    elseif isequal(units,{'mm' 'mm' '%'})
        warning('Handling of data units changed: please re-estimate model.');
        units   = {'mm' 'mm' 'ms'};
        Md      = eye(4);
        Md(3,3) = 80 / (M(3,3)*DIM(3));
        Md(3,4) = -80 * M(3,4) / (M(3,3)*DIM(3));
        Ms      = diag([scale(1:2) 1 1]);
    else
        error('Unknown data units.');
    end
    
    %-Display MIP
    %----------------------------------------------------------------------
    Funits = get(F,'Units');
    set(F,'Units','normalized')
    if isempty(hMIPax)
        hMIPax = axes('Position',[0.24 0.54 0.62 0.42],'Parent',F);
    else
        axes(hMIPax), cla reset
    end
    
    %-NB: spm_mip's `image` uses a newplot, & screws stuff without the figure.
    figure(F)
    pXYZ = Ms*Md*[XYZ;ones(1,size(XYZ,2))];
    spm_mip(Z,pXYZ(1:3,:),Ms*Md*M,units);
    hMIPim = get(gca,'Children');
    
    
    %-Print coordinates
    %----------------------------------------------------------------------
    hMIPxyz = text(0,max(get(hMIPax,'YLim'))/2,...
        {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz(1:3))},...
        'Interpreter','TeX','FontName',spm_platform('font','times'),...
        'Color',[1,1,1]*.7,...
        'HorizontalAlignment','Center',...
        'VerticalAlignment','Bottom',...
        'Rotation',90,...
        'Tag','hMIPxyz',...
        'UserData',xyz);
    
    %-Create point markers
    %----------------------------------------------------------------------
    xyz = Ms*Md*[xyz(:);1];
    xyz = xyz(1:3);
    hX1r  = text(Po(1)+xyz(2),Po(2)+xyz(1),'<',...
        'Color','r','Fontsize',16,...
        'Tag','hX1r',...
        'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
    hX2r  = text(Po(1)+xyz(2),Po(3)-xyz(3),'<',...
        'Color','r','Fontsize',16,...
        'Tag','hX2r',...
        'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
    hX3r  = text(Po(4)+xyz(1),Po(3)-xyz(3),'<',...
        'Color','r','Fontsize',16,...
        'Tag','hX3r',...
        'ButtonDownFcn','spm_mip_ui(''MoveStart'')');
    hXr   = [hX1r,hX2r,hX3r];
    
    
    if DIM(3) == 1
        %-2 dimensional data
        %------------------------------------------------------------------
        set(hXr(3),'Visible','off');
        set(hXr(2),'Visible','off');
        
    end
    
    %-Create UIContextMenu for marker jumping
    %-----------------------------------------------------------------------
    h = uicontextmenu('Tag','MIPconmen','UserData',hMIPax,'Parent',F);
    if isempty(XYZ), str='off'; else str='on'; end
    uimenu(h,'Separator','off','Label','goto nearest suprathreshold voxel',...
        'CallBack',['spm_mip_ui(''Jump'',',...
        'get(get(gcbo,''Parent''),''UserData''),''nrvox'');'],...
        'Interruptible','off','BusyAction','Cancel','Enable',str);
    uimenu(h,'Separator','off','Label','goto nearest local maximum',...
        'CallBack',['spm_mip_ui(''Jump'',',...
        'get(get(gcbo,''Parent''),''UserData''),''nrmax'');'],...
        'Interruptible','off','BusyAction','Cancel','Enable',str);
    uimenu(h,'Separator','off','Label','goto global maximum',...
        'CallBack',['spm_mip_ui(''Jump'',',...
        'get(get(gcbo,''Parent''),''UserData''),''glmax'');'],...
        'Interruptible','off','BusyAction','Cancel','Enable',str);
    % ====Edition adapted by Changming Chen starts here,20181001====
    try
        SPM  = evalin('base','SPM');
        h1=uimenu(h,'Separator','off','Label','Neural-Behavior Scattering and Fit voxel');
        if isfield(SPM.xsDes,'Basis_functions') && strcmp(SPM.xsDes.Basis_functions,'hrf')
            uimenu(h1,'Label','VOXEL: plot BOLD signal against design matrix',...
                'CallBack',{@neubehfit,'voxel','cc'},...
                'Interruptible','off','BusyAction','Cancel');
        elseif strcmp(SPM.xsDes.Design,'Multiple regression')
            uimenu(h1,'Label','VOXEL: Behavioral data: select the index of interest in design matrix, regression supported)',...
                'CallBack',{@neubehfit,'voxel','xx'},...
                'Interruptible','off','BusyAction','Cancel');
            uimenu(h1,'Label','VOXEL: Behavioral data (variables from workspace, simple/multiple regression supported)',...
                'CallBack',{@neubehfit,'voxel','yy'},...
                'Interruptible','off','BusyAction','Cancel');
        elseif strcmp(SPM.xsDes.Design,'Full factorial') ||  strcmp(SPM.xsDes.Design,'Two-sample t-test')
            uimenu(h1,'Label','VOXEL: Plot individual datapoints',...
                'CallBack',{@neubehfit,'voxel','zz'},...
                'Interruptible','off','BusyAction','Cancel');
            uimenu(h1,'Label','VOXEL: Behavioral data (variables from workspace, simple/multiple regression supported)',...
                'CallBack',{@neubehfit,'voxel','kk'},...
                'Interruptible','off','BusyAction','Cancel');
        end
        
        h2=uimenu(h,'Separator','off','Label','Neu-Beh Scattering Fit Cluster',...
            'CallBack',['spm_mip_ui(''Jump'',',...
            'get(get(gcbo,''Parent''),''UserData''),''glmax'');'],...
            'Interruptible','off','BusyAction','Cancel','Enable',str);
        if isfield(SPM.xsDes,'Basis_functions') && strcmp(SPM.xsDes.Basis_functions,'hrf')
            uimenu(h2,'Label','Cluster: plot BOLD signal against design matrix',...
                'CallBack',{@neubehfit,'cluster','cc'},...
                'Interruptible','off','BusyAction','Cancel');
        elseif strcmp(SPM.xsDes.Design,'Multiple regression')
            uimenu(h2,'Label','Cluster: Behavioral data: indices of interest in design matrix, regression supported)',...
                'CallBack',{@neubehfit,'cluster','xx'},...
                'Interruptible','off','BusyAction','Cancel');
            uimenu(h2,'Label','Cluster: Behavioral data: variables from workspace, simple/multiple regression supported)',...
                'CallBack',{@neubehfit,'cluster','yy'},...
                'Interruptible','off','BusyAction','Cancel');
        elseif strcmp(SPM.xsDes.Design,'Full factorial') ||  strcmp(SPM.xsDes.Design,'Two-sample t-test')
            uimenu(h2,'Label','Cluster: Plot individual datapoints',...
                'CallBack',{@neubehfit,'cluster','zz'},...
                'Interruptible','off','BusyAction','Cancel');
            uimenu(h2,'Label','Cluster: Behavioral data: variables from workspace, simple/multiple regression supported)',...
                'CallBack',{@neubehfit,'cluster','kk'},...
                'Interruptible','off','BusyAction','Cancel');
        end
    catch
    end
    % ==Edition adapted by Changming Chen finished here====
    
    uimenu(h,'Separator','on','Label','save MIP as...',...
        'CallBack',['spm_mip_ui(''Save'', ',...
        'get(get(gcbo,''Parent''),''UserData''));'],...
        'Interruptible','off','BusyAction','Cancel');
    h1 = uimenu(h,'Separator','on','Label','Extract data');
    h2 = uimenu(h1,'Label','raw y');
    uimenu(h2,'Label','This voxel',...
        'CallBack','y=spm_mip_ui(''Extract'', ''Y'', ''voxel'')',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h2,'Label','This cluster',...
        'CallBack','y=spm_mip_ui(''Extract'', ''Y'', ''cluster'')',...
        'Interruptible','off','BusyAction','Cancel');
    h2 = uimenu(h1,'Label','whitened and filtered y');
    uimenu(h2,'Label','This voxel',...
        'CallBack','y=spm_mip_ui(''Extract'', ''KWY'', ''voxel'')',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h2,'Label','This cluster',...
        'CallBack','y=spm_mip_ui(''Extract'', ''KWY'', ''cluster'')',...
        'Interruptible','off','BusyAction','Cancel');
    h2 = uimenu(h1,'Label','betas');
    uimenu(h2,'Label','This voxel',...
        'CallBack','beta=spm_mip_ui(''Extract'', ''beta'', ''voxel'')',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h2,'Label','This cluster',...
        'CallBack','beta=spm_mip_ui(''Extract'', ''beta'', ''cluster'')',...
        'Interruptible','off','BusyAction','Cancel');
    
    % overlay sensor positions for M/EEG
    %----------------------------------------------------------------------
    if strcmp(spm('CheckModality'), 'EEG') && ...
            (isequal(units,{'mm' 'mm' 'ms'}) || isequal(units,{'mm' 'mm' 'Hz'}))
        uimenu(h,'Separator','on','Label','display/hide channels',...
            'CallBack',['spm_mip_ui(''Channels'', ',...
            'get(get(gcbo,''Parent''),''UserData''));'],...
            'Interruptible','off','BusyAction','Cancel','Enable',str);
        uimenu(h,'Separator','off','Label','goto nearest suprathreshold channel',...
            'CallBack',['spm_mip_ui(''Jump'',',...
            'get(get(gcbo,''Parent''),''UserData''),''nrchan'');'],...
            'Interruptible','off','BusyAction','Cancel','Enable',str);
    end
    
    set(hMIPim,'UIContextMenu',h)
    
    %-Save handles and data
    %----------------------------------------------------------------------
    set(hMIPax,'Tag','hMIPax','UserData',...
        struct(...
        'hReg',     [],...
        'hMIPxyz',  hMIPxyz,...
        'XYZ',      XYZ,...
        'Z',        Z,...
        'M',        M,...
        'Md',       Md,...
        'Ms',       Ms,...
        'DIM',      DIM,...
        'units',    {units},...
        'hXr',      hXr))
    
    varargout = {hMIPax};
    
    
    
    %======================================================================
    case 'getcoords'
        %======================================================================
        % xyz = spm_mip_ui('GetCoords',h)
        if nargin<2, h=spm_mip_ui('FindMIPax'); else h=varargin{2}; end
        varargout = {get(getfield(get(h,'UserData'),'hMIPxyz')  ,'UserData')};
        
        
        
        %======================================================================
    case 'setcoords'
        %======================================================================
        % [xyz,d] = spm_mip_ui('SetCoords',xyz,h,hC)
        if nargin<4, hC=0; else hC=varargin{4}; end
        if nargin<3, h=spm_mip_ui('FindMIPax'); else h=varargin{3}; end
        if nargin<2, error('Set coords to what?'), else xyz=varargin{2}; end
        
        MD  = get(h,'UserData');
        
        %-Check validity of coords only when called without a caller handle
        %------------------------------------------------------------------
        if hC<=0
            [xyz,d] = spm_XYZreg('RoundCoords',xyz,MD.M,MD.DIM);
            if d>0 && nargout<2, warning(sprintf(...
                    '%s: Co-ords rounded to nearest voxel center: Discrepancy %.2f',...
                    mfilename,d));
            end
        else
            d = [];
        end
        
        %-Move marker points, update internal cache in hMIPxyz
        %------------------------------------------------------------------
        spm_mip_ui('PosnMarkerPoints',MD.Md(1:3,:)*[xyz;1],h);
        set(MD.hMIPxyz,'UserData',reshape(xyz(1:3),3,1))
        set(MD.hMIPxyz,'String',{'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz(1:3))})
        
        %-Tell the registry, if we've not been called by the registry...
        %------------------------------------------------------------------
        if ~isempty(MD.hReg) && MD.hReg~=hC, spm_XYZreg('SetCoords',xyz,MD.hReg,h); end
        
        %-Return arguments
        %------------------------------------------------------------------
        varargout = {xyz,d};
        
        
        
        %======================================================================
    case 'posnmarkerpoints'
        %======================================================================
        % spm_mip_ui('PosnMarkerPoints',xyz,h)
        if nargin<3, h=spm_mip_ui('FindMIPax'); else h=varargin{3}; end
        if nargin<2, xyz = spm_mip_ui('GetCoords',h); else xyz = varargin{2}; end
        MD = get(h,'UserData');
        
        %-Get handles of marker points from UserData of hMIPax
        %------------------------------------------------------------------
        hX = MD.hXr;
        
        %-Set marker points
        %------------------------------------------------------------------
        set(hX,'Units','Data')
        if length(hX)==1
            vx  = sqrt(sum(MD.M(1:3,1:3).^2));
            tmp = MD.M\[xyz ; 1];
            tmp = tmp(1:2).*vx(1:2)';
            set(hX,'Position',[tmp(1), tmp(2), 1])
        else
            pxyz = MD.Ms*[xyz(1:3);1];
            set(hX(1),'Position',[ Po(1) + pxyz(2), Po(2) + pxyz(1), 0])
            set(hX(2),'Position',[ Po(1) + pxyz(2), Po(3) - pxyz(3), 0])
            set(hX(3),'Position',[ Po(4) + pxyz(1), Po(3) - pxyz(3), 0])
        end
        
        
        %======================================================================
    case 'jump'
        %======================================================================
        % [xyz,d] = spm_mip_ui('Jump',h,loc)
        if nargin<3, loc='nrvox'; else loc=varargin{3}; end
        if nargin<2, h=spm_mip_ui('FindMIPax'); else h=varargin{2}; end
        
        %-Get current location & MipData
        %------------------------------------------------------------------
        oxyz = spm_mip_ui('GetCoords',h);
        MD   = get(h,'UserData');
        
        
        %-Compute location to jump to
        %------------------------------------------------------------------
        if isempty(MD.XYZ), loc='dntmv'; end
        switch lower(loc), case 'dntmv'
            spm('alert!','No suprathreshold voxels to jump to!',mfilename,0);
            varargout = {oxyz, 0};
            return
            case 'nrvox'
                str       = 'nearest suprathreshold voxel';
                [xyz,i,d] = spm_XYZreg('NearestXYZ',oxyz,MD.XYZ);
            case 'nrmax'
                str       = 'nearest local maximum';
                iM        = inv(MD.M);
                XYZvox    = iM(1:3,:)*[MD.XYZ; ones(1,size(MD.XYZ,2))];
                [null,null,XYZvox] = spm_max(MD.Z,XYZvox);
                XYZ       = MD.M(1:3,:)*[XYZvox; ones(1,size(XYZvox,2))];
                [xyz,i,d] = spm_XYZreg('NearestXYZ',oxyz,XYZ);
            case 'glmax'
                str       = 'global maximum';
                [null, i] = max(MD.Z); i = i(1);
                xyz       = MD.XYZ(:,i);
                d         = sqrt(sum((oxyz-xyz).^2));
            case 'nrchan'
                str       = 'nearest suprathreshold channel';
                if ~isfield(MD, 'hChanPlot'), spm_mip_ui('Channels',h); end
                MD        = get(h,'UserData');
                if ~isfield(MD, 'hChanPlot'), return; end
                [xyz,i,d] = spm_XYZreg('NearestXYZ',[oxyz(1); oxyz(2); 0],MD.Channels.pos);
                xyz(3)    = oxyz(3);
                str       = [str sprintf(' (%s)',MD.Channels.name{i})];
            otherwise
                warning('Unknown jumpmode')
                varargout = {xyz,0};
                return
        end
        
        %-Write jump report, jump, and return arguments
        %------------------------------------------------------------------
        fprintf(['\n\t%s:\tJumped %0.2fmm from [%3.0f, %3.0f, %3.0f],\n\t\t\t',...
            'to %s at [%3.0f, %3.0f, %3.0f]\n'],...
            mfilename, d, oxyz, str, xyz)
        
        spm_mip_ui('SetCoords',xyz,h,h);
        varargout = {xyz, d};
        
        
        %======================================================================
    case 'findmipax'
        %======================================================================
        % hMIPax = spm_mip_ui('FindMIPax',h)
        % Checks / finds hMIPax handles
        %-**** h is handle of hMIPax, or figure containing MIP (default gcf)
        if nargin<2, h=get(0,'CurrentFigure'); else h=varargin{2}; end
        if ischar(h), h=spm_figure('FindWin',h); end
        if ~ishandle(h), error('invalid handle'), end
        if ~strcmp(get(h,'Tag'),'hMIPax'), h=findobj(h,'Tag','hMIPax'); end
        if isempty(h), error('MIP axes not found'), end
        if length(h)>1, error('Multiple MIPs in this figure'), end
        varargout = {h};
        
        
        
        %======================================================================
    case 'movestart'
        %======================================================================
        % spm_mip_ui('MoveStart')
        [cO,cF] = gcbo;
        hMIPax  = get(cO,'Parent');
        MD      = get(hMIPax,'UserData');
        
        %-Store useful quantities in UserData of gcbo, the object to be dragged
        %------------------------------------------------------------------
        set(hMIPax,'Units','Pixels')
        set(cO,'UserData',struct(...
            'hReg',     MD.hReg,...
            'xyz',      spm_mip_ui('GetCoords',hMIPax),...
            'MIPaxPos', get(hMIPax,'Position')*[1,0;0,1;0,0;0,0],...
            'hMIPxyz',  MD.hMIPxyz,...
            'M',        MD.M,...
            'Md',       MD.Md,...
            'Ms',       MD.Ms,...
            'DIM',      MD.DIM,...
            'hX',       MD.hXr))
        
        %-Initiate dragging
        %------------------------------------------------------------------
        if strcmp(get(cF,'SelectionType'),'normal') || isempty(MD.XYZ)
            %-Set Figure callbacks for drop but no drag (DragType 0)
            %--------------------------------------------------------------
            set(MD.hMIPxyz,'Visible','on','String',...
                {'{\bfSPM}{\itmip}','\itPoint & drop...'})
            set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''Move'',0)',...
                'Interruptible','off')
            set(cF,'Pointer','CrossHair')
            
            %-Set Figure callbacks for drag'n'drop (DragType 1)
            %-------------------------------------------------------------------
        elseif strcmp(get(cF,'SelectionType'),'extend')
            set(MD.hMIPxyz,'Visible','on','String',...
                {'{\bfSPM}{\itmip}','\itDynamic drag & drop...'})
            set(cF,'WindowButtonMotionFcn','spm_mip_ui(''Move'',1)',...
                'Interruptible','off')
            set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''MoveEnd'')',...
                'Interruptible','off')
            set(cF,'Pointer','Fleur')
            
            %-Set Figure callbacks for drag'n'drop with co-ord updating (DragType 2)
            %------------------------------------------------------------------
        elseif strcmp(get(cF,'SelectionType'),'alt')
            set(MD.hMIPxyz,'Visible','on','String',...
                {'{\bfSPM}{\itmip}','\itMagnetic drag & drop...'})
            set(cF,'WindowButtonMotionFcn','spm_mip_ui(''Move'',2)',...
                'Interruptible','off')
            set(cF,'WindowButtonUpFcn',    'spm_mip_ui(''MoveEnd'')',...
                'Interruptible','off')
            set(cF,'Pointer','Fleur')
        end
        
        
        
        %======================================================================
    case 'move'
        %======================================================================
        % spm_mip_ui('Move',DragType)
        if nargin<2, DragType = 2; else DragType = varargin{2}; end
        cF = gcbf;
        cO = gco(cF);
        
        %-Get useful data from UserData of gcbo, the object to be dragged
        %------------------------------------------------------------------
        MS  = get(cO,'UserData');
        
        %-Work out where we are moving to - Use HandleGraphics to give position
        %------------------------------------------------------------------
        set(cF,'Units','pixels')
        d = get(cF,'CurrentPoint') - MS.MIPaxPos;
        set(cO,'Units','pixels')
        set(cO,'Position',d)
        set(cO,'Units','data')
        d = get(cO,'Position');
        
        %-Work out xyz, depending on which view is being manipulated
        %------------------------------------------------------------------
        sMarker = get(cO,'Tag');
        pxyz = MS.Ms*[MS.xyz;1];
        if strcmp(sMarker,'hX1r')
            xyz = [d(2) - Po(2); d(1) - Po(1); pxyz(3)];
        elseif strcmp(sMarker,'hX2r')
            xyz = [pxyz(1); d(1) - Po(1); Po(3) - d(2)];
        elseif strcmp(sMarker,'hX3r')
            xyz = [d(1) - Po(4); pxyz(2); Po(3) - d(2)];
        else
            error('Can''t work out which marker point')
        end
        xyz = inv(MS.Ms*MS.Md) * [xyz(:);1]; xyz = xyz(1:3);
        
        %-Round coordinates according to DragType & set in hMIPxyz's UserData
        %------------------------------------------------------------------
        if DragType==0
            xyz    = spm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
        elseif DragType==1
            xyz    = spm_XYZreg('RoundCoords',xyz,MS.M,MS.DIM);
            hMIPax = get(cO,'Parent');
            MD     = get(hMIPax,'UserData');
            i      = spm_XYZreg('FindXYZ',xyz,MD.XYZ);
        elseif DragType==2
            hMIPax = get(cO,'Parent');
            MD     = get(hMIPax,'UserData');
            [xyz,i,d] = spm_XYZreg('NearestXYZ',xyz,MD.XYZ);
        end
        set(MS.hMIPxyz,'UserData',xyz)
        
        %-Move marker points
        %------------------------------------------------------------------
        set(MS.hX,'Units','Data')
        xyz2 = MS.Ms * MS.Md * [xyz(:);1];
        xyz2 = xyz2(1:3);
        set(MS.hX(1),'Position',[ Po(1) + xyz2(2), Po(2) + xyz2(1), 0])
        set(MS.hX(2),'Position',[ Po(1) + xyz2(2), Po(3) - xyz2(3), 0])
        set(MS.hX(3),'Position',[ Po(4) + xyz2(1), Po(3) - xyz2(3), 0])
        
        
        %-Update dynamic co-ordinate strings (if appropriate DragType)
        %------------------------------------------------------------------
        if DragType==0
            spm_mip_ui('MoveEnd')
        elseif DragType==1
            if isempty(i)
                str = {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',xyz)};
            else
                str = {'{\bfSPM}{\itmip}: ',...
                    sprintf('[%g, %g, %g]: %.4f',xyz,MD.Z(i))};
            end
            set(MD.hMIPxyz,'String',str)
        elseif DragType==2
            set(MD.hMIPxyz,'String',...
                {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]: %.4f',xyz,MD.Z(i))})
        else
            error('Illegal DragType')
        end
        
        
        
        %======================================================================
    case 'moveend'
        %======================================================================
        % spm_mip_ui('MoveEnd')
        cF = gcbf;
        cO = gco(cF);
        hMIPax  = get(cO,'Parent');
        MS      = get(cO,'UserData');
        
        %-Reset WindowButton functions, pointer & SPMmip info-string
        %------------------------------------------------------------------
        set(gcbf,'WindowButtonMotionFcn',' ')
        set(gcbf,'WindowButtonUpFcn',' ')
        set(gcbf,'Pointer','arrow')
        set(MS.hMIPxyz,'String',...
            {'{\bfSPM}{\itmip}',sprintf('[%g, %g, %g]',get(MS.hMIPxyz,'UserData'))})
        
        %-Set coordinates after drag'n'drop, tell registry
        %------------------------------------------------------------------
        % don't need to set internal coordinates 'cos 'move' does that
        if ~isempty(MS.hReg)
            spm_XYZreg('SetCoords',get(MS.hMIPxyz,'UserData'),MS.hReg,hMIPax);
        end
        
        %======================================================================
    case 'channels'
        %======================================================================
        % spm_mip_ui('Channels',h)
        % M/EEG specific to display channels on 1-slice MIP
        
        if nargin<2, h=spm_mip_ui('FindMIPax'); else h=varargin{2}; end
        MD  = get(h,'UserData');
        
        if ~isfield(MD, 'hChanPlot')
            % first time call
            D = spm_eeg_load;
            if isempty(D), return; end
            
            DIM = get(findobj('Tag','hFxyz'), 'UserData');
            
            [mod, Cind] = spm_eeg_modality_ui(D, 1, 1);
            otherind = setdiff(1:nchannels(D), Cind);
            if ~isempty(otherind)
                D = chantype(D, otherind, 'Other');
            end
            [Cel x, y] = spm_eeg_locate_channels(D, DIM.DIM(1), Cind);
            Cel = DIM.M * [Cel'; ones(2,size(Cel,1))];
            Cel = Cel(1:2,:)';
            pos = [Cel'; zeros(1,size(Cel,1))];
            Cel(:,1) = Cel(:,1) + Po(2);
            Cel(:,2) = Cel(:,2) + Po(1);
            hold(h,'on'), hChanPlot = plot(h,Cel(:, 2), Cel(:, 1), 'b*');
            
            hChanText = cell(1,size(Cel,1));
            name = D.chanlabels(Cind);
            figure(spm_figure('FindWin'));
            for i = 1:size(Cel, 1)
                hChanText{i} = text(Cel(i, 2)+0.5, Cel(i, 1), name{i}, 'Color', 'b');
            end
            
            MD.hChanPlot     = hChanPlot;
            MD.hChanText     = hChanText;
            MD.Channels.pos  = pos;
            MD.Channels.name = D.chanlabels(Cind);
            set(h, 'UserData', MD);
        else
            if strcmp(get(MD.hChanPlot, 'Visible'), 'on');
                % switch it off
                set(MD.hChanPlot, 'Visible', 'off');
                for i = 1:length(MD.hChanText)
                    set(MD.hChanText{i}, 'Visible', 'off');
                end
            else
                % switch it on
                set(MD.hChanPlot, 'Visible', 'on');
                for i = 1:length(MD.hChanText)
                    set(MD.hChanText{i}, 'Visible', 'on');
                end
            end
            
        end
        
        
        %======================================================================
    case 'save'
        %======================================================================
        % spm_mip_ui('Save',h)
        if nargin<2, h=spm_mip_ui('FindMIPax'); else h=varargin{2}; end
        %MD  = get(h,'UserData');
        %pXYZ = MD.Ms*MD.Md*[MD.XYZ;ones(1,size(MD.XYZ,2))];
        %mip  = spm_mip(MD.Z,pXYZ(1:3,:),MD.Ms*MD.Md*MD.M,MD.units);
        mip = get(findobj(h,'Type','image'),'CData');
        [f, p] = uiputfile('*.png', 'Save as', 'mip.png');
        if ~isequal(f,0) && ~isequal(p,0)
            if ndims(mip) == 3
                imwrite(mip,fullfile(p,f),'png');
            else
                imwrite(mip,gray(64),fullfile(p,f),'png');
            end
            fprintf('Saving SPM MIP in %s\n',fullfile(p,f));
        end
        
        
        %======================================================================
    case 'extract'
        %======================================================================
        % data = spm_mip_ui('Extract',what,where)
        if nargin<2, what='beta'; else what=varargin{2}; end
        if nargin<3, where='voxel'; else where=varargin{3}; end
        xSPM = evalin('base','xSPM');
        SPM  = evalin('base','SPM');
        
        XYZmm = spm_results_ui('GetCoords');
        [XYZmm,i] = spm_XYZreg('NearestXYZ',XYZmm,xSPM.XYZmm);
        if isempty(i), varargout = {[]}; return; end
        spm_results_ui('SetCoords',xSPM.XYZmm(:,i));
        
        switch lower(where)
            case 'voxel'
                % current voxel
                XYZ = SPM.xVol.iM(1:3,:)*[XYZmm;1];
            case 'cluster'
                % current cluster
                A   = spm_clusters(xSPM.XYZ);
                j   = find(A == A(i));
                XYZ = xSPM.XYZ(:,j);
            otherwise
                error('Unknown action.');
        end
        
        switch lower(what)
            case {'y','kwy'}
                data = spm_get_data(SPM.xY.VY,XYZ);
                if strcmpi(what,'kwy')
                    data = spm_filter(SPM.xX.K,SPM.xX.W*data);
                end
            case 'beta'
                data = spm_get_data(SPM.Vbeta,XYZ);
            otherwise
                error('Unknown action.');
        end
        
        varargout = {data};
        
        
        %======================================================================
    otherwise
        %======================================================================
        error('Unknown action string')
        
end
end


%=======below added by Changming Chen, for better visualization=========
function neubehfit(src,event,where,behinput)
handles=guidata(src);
global ms;
%============call extract,get y===========
xSPM = evalin('base','xSPM');
SPM  = evalin('base','SPM');
if strcmp(where, 'cluster')
    y=spm_mip_ui('Extract', 'Y', where);
    y=mean(y,2);
elseif strcmp(where, 'voxel')
    XYZmm = spm_results_ui('GetCoords');
    %     [XYZmm,~] = spm_XYZreg('RoundCoords',XYZmm,xSPM.XYZmm);
    XYZ = SPM.xVol.iM(1:3,:)*[XYZmm;1];
    try
        y = spm_get_data(SPM.xY.VY,XYZ);
    catch
        gunzip([SPM.xY.VY(1).fname,'.gz']);
        y = spm_get_data(SPM.xY.VY,XYZ);
    end
end
if find(isnan(y))
    warndlg('Some participant has image data in this location, mission aborts!');
    return;
end
%% =============================
if strcmp(behinput,'xx') || strcmp(behinput,'yy')
    if strcmp(behinput,'xx')     % regression, ask the user to input variable indices in order, and find values accordingly in design
        answer = inputdlg('Please tell me in order the indices of variables in the design matrix that you want to explore the data','Input the variable index in order, separated by one space',1);
        answer=answer{1};
        try
            temp=regexp(answer,' ','split');
            vindx=str2double(temp);
        catch
            vindx=str2double(answer);
        end
        beh=SPM.xX.X(:,vindx);
        plotregression_ccm(y,beh,vindx(1),1);
        
    elseif strcmp(behinput,'yy')
        choosevariable(y,'yy');
    end
    
elseif   strcmp(behinput,'zz')     % full-factorial
    plotfactorial_ccm(y);
    
elseif  strcmp(behinput,'kk')
    choosevariable(y,'kk');
    
elseif isfield(SPM.xsDes,'Basis_functions') && strcmp(SPM.xsDes.Basis_functions,'hrf')
    figure('tag','plot_ccm','Name','Plot Individuals Cross Conditions, Changming Chen','visible','on');
    numcons=size(SPM.Sess.U,2);
    a=jet(numcons);
    matricon=SPM.xX.X(:,(1:numcons));
    rc=range(matricon(:));
    ry=range(y);
    plot((1:length(y)),y);
    for i=1:numcons
        matricon(:,i)=matricon(:,i).*ry/rc+min(y)-0.1*ry;
        hold on;
        area(matricon(:,i),'FaceColor',a(i,:),'FaceAlpha',0.20,'EdgeColor','none');
    end
    legs={'BOLD signal'};
    for i=1:numcons
        legs(i+1)=SPM.Sess.U(i).name;
    end
    legend(legs,'EdgeColor','none','Color','none','Position',[0.80 0.82 0.16 0.12]);
    set(gca,'YLim',[min(y)-0.1*ry,max(y)+0.1*ry]);
    box off;
    ylb=ylabel('MRI signal','FontWeight','bold','tag','yl_ccm','Color',[1 0 0]);
    ylpos=get(ylb,'Position');    yticks=get(gca,'YTick'); ylpos(2)=max(yticks)-range(yticks)*0.15; set(ylb,'Position',ylpos);
    xlb=xlabel('Time','FontWeight','bold','tag','xl_ccm','Color',[1 0 0]);
    xlmt=get(gca,'XLim'); pos=get(xlb,'Position');    pos(1)=xlmt(2)-range(xlmt)*0.15;  set(xlb,'Position',pos);
end
guidata(src, handles);
end

function choosevariable(y,type)
global ms;
if strcmp(type,'yy')
    ms.type='yy';
elseif strcmp(type,'kk')
    ms.type='kk';
end
ms.y=y;
f01=findall(0,'tag','choosevariable');
if ~isempty(f01)
    close(f01);
end
figureBKcolor=[0.8314 0.8157 0.8235];
f01=figure('Units','pixels','position',[8 62 390 334],'Color',figureBKcolor,'defaultuicontrolBackgroundColor',figureBKcolor,'Name','Chooose Variable', 'NumberTitle','on','resize','off','visible','on', 'DoubleBuffer','on','MenuBar','none','tag','choosevariable','CloseRequestFcn',@figclose);
handles.f01=f01;
position.f01text01=[6 300 130 18];                 string='Where is the Variable';
handles.f01text01 = uicontrol('Parent',f01,'Style','text','Units','Pixels','position',position.f01text01,'string',string,'ForegroundColor',[0 0 0],'horizontal','left','fontsize',10);
position.f01pop01=[137 308 150 15];               string={'Workspace','Load from *.Mat file','load from *.XLS file'};
handles.f01pop01 = uicontrol('Parent',f01,'Style','popup','Units','Pixels','position',position.f01pop01,'string',string,'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],'horizontal','left','fontsize',10,'callback',@wherevars);
handles.f01=f01;
guidata(f01, handles);
end

function wherevars(hObject, eventdata)
global ms;
handles=guidata(hObject);
a=get(handles.f01pop01,'value');
varnames=cell(1);
if a==1
    % set variable settings visible
    temp=whos;
    for i=1:size(temp,1)
        varnames{i,1}=temp(i).name;
    end
elseif a==2
    % use spm_get to get file,
    [fn,pn,~] =uigetfile('*.mat','Choose the .mat file');
    temp=deblank([pn,fn]);
    t1=load(temp);
    varnames=fieldnames(t1);
    for i=1:length(varnames)
        eval(['temp=t1.',varnames{i}]);
        eval(['assignin(''base'',''',varnames{i},''',temp)']);
    end
elseif a==3
    [temp,path]=uigetfile({'*.xls;*.xlsx'},'Choose the .xls file');
    [~,filename,~]=fileparts(fullfile(path,temp));
    assignin('base',filename,eval(['xlsread(''',fullfile(path,temp),''');']));
    varnames{1}=filename;
end
handles.f01text02 = uicontrol('Parent',handles.f01,'Style','text','Units','Pixels','position',[6 260 200 18],'string','Please Select the Variable','ForegroundColor',[0 0 0],'horizontal','left','fontsize',10,'visible','on');
handles.f01pop02 = uicontrol('Parent',handles.f01,'Style','popup','Units','Pixels','position',[162 269 220 15],'string',strcat('-',varnames),'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],'fontsize',10,'horizontal','left','visible','on','callback',@fieldsorcolumn);
guidata(hObject, handles);
end

function fieldsorcolumn(hObject,eventdata)
global ms;
handles=guidata(hObject);
if isfield(handles,'f01text03')
    handles=rmfield(handles,'f01text03');
end
if isfield(handles,'f01edit03')
    handles=rmfield(handles,'f01edit03');
end
a=get(handles.f01pop02,'value');
strs=get(handles.f01pop02,'string');
idx=cellfun(@max,strfind(strs,'-'))+1;
varnames=cell(1);
for i=1:length(idx)
    varnames{i,1}=strs{i}(idx(i):end);
end
varcu=varnames{a};
varcu=evalin('base',varcu);
if isstruct(varcu)
    temp=fieldnames(varcu);
    tempfields=strcat('---',varnames{a},'.',temp);
    str=[strs(1:a-1),tempfields,strs(a+1:end)];
    handles.f01pop02 = uicontrol('Parent',handles.f01,'Style','popup','Units','Pixels','position',[162 269 220 15],'string',str,'ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],'fontsize',10,'horizontal','left','visible','on','callback',@fieldsorcolumn);
else
    handles.f01text03 = uicontrol('Parent',handles.f01,'Style','text','Units','Pixels','position',[6 220 150 15],'string','Input the column indices as the example','ForegroundColor',[0 0 0],'fontsize',10,'horizontal','left','visible','on');
    handles.f01edit03 = uicontrol('Parent',handles.f01,'Style','edit','Units','Pixels','position',[200 225 100 18],'string','[1,2]','ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],'fontsize',10,'horizontal','left','visible','on');
end
uicontrol('Parent',handles.f01,'Style','push','Units','Pixels','position',[250 10 100 18],'string','goooooo','ForegroundColor',[0 0 0],'BackgroundColor',[1 1 1],'fontsize',10,'horizontal','left','callback',@goooooo);
ms.data=varcu;
guidata(hObject, handles);
end

function goooooo(hObject,eventdata)
global ms;
y=ms.y;
if sum(double(isnan(y)))~=0
    warndlg('Some participants do not have data in this voxel/cluster, mission aborts!');
    return;
end
handles=guidata(hObject);
a=get(handles.f01edit03,'String');
eval(['data=ms.data(:,',a,')']);
beh=data;
ms.beh=beh;
assignin('base','ms',ms);
if strcmp(ms.type,'yy')
    plotregression_ccm(y,beh,1,2);
elseif strcmp(ms.type,'kk')
    plotfactorial_ccm(y,beh,1,2);
end
handles=guidata(hObject);
end

function plotregression_ccm(y,beh,bhp,type)
xSPM = evalin('base','xSPM');
SPM  = evalin('base','SPM');
h = findall(0, 'Tag', 'plot_ccm');
if ~isempty(h)
    ct=findall(h,'Tag','colortitle');
    ctc=get(ct,'Color');
    yt=findall(h,'Tag','yl_ccm');
    ytc=get(yt,'Color');
    xt=findall(h,'Tag','xl_ccm');
    xtc=get(xt,'Color');
    close(h);
end
h=figure('tag','plot_ccm','Name','Plot Individuals Cross Conditions, Changming Chen','visible','on');
% scatter(beh(:,1),y);
mdns={'purequadratic','quadratic','linear'};
for i=1:numel(mdns)
    model=mdns{i};
    stats=regstats(y,beh,model);
    text(mean(beh(:,1))+0.3*i*std(beh(:,1)),min(y)+.5*i*std(y),['p=',num2str(stats.fstat.pval),'  (',mdns{i},')']);
end
[population2, gof] = fit( beh(:,1), y, 'poly2' );
hold on
plot( population2, beh(:,1), y);
hold on;
lsline;
legend({'data','curve fit with poly2','least-squares line'},'EdgeColor',[0.5 0 0],'Color','none','Position',[0.75 0.82 0.16 0.12]);
ylb=ylabel('Neural');
ylpos=get(ylb,'Position');    yticks=get(gca,'YTick'); ylpos(2)=max(yticks)-range(yticks)*0.15; set(ylb,'Position',ylpos,'Color','r');
if type==1
    xlb=xlabel(SPM.xX.name{bhp});
else
    xlb=xlabel('First Behavior');
end
xt=get(gca,'XLim');  ylmt=get(gca,'YLim');  pos=get(xlb,'Position');    pos(1)=xt(2);  pos(2)=ylmt(1)-abs(range(ylmt))*0.07;  set(xlb,'Position',pos,'Color','r');

if size(beh,2)>1
    nx4=colorbar;
    a1=zscore(beh(:,2));
    a2=ceil(a1-min(a1))+1;
    numa2=numel(unique(a2));
    a=jet(numa2);
    for xx=1:numel(a1)
        colors(xx,:)=a(a2(xx),:);
    end
    temp=get(nx4,'Ticks');
    set(nx4,'Ticks',linspace(min(temp),max(temp),numa2),'TickLabels',linspace(min(beh(:,2)),max(beh(:,2)),numa2)');
    title(nx4,'behavioral2');
end
assignin('base','neuraldata_ccm',y);
assignin('base','behavioraldata_ccm',beh);
assignin('base','goodnessoffit_poly2_ccm',gof);
end

function plotfactorial_ccm(y,varargin)
% just individual data
xSPM = evalin('base','xSPM');
SPM  = evalin('base','SPM');
h = findall(0, 'Tag', 'plot_ccm');
if ~isempty(h)
    ct=findall(h,'Tag','colortitle');
    ctc=get(ct,'Color');
    yt=findall(h,'Tag','yl_ccm');
    ytc=get(yt,'Color');
    xt=findall(h,'Tag','xl_ccm');
    xtc=get(xt,'Color');
    close(h);
end
h=figure('tag','plot_ccm','Name','Plot Individuals Cross Conditions, Changming Chen','visible','on');
ctc=[0 0 1];
ytc=[0 0 1];
xtc=[0 0 1];
if isequal(ctc,[1 0 0]);  ctc=[0 0 1];    else      ctc=[1 0 0];     end
if isequal(ytc,[1 0 0]);  ytc=[0 0 1];    else      ytc=[1 0 0];     end
if isequal(xtc,[1 0 0]);  xtc=[0 0 1];    else      xtc=[1 0 0];     end
    
conds=SPM.xX.X*(1:size(SPM.xX.X,2))';
temp=repmat(y,1,size(SPM.xX.X,2));
temp=temp.*SPM.xX.X;
size2=size(temp,2);
for i=1:size2
    means(i)=mean(temp(SPM.xX.X(:,i)~=0,i));
    stds(i)=std(temp(SPM.xX.X(:,i)~=0,i));
    nsams(i)=numel(temp(SPM.xX.X(:,i)~=0,i));
    nsams(i)=sqrt(nsams(i));
end
bar((1:size(temp,2)),means,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'BarWidth',0.5);
hold on;
if  numel(varargin)==0
    n3=scatter(conds,y);
else
    beh=varargin{1};
    a1=zscore(beh(:,1));
    a2=ceil(a1-min(a1))+1;
    numa2=numel(unique(a2));
    a=jet(numa2);
    for xx=1:numel(a1)
        colors(xx,:)=a(a2(xx),:);
    end
    n3=scatter(conds,y,[],colors);   % attention here
    nx4=colorbar;
    clrbrtks=get(nx4,'Ticks');
    set(nx4,'Ticks',linspace(min(clrbrtks),max(clrbrtks),numa2),'TickLabels',linspace(min(beh(:,1)),max(beh(:,1)),numa2)');
    title(nx4,'behavioral');
end
hold on;
nx=errorbar((1:size(temp,2)),means,stds./nsams(i),'Color','r');
set(nx,'LineStyle','none');
box off;
ylb=ylabel('Neural','Color','r');
xlb=xlabel('Conditions','Color','r');
xticks(1:size2);
xticklabels(SPM.xX.name');
ylpos=get(ylb,'Position');    yticks=get(gca,'YTick'); ylpos(2)=max(yticks)-range(yticks)*0.15; set(ylb,'Position',ylpos,'Color','r');
xt=get(gca,'XLim');  ylmt=get(gca,'YLim');  pos=get(xlb,'Position');    pos(1)=xt(2);  pos(2)=ylmt(1)-abs(range(ylmt))*0.07;  set(xlb,'Position',pos,'Color','r');
assignin('base','neuraldata_ccm',y);
assignin('base','behavioraldata_ccm',conds);
end

function figclose(hObject,eventdata)
handles=guidata(hObject);
delete(gcf);
try
    h = findall(0, 'tag', 'plot_ccm');
    if ~isempty(h)
        close(h);
    end
catch
end
end