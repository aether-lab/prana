function [newcontrolpts] = calibration_point_select(imgmat,data)

% data.X and data.x are matrices

%scrsz = get(0,'ScreenSize');        % goes [left bottom width height]
scrsz=[1 1 1280 960];

han.hfig=figure('Position',[scrsz(3)/20 scrsz(4)/20 18*scrsz(3)/20 17.9*scrsz(4)/20]);

han.hax = axes('Parent',han.hfig,'Units','normalized','Position',[0 .333 1 .6666]);

hold(han.hax,'on');
han.himage=imshow(imgmat,'InitialMagnification','fit');
%[xdat, ydat, A] = getimage(han.himage);
%set(han.himage,'Xdata',[xdat-0.5],'Ydata',[fliplr(ydat)-0.5]);   

han.hpanel=imscrollpanel(han.hfig,han.himage);
set(han.hpanel,'Units','normalized','Position',[0 0.3333 1 0.6666]);  % position is [xmin ymin xmax ymax]

han.hMagBox = immagbox(han.hfig,han.himage);

set(han.hMagBox,'Units','normalized','Position',[0 0.965 0.06 0.035]);

%han.hpixinfo=impixelinfo(han.hfig,han.himage);
%set(han.hpixinfo,'Units','Normalized','Position',[0.75 0.355 0.1 0.02]);   

guipropsimg = guidata(han.hfig);

guipropsimg.imgmat=imgmat;
guipropsimg.imgcalpts.X=num2cell(data.X);
guipropsimg.imgcalpts.x=num2cell(data.x);
guipropsimg.markerd=data.markerd;
guipropsimg.horspace=data.horspace;
guipropsimg.verspace=data.verspace;
guipropsimg.targcolor=data.targcolor;
guipropsimg.elementsize=[];
guipropsimg.ptsselected=[];
guipropsimg.findxrange=[];
guipropsimg.findyrange=[];

% display points
guidata(han.hfig,guipropsimg);
defaultcolor=get(0,'defaultUicontrolBackgroundColor');
figcolor=get(han.hfig,'Color');     % figure background color

sztmp=size(guipropsimg.imgcalpts.X,1);
sztmp2=size(guipropsimg.imgcalpts.x,1);
tdata=cell(sztmp,5);

if ~isempty(guipropsimg.imgcalpts.X)
    tdata(:,1:2)=guipropsimg.imgcalpts.X;
end
ind=0;
if ~isempty(guipropsimg.imgcalpts.x)
    if sztmp2<=sztmp
        ind=sztmp2;
    else
        ind=sztmp;
    end

    tdata(1:ind,3:5)=guipropsimg.imgcalpts.x(1:ind,:);
end

ts=get(han.hfig,'Position');
colname={'X (pix)' 'Y (pix)' 'x (mm)' 'y (mm)' 'z (mm)'};
columneditable =  [false false true true true];
han.htable = uitable('Data',tdata,'ColumnNames',colname,... 
            'Parent',han.hfig,'Position',[ts(3)*0.63 ts(4)*0.01 ts(3)*0.36 ts(4)*0.29]);

% han.htable = uitable('Data',tdata,'ColumnNames',colname,... 
%              'Parent',han.hfig,'Units','normalized','Position',[0.5 0.01 0.36 0.29]);

han.hpanel1 = uipanel('Parent',han.hfig,'Title','','Units','normalized',...
             'Position',[.35 0.085 .12 0.19],'BackgroundColor',figcolor);

han.hpanel2 = uipanel('Parent',han.hfig,'Title','','Units','normalized',...
             'Position',[.49 0.085 .12 0.19],'BackgroundColor',figcolor);

han.htext = uicontrol(han.hfig,'Style','text','Units','normalized','String',...
    'Image Plane','Position',[0.66 0.30 0.1 0.02],'FontSize',10,...
    'FontWeight','bold','BackgroundColor',figcolor);
han.htext2 = uicontrol(han.hfig,'Style','text','Units','normalized','String',...
    'Object Plane','Position',[0.82 0.30 0.1 0.02],'FontSize',10,...
    'FontWeight','bold','BackgroundColor',figcolor);

han.hfillxyz = uicontrol(han.hpanel2,'Style','pushbutton',...
    'Units','normalized','String','Fill (x,y,z)',...
    'Position',[0.15 0.55 0.7 0.17],'Callback',{@Fillxyz_callback,han});

han.hfindXYxyz = uicontrol(han.hpanel2,'Style','pushbutton',...
    'Units','normalized','String','Find XY & xyz',...
    'Position',[0.15 0.35 0.7 0.17],'Callback',{@FindXYxyz_callback,han});

han.hfindxrange = uicontrol(han.hpanel2,'Style','edit','Units','normalized',...
    'Position',[0.06 0.22 0.75 0.1],'BackgroundColor',[1 1 1],...
    'Callback',{@Findxrange_callback,han},'String','xmin : xstep : xmax');

han.hfindyrange = uicontrol(han.hpanel2,'Style','edit','Units','normalized',...
    'Position',[0.06 0.10 0.75 0.1],'BackgroundColor',[1 1 1],...
    'Callback',{@Findyrange_callback,han},'String','ymin : ystep : ymax');

han.htextxrange = uicontrol(han.hpanel2,'Style','text','Units','normalized','String',...
    'mm','Position',[0.82 0.22 0.15 0.08],...
    'HorizontalAlignment','left','BackgroundColor',figcolor);

han.htextyrange = uicontrol(han.hpanel2,'Style','text','Units','normalized','String',...
    'mm','Position',[0.82 0.10 0.15 0.08],...
    'HorizontalAlignment','left','BackgroundColor',figcolor);

han.happend = uicontrol(han.hpanel2,'Style','pushbutton','Units','normalized','String',...
    'Append new pts','Position',[0.15 0.75 0.7 0.17],...
    'Callback',{@Append_callback,han});

han.hpts = uicontrol(han.hpanel1,'Style','edit','Units','normalized',...
    'Position',[0.50 0.35 0.4 0.1],'BackgroundColor',[1 1 1],...
    'Callback',{@Pts_callback,han});

han.htextsnap = uicontrol(han.hpanel1,'Style','text','Units','normalized','String',...
    'pts.','Position',[0.20 0.33 0.25 0.1],...
    'HorizontalAlignment','right','BackgroundColor',figcolor);

han.hsnap = uicontrol(han.hpanel1,'Style','pushbutton','Units','normalized','String',...
    'Snap to center','Position',[0.15 0.55 0.7 0.17],...
    'Callback',{@Snap_callback,han});

han.hselect = uicontrol(han.hpanel1,'Style','pushbutton','Units','normalized','String',...
    'Select new pts','Position',[0.15 0.75 0.7 0.17],...
    'Callback',{@Select_callback,han});

han.hquit = uicontrol(han.hfig,'Style','pushbutton','Units','normalized','String',...
    'Save & Quit','Position',[0.36 0.03 0.09 0.035],'FontWeight','bold','FontSize',10,...
    'Callback','uiresume(gcbf)');

han.helementtext = uicontrol(han.hpanel1,'Style','text','Units','normalized','String',...
    'element size (pix)','Position',[0.08 0.10 0.4 0.2],...
    'HorizontalAlignment','right','BackgroundColor',figcolor);

han.helement= uicontrol(han.hpanel1,'Style','edit','Units','normalized',...
    'Position',[0.5 0.15 0.2 0.1],'BackgroundColor',[1 1 1],...
    'Callback',{@Element_callback,han});


        %'Units','normalized',

directionstring={'1. Select at least six points.  The first point should be near the middle of';...
                '   the image.  The second point should be one grid spacing over horizontally.';...
                '   The points cannot be colinear, they should generally form a box pattern';...
                '';...
                '2. Snap the points to the middle of the markers.  If the pts. box is empty,';...
                '   then all points will be snapped to center.';...
                '';...                
                '3. Fill in the (x,y,z) for the points.  All z values should be the same.';...
                '   A convenient choice for the (x,y) origin is the first input point in the';...
                '   first image.  The z-origin is the laser sheet.  (x,y,z) must be consistent';...
                '   for all images!';
                '';...                
                '4. Append more points and fill in the corresponding (x,y,z) or use fill xyz.';...
                '   Find XY & xyz may also be used to find points in a grid pattern.';...
                '';...                
                '6. Delete points by deleting an entire table row if needed.  The points';...
                '   should take up the entire target.  Minimum number of points for all planes';...
                '   in a camera is 16 for "cubic xy,linear z," 19 for "cubic xy, quadratic z."';...
                '   More should be selected though.'};

han.hdirections=uicontrol(han.hfig,'Style','text','Units','normalized','String',...
    directionstring,'Position',[0.01 0.01 0.32 0.31],'FontSize',8,...
    'HorizontalAlignment','left','BackgroundColor',figcolor);                
                


han=displaypts(data.X,han);

uiwait(han.hfig);           % this will pause execution until quit button is pressed
%disp('after wait');
guipropsimg = guidata(han.hfig);

[newdata2]=gettabledata(han);
for j=1:size(newdata2,1)
    for k=1:size(newdata2,2)
        if isempty(newdata2{j,k})
            newdata2{j,k}=NaN;
        end
    end
end

newcontrolpts.X=cell2mat(newdata2(:,1:2));

newcontrolpts.x=cell2mat(newdata2(:,3:5));

guidata(han.hfig,guipropsimg);
close(han.hfig);

function Select_callback(hObject,eventdata,han)
    guipropsimg = guidata(han.hfig);
    displaypts([],han);
    figure(han.hfig);
    [c,r,P] = impixel;


    guipropsimg.imgcalpts.X=num2cell([c r]);
    guipropsimg.imgcalpts.x=cell(length(r),3);
    guidata(han.hfig,guipropsimg);
    displaypts([c r],han);
    updatetable(han);

function Append_callback(hObject, eventdata,han)

    [newdata2]=gettabledata(han);
    guipropsimg = guidata(han.hfig);        % check for deleted line
    guipropsimg.imgcalpts.X=newdata2(:,1:2);
    guipropsimg.imgcalpts.x=newdata2(:,3:5);
    %gettabvals(han);
    figure(han.hfig);
    hold on;
    [c,r,P] = impixel;
    guipropsimg.imgcalpts.X=[guipropsimg.imgcalpts.X;num2cell([c r])];

    displaypts(cell2mat(guipropsimg.imgcalpts.X),han);

    guidata(han.hfig,guipropsimg);
    updatetable(han);
    
function Element_callback(hObject,eventdata,han)
    guipropsimg = guidata(han.hfig);
    guipropsimg.elementsize=str2num(get(hObject,'String')); %get element size
    guidata(han.hfig,guipropsimg);
    
function Pts_callback(hObject,eventdata,han)
    guipropsimg = guidata(han.hfig);
    guipropsimg.ptsselected=str2num(get(hObject,'String')); %get element size
    guidata(han.hfig,guipropsimg);
    
function Findxrange_callback(hObject,eventdata,han)
    guipropsimg = guidata(han.hfig);
    guipropsimg.findxrange=str2num(get(hObject,'String')); %get element size
    guidata(han.hfig,guipropsimg);
    
function Findyrange_callback(hObject,eventdata,han)
    guipropsimg = guidata(han.hfig);
    guipropsimg.findyrange=str2num(get(hObject,'String')); %get element size
    guidata(han.hfig,guipropsimg);    
  
function Snap_callback(hObject, eventdata,han)
    guipropsimg = guidata(han.hfig);
    [newdata2]=gettabledata(han);       % always do this at beginning of function in case of user input

    props.d=guipropsimg.markerd;   % marker diameter in mm
    props.ptsselected=guipropsimg.ptsselected;      % %#ok<ST2NM>
    props.elementsize=guipropsimg.elementsize; %get element size
    props.imgmat=guipropsimg.imgmat;
    props.horspace=guipropsimg.horspace;
    props.targcolor=guipropsimg.targcolor;
    
    newdata3=snaptocenter(newdata2,props);
 
    guipropsimg.imgcalpts.X=newdata3(:,1:2);
    guipropsimg.imgcalpts.x=newdata3(:,3:5);
    guidata(han.hfig,guipropsimg);
    updatetable(han);
    displaypts(cell2mat(newdata3(:,1:2)),han);

function Fillxyz_callback(hObject,eventdata,han)
    guipropsimg = guidata(han.hfig);
    [newdata2]=gettabledata(han);       % always do this at beginning of function in case of user input
    [newsize1,newsize2]=size(newdata2);
    ind=[];
    indfind=[];
    for j=1:newsize1
        val3=newdata2{j,3};
        val4=newdata2{j,4};
        validvalue = (~isnan(val3) & ~isnan(val4)) & (~isempty(val3) & ~isempty(val4));
        if validvalue
            ind=[ind j];
        else
            indfind=[indfind j];
        end
    end
    
    fitdatax=cell2mat(newdata2(ind,3:4));
    fitdataX=cell2mat(newdata2(ind,1:2));
    
    finddataX=cell2mat(newdata2(indfind,1:2));
    
    if ~isempty(indfind)
        if length(ind)>=10
            TFORM_Xtot = cp2tform(fitdatax,fitdataX,'polynomial',3);        % make transformation matrix            
        else
            TFORM_Xtot = cp2tform(fitdatax,fitdataX,'polynomial',2);        % make transformation matrix
        end

        [xfind, yfind] = tforminv(TFORM_Xtot,finddataX(:,1),finddataX(:,2)); 
        
        % flag the ones that are 4% away from a real grid point and throw away
        r=find((abs(round(xfind/guipropsimg.horspace)-xfind/guipropsimg.horspace)>0.04) | ...
            ((abs(round(yfind/guipropsimg.verspace)-yfind/guipropsimg.verspace))>0.04));
        
        xfindnorm=round(xfind/guipropsimg.horspace)*guipropsimg.horspace;
        yfindnorm=round(yfind/guipropsimg.verspace)*guipropsimg.verspace;        
        
        indfind(r)=[];      % get rid of the found points so that they wont be saved
        xfindnorm(r)=[];
        yfindnorm(r)=[];
        xfind(r)=[];
        yfind(r)=[];

        if ~isempty(indfind)
            newdata2(indfind,3:4)=num2cell([xfindnorm yfindnorm]);
        end
        
        guipropsimg.imgcalpts.X=newdata2(:,1:2);
        guipropsimg.imgcalpts.x=newdata2(:,3:5);
    
    end
    
    guidata(han.hfig,guipropsimg);
    
    updatetable(han);
 
function FindXYxyz_callback(hObject,eventdata,han)
    guipropsimg = guidata(han.hfig);

    [newdata2]=gettabledata(han);       % always do this at beginning of function in case of user input
    [newsize1,newsize2]=size(newdata2);
    ind=[];
    indfind=[];
    for j=1:newsize1
        val3=newdata2{j,3};
        val4=newdata2{j,4};
        validvalue = (~isnan(val3) & ~isnan(val4)) & (~isempty(val3) & ~isempty(val4));
        if validvalue
            ind=[ind j];
        else
            indfind=[indfind j];
        end
    end
    
    fitdatax=cell2mat(newdata2(ind,3:4));
    fitdataX=cell2mat(newdata2(ind,1:2));
    

    [finddatax finddatay]=meshgrid(guipropsimg.findxrange,guipropsimg.findyrange);
    [rows,cols]=size(finddatax);
    newrows=rows*cols;
    finddatax=reshape(finddatax,rows*cols,1);
    finddatay=reshape(finddatay,rows*cols,1);   

    % find where the new grid points are the same as points already
    % input.  if so, then remove them and keep the old points
    indsame=[];
    for j=1:newrows
        for k=1:size(fitdatax,1)
            if ((finddatax(j)==fitdatax(k,1)) && (finddatay(j)==fitdatax(k,2)))
                indsame=[indsame j];
            end
        end
    end
     
        finddatax(indsame)=[];
        finddatay(indsame)=[];
    
    
    if ~isempty(finddatax)
        if length(ind)>=10
            TFORM_Xtot = cp2tform(fitdataX,fitdatax,'polynomial',3);        % make transformation matrix            
        else
            TFORM_Xtot = cp2tform(fitdataX,fitdatax,'polynomial',2);        % make transformation matrix
        end
        

        [Xfind, Yfind] = tforminv(TFORM_Xtot,finddatax,finddatay); 
                
%         guipropsimg.imgcalpts.X=[guipropsimg.imgcalpts.X;num2cell([Xfind Yfind])];
%         guipropsimg.imgcalpts.x=[guipropsimg.imgcalpts.x;[num2cell([finddatax finddatay]) cell(length(Xfind),1)]];
        guipropsimg.imgcalpts.X=[newdata2(:,1:2);num2cell([Xfind Yfind])];
        guipropsimg.imgcalpts.x=[newdata2(:,3:5);[num2cell([finddatax finddatay]) cell(length(Xfind),1)]];
    
    end

    displaypts(cell2mat(guipropsimg.imgcalpts.X),han);
    guidata(han.hfig,guipropsimg);
    updatetable(han);
    
function han=displaypts(pts,han)

    guipropsimg = guidata(han.hfig);
    
%      names=fieldnames(han);
%      r=max(strcmp(names,'hpanel'));
%       magval=1;
%      if r
%        apiSP = iptgetapi(han.hpanel);
%        magval = apiSP.getMagnification()
%      end

    cla(han.hax);
    axis off;
    han.hax = axes('Parent',han.hfig,'Units','normalized','Position',[0 .333 1 .6666]);
    hold(han.hax,'on');
    han.himage=imshow(guipropsimg.imgmat,'InitialMagnification','fit');
    %[xdat, ydat, A] = getimage(han.himage);
    %set(han.himage,'Xdata',[xdat-0.5],'Ydata',[fliplr(ydat)-0.5]);

    han.hpanel=imscrollpanel(han.hfig,han.himage);
    set(han.hpanel,'Units','normalized','Position',[0 0.3333 1 0.6666]);  % position is [xmin ymin xmax ymax]

    han.hMagBox = immagbox(han.hfig,han.himage);
    

    
    set(han.hMagBox,'Units','normalized','Position',[0 0.965 0.06 0.035]);
    
    %han.hpixinfo=impixelinfo(han.hfig,han.himage);
    %set(han.hpixinfo,'Units','Normalized','Position',[0.75 0.355 0.1 0.02]);   
    
     apiSP = iptgetapi(han.hpanel);
     apiSP.setMagnification(1);



    r=size(pts,1);
    hold(han.hax,'on');
    figure(han.hfig);    
    for j=1:r

         plot(han.hax,pts(j,1),pts(j,2),'pr','MarkerSize',7,'MarkerFaceColor','r');
%         plot(han.hax,pts(j,1),pts(j,2),'pr','MarkerSize',12,'MarkerFaceColor','r');

%         text(pts(j,1)+3,pts(j,2)-2,num2str(j),'FontSize',12,'BackgroundColor',[1 1 1],...
%             'VerticalAlignment','bottom','FontUnits','pixel');
        text(pts(j,1)+10,pts(j,2)-6,num2str(j),'FontSize',7,'BackgroundColor',[1 1 1],...
            'VerticalAlignment','bottom','FontUnits','pixel');
        %text(pts(j,1)+7,pts(j,2)-7,num2str(j),'FontSize',8,'Color','r');
    end
    %end
    guidata(han.hfig,guipropsimg);

function updatetable(han)
    guipropsimg = guidata(han.hfig);


    sztmp=size(guipropsimg.imgcalpts.X,1);
    sztmp2=size(guipropsimg.imgcalpts.x,1);
    tdata=cell(sztmp,5);
    %seeit='off';
    % data=guipropsimg.imgcalpts.X;
    if ~isempty(guipropsimg.imgcalpts.X)
        tdata(:,1:2)=guipropsimg.imgcalpts.X;
    end
    ind=0;
    if ~isempty(guipropsimg.imgcalpts.x)
        if sztmp2<=sztmp
            ind=sztmp2;
        else
            ind=sztmp;
        end
        tdata(1:ind,3:5)=guipropsimg.imgcalpts.x(1:ind,:);
    end
    %tdata;
    [rows,cols]=size(tdata);        
    if rows==1          % this statement fails if only 1 point selected (1 row in tdata)
%         dblArray = javaArray ('java.lang.Double', rows, cols);
%         for m = 1:rows
%            for n = 1:cols
%            dblArray(m,n) = tdata{m,n};
%            end
%         end
%         set(han.htable,'Data',dblArray);        
    else
        set(han.htable,'Data',tdata);
    end
    guidata(han.hfig,guipropsimg);

function [newdata2]=gettabledata(han)
    newdata=get(han.htable,'Data');
    r=size(newdata,1);
    newdata1=cell(r,5);

    for j=1:r       % convert a java object to cell, and convert '' to []
        for k=1:5
            if ischar(newdata(j,k))
                newdata1{j,k}=str2num(newdata(j,k));
            else
                newdata1{j,k}=newdata(j,k);
            end
        end
    end
    flag=[];
    newdata2=newdata1;
    for j=1:r
        logic1=isempty(newdata1{j,1}) || isempty(newdata1{j,2});
        if logic1
            flag=[flag j];
        end
    end

    newdata2(flag,:)=[];        % remove the lines that were deleted

    if isempty(newdata2)
        set(han.htable,'Data',cell(0,5)); 
    elseif size(newdata2,1)==1

        set(han.htable,'Data',newdata2); 
    else
        set(han.htable,'Data',newdata2);
    end