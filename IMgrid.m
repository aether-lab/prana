function [X,Y]=IMgrid(L,S,G)
% --- Grid Generation Subfunction ---

%grid buffer
if nargin<3
    G=[0 0 0 0];
end

S=[S(2) S(1)];
G=[G(2) G(1) L(1)-G(2)+1 L(2)-G(1)+1];

%form grid
if max(S)==0
    %pixel grid
    y=(1:L(1))';
    x=1:L(2);
else
    if G(1)==0
        %buffers 1/2 grid spacing
        y=(ceil((L(1)-(floor(L(1)/S(1))-2)*S(1))/2):S(1):(L(1)-S(1)))';
    else
        %predefined grid buffer
        y=(G(1):S(1):G(3))';
    end
    if G(2)==0
        %buffers 1/2 grid spacing
        x=ceil((L(2)-(floor(L(2)/S(2))-2)*S(2))/2):S(2):(L(2)-S(2));
    else
        %predefined grid buffer
        x=(G(2):S(2):G(4));
    end
end

%vector2matrix conversion
X=x(ones(length(y),1),:);
Y=y(:,ones(1,length(x)));