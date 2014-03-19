function [xh]=wlsq(y,H,W)
% --- Weighted Least Squares Fit for Phase Correlation ---
tempmat=sortrows([y',H',W'],2);
y=tempmat(:,1);
H=tempmat(:,2);
W=diag(tempmat(:,3));

% xh=inv(H'*W*H)*H'*W*y;
xh=(H'*W*H)\(H'*W*y);

end
