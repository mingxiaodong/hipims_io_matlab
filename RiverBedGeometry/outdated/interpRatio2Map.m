function hq=interpRatio2Map(lwh_line1,lwh_line2, lq,wq)
hq = nan(size(lq));
for i=1:length(lq)
    l_value = lq(i);
    w_value = wq(i);
    % 1D interp h based on w coords
    [~,ia,~]=unique(lwh_line1(:,2));
    h1 = interp1(lwh_line1(ia,2),lwh_line1(ia,3),w_value,'PCHIP');
    [~,ia,~]=unique(lwh_line2(:,2));
    h2 = interp1(lwh_line2(ia,2),lwh_line2(ia,3),w_value,'PCHIP');
    hq(i) = h1+(h2-h1)*l_value;
    
end
end
    