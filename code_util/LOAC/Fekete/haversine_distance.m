function d = haversine_distance(lat1,lon1,lat2,lon2)

    R = 6371e3;

    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    dlat = lat2 - lat1;
    dlon = lon2 - lon1;

    a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*sin(dlon/2).^2;
    % ---- NUMERICAL SAFETY ----
    a = max(min(a,1),0);
    c = 2*atan2(sqrt(a),sqrt(1-a));

    d = R*c;

end
