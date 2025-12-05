function area = Area(TRI,Vertice_Location)

    Tri_Num = size(TRI,1);
    area = zeros(Tri_Num,1);
    %area2 = zeros(Tri_Num,1);
    for i=1:Tri_Num
        v1 = Vertice_Location(:,TRI(i,1)) - Vertice_Location(:,TRI(i,2));
        v2 = Vertice_Location(:,TRI(i,1)) - Vertice_Location(:,TRI(i,3));
        sinus = sqrt(1-((v1.'*v2)/norm(v1)/norm(v2))^2);
        area(i) = 0.5*norm(v1)*norm(v2)*sinus;
        
        %c = cross(v1,v2);
        %area2(i) = 0.5* norm(c);
    end
    
end