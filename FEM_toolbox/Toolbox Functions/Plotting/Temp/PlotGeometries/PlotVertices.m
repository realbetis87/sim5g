function [] = PlotVertices(id),global Vertices;Colors=['k','b','r','g','y','m'];cc=0;
    for ii=1:numel(id)
        for jj=1:numel(Vertices),vertex=Vertices(jj);if(vertex.Id==id(ii)),plot3([vertex.X],[vertex.Y],[vertex.Z],'*','Color',Colors(ii));hold on;cc=cc+1;end,end
    end
    disp(int2str(cc));
end

