function [] = PlotEdges(varargin),global Edges;global Vertices;
if (nargin==1),Colors=['k','b','r','g','y','m'];cc=0;id=varargin{1};
    for ii=1:numel(id)
        for jj=1:numel(Edges),edge=Edges(jj);if(edge.Id==id(ii)),nodes=edge.Vertices;vertices=Vertices(nodes);plot3([vertices(1).X vertices(2).X],[vertices(1).Y vertices(2).Y],[vertices(1).Z vertices(2).Z],'Color',Colors(ii));hold on;cc=cc+1;end,end
    end
    disp(int2str(cc));
elseif(nargin==2),LocalEdges=varargin{1};id=varargin{2};
    for ii=1:numel(id)
        for jj=1:numel(LocalEdges),edge=LocalEdges(jj);if(edge.Id==id(ii)),nodes=edge.Vertices;vertices=Vertices(nodes);plot3([vertices(1).X vertices(2).X],[vertices(1).Y vertices(2).Y],[vertices(1).Z vertices(2).Z],'Color',Colors(ii));hold on;cc=cc+1;end,end
    end
    disp(int2str(cc));
end
end

