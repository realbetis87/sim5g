function [Structure] = InitializeElementStructure(Structure)
    model=Structure.model;Mesh=model.Mesh;
    Structure.NumberOfElements=size(Mesh.Elements,2);Structure.Elements=Element.empty(0,Structure.NumberOfElements);CE=zeros(model.Geometry.NumCells,Structure.NumberOfElements);
    for ii=1:model.Geometry.NumCells,elements=findElements(model.Mesh,"region","Cell",ii);CE(ii,elements)=1;Structure.Domains(ii).Elements=elements;end
    for kk=1:Structure.NumberOfElements,inDomain=find(CE(:,kk));Structure.Elements(kk)=Element(kk,model.Mesh.Elements(:,kk),inDomain);end
end