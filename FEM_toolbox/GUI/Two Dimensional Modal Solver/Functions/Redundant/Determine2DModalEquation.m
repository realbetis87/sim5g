%--------------------------------------------------------------------------
%{
            Formulation selection for the TwoDModalSolver app

    1. if any of the app.BoundaryIndices - Boundaries contains
    Bianisotropic Media -> Equation=3 (Vector Wave Equation E)
    2. if any of the app.BoundaryIndices - Boundaries contains 
    Non Orthotropic Media -> Equation=2 (TFTF EH)
    3. else Equation=1 (TFLF EH)
%}
%--------------------------------------------------------------------------
function [] = Determine2DModalEquation(app),bianFlag=false;orthFlag=false;
    for ii=1:numel(app.BoundaryIndices),boundary=app.TModel.Boundaries(app.BoundaryIndices(ii));element=app.TModel.Facets(boundary.Facets(1));medium=element.Medium2D;
        if(medium.Type~="Iso"),if(medium.Type=="Bian"),bianFlag=true;
                               elseif(CheckNonOrthotropicMedium(boundary,medium)),orthFlag=true;
                               end
        end
    end
    if(bianFlag==true),app.Equation=3;
    elseif(bianFlag==false && orthFlag==true),app.Equation=2;
    else,app.Equation=1;
    end
end

function [res] = CheckNonOrthotropicMedium(boundary,medium),res=false;
    if(medium.IsDispersive)
        for ii=1:medium.FRange.NF,r1=CheckNonOrthotropicTensor(boundary,medium.Epsilon{ii});r2=CheckNonOrthotropicTensor(boundary,medium.Mu{ii});
            if(r1 || r2),res=true;end
        end
    else,r1=CheckNonOrthotropicTensor(boundary,medium.Epsilon);r2=CheckNonOrthotropicTensor(boundary,medium.Mu);if(r1 || r2),res=true;end
    end
end

function [res] = CheckNonOrthotropicTensor(boundary,tensor),res=false;
    switch boundary.Axis
        case  1,if(tensor(1,2)~=0 || tensor(1,3)~=0 || tensor(3,1)~=0 || tensor(2,1)~=0),res=true;end
        case -1,if(tensor(1,2)~=0 || tensor(1,3)~=0 || tensor(3,1)~=0 || tensor(2,1)~=0),res=true;end
        case  2,if(tensor(2,1)~=0 || tensor(2,3)~=0 || tensor(1,2)~=0 || tensor(3,2)~=0),res=true;end
        case -2,if(tensor(2,1)~=0 || tensor(2,3)~=0 || tensor(1,2)~=0 || tensor(3,2)~=0),res=true;end
        case  3,if(tensor(3,1)~=0 || tensor(3,2)~=0 || tensor(2,3)~=0 || tensor(1,3)~=0),res=true;end
        case -3,if(tensor(3,1)~=0 || tensor(3,2)~=0 || tensor(2,3)~=0 || tensor(1,3)~=0),res=true;end
    end
end