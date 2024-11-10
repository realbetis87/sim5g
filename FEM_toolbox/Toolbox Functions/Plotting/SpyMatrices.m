function [] = SpyMatrices(ToolboxModel)
    if(~isempty(ToolboxModel.Assembled)),assembly=ToolboxModel.Assembled;
        switch assembly.Type
            case "Excitation"
                if(assembly.Dispersive)
                    if(assembly.IsDir)
                        if(assembly.Bian)
                           spy(assembly.TE{1},'b');hold on;spy(assembly.TB{1},'r');spy(assembly.AM{1},'m');spy(assembly.FM{1},'m'); 
                           spy(assembly.TS{1},'g'); spy(assembly.TBC{1},'g');spy(assembly.TG{1},'g');
                           spy(assembly.P{1},'c');spy(assembly.TA{1},'y');spy(assembly.TC{1},'k');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG","P","TA","TC");
                           figure;
                           spy(assembly.TEB{1},'b');hold on;spy(assembly.TBB{1},'r');spy(assembly.AB{1},'g');spy(assembly.FB{1},'m'); 
                           spy(assembly.PB{1},'c');spy(assembly.TAB{1},'y');spy(assembly.TCB{1},'k');
                           title("Matrix B");legend("TE","TB","AM","FM","P","TA","TC");
                        else
                           spy(assembly.TE{1},'b');hold on;spy(assembly.TB{1},'r');spy(assembly.AM{1},'m');spy(assembly.FM{1},'m'); 
                           spy(assembly.TS{1},'g'); spy(assembly.TBC{1},'g');spy(assembly.TG{1},'g');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG");
                           figure;
                           spy(assembly.TEB{1},'b');hold on;spy(assembly.TBB{1},'r');spy(assembly.AB{1},'g');spy(assembly.FB{1},'m'); 
                           title("Matrix B");legend("TE","TB","AM","FM");
                        end
                    else
                        if(assembly.Bian)
                           spy(assembly.TE{1},'b');hold on;spy(assembly.TB{1},'r');spy(assembly.AM{1},'m');spy(assembly.FM{1},'m'); 
                           spy(assembly.TS{1},'g'); spy(assembly.TBC{1},'g');spy(assembly.TG{1},'g');
                           spy(assembly.P{1},'c');spy(assembly.TA{1},'y');spy(assembly.TC{1},'k');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG","P","TA","TC");
                           figure;
                           spy(assembly.TPV{1},'b');
                           title("Matrix B");legend("TPV");
                        else
                           spy(assembly.TE{1},'b');hold on;spy(assembly.TB{1},'r');spy(assembly.AM{1},'m');spy(assembly.FM{1},'m'); 
                           spy(assembly.TS{1},'g'); spy(assembly.TBC{1},'g');spy(assembly.TG{1},'g');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG");
                           figure;
                           spy(assembly.TPV{1},'b');hold on;
                           title("Matrix B");legend("TPV");
                        end
                    end
                else
                    if(assembly.IsDir)
                        if(assembly.Bian)
                           spy(assembly.TE,'b');hold on;spy(assembly.TB,'r');spy(assembly.AM,'m');spy(assembly.FM,'m'); 
                           spy(assembly.TS,'g'); spy(assembly.TBC,'g');spy(assembly.TG,'g');
                           spy(assembly.P,'c');spy(assembly.TA,'y');spy(assembly.TC,'k');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG","P","TA","TC");
                           figure;
                           spy(assembly.TEB,'b');hold on;spy(assembly.TBB,'r');spy(assembly.AB,'g');spy(assembly.FB,'m'); 
                           spy(assembly.PB,'c');spy(assembly.TAB,'y');spy(assembly.TCB,'k');
                           title("Matrix B");legend("TE","TB","AM","FM","P","TA","TC");
                        else
                           spy(assembly.TE,'b');hold on;spy(assembly.TB,'r');spy(assembly.AM,'m');spy(assembly.FM,'m'); 
                           spy(assembly.TS,'g'); spy(assembly.TBC,'g');spy(assembly.TG,'g');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG");
                           figure;
                           spy(assembly.TEB,'b');hold on;spy(assembly.TBB,'r');spy(assembly.AB,'g');spy(assembly.FB,'m'); 
                           title("Matrix B");legend("TE","TB","AM","FM");
                        end
                    else
                        if(assembly.Bian)
                           spy(assembly.TE,'b');hold on;spy(assembly.TB,'r');spy(assembly.AM,'m');spy(assembly.FM,'m'); 
                           spy(assembly.TS,'g'); spy(assembly.TBC,'g');spy(assembly.TG,'g');
                           spy(assembly.P,'c');spy(assembly.TA,'y');spy(assembly.TC,'k');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG","P","TA","TC");
                           figure;
                           spy(assembly.TPV,'b');
                           title("Matrix B");legend("TPV");
                        else
                           spy(assembly.TE,'b');hold on;spy(assembly.TB,'r');spy(assembly.AM,'m');spy(assembly.FM,'m'); 
                           spy(assembly.TS,'g'); spy(assembly.TBC,'g');spy(assembly.TG,'g');
                           title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG");
                           figure;
                           spy(assembly.TPV,'b');hold on;
                           title("Matrix B");legend("TPV");
                        end
                    end
                end
            case "EigenMode"
                 if(assembly.Bian)
                 else
                     if(assembly.Dispersive)
                         figure;spy(assembly.TE{1},'b');hold on;spy(assembly.TB{1},'r');spy(assembly.AM{1},'k');spy(assembly.FM{1},'m'); 
                          spy(assembly.TS{1},'g'); spy(assembly.TBC{1},'g');spy(assembly.TG{1},'g');
                          title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG");
                          figure;spy(assembly.AW{1},'b');hold on;spy(assembly.FW{1},'r');
                          title("Matrix B");legend("AW","FW");
                     else,figure;spy(assembly.TE,'b');hold on;spy(assembly.TB,'r');spy(assembly.AM,'k');spy(assembly.FM,'m'); 
                          spy(assembly.TS,'g'); spy(assembly.TBC,'g');spy(assembly.TG,'g');
                          title("Matrix A");legend("TE","TB","AM","FM","TS","TBC","TG");
                          figure;spy(assembly.AW,'b');hold on;spy(assembly.FW,'r');
                          title("Matrix B");legend("AW","FW");
                     end
                 end
            case "EigenFrequency"
        end
    end
end

