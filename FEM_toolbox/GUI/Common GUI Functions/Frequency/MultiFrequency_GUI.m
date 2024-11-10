function [] = MultiFrequency_GUI(app)
     app.SingleFrequencyPanel.Enable=false;app.SingleFrequencyPanel.Visible=false;
     app.MFPanel.Enable=true;app.MFPanel.Visible=true;val=get(app.StartStepStopButton,'value');
     if(val==1),app.StartStopPanel.Visible=true;app.StartStopPanel.Enable=true;app.FrequencyMode=2;
                app.VectorPanel.Visible=false;app.VectorPanel.Enable=false;
                app.FrequencyGrid.RowHeight={'1x','0x','1x','2x','1x','1x','1x'};
     else,app.VectorPanel.Visible=true;app.VectorPanel.Enable=true;app.FrequencyMode=3;
            app.VectorPanel.Visible=true;app.VectorPanel.Enable=true;
            app.FrequencyGrid.RowHeight={'1x','0x','1x','0x','1x','3x','1x'};
     end
end

