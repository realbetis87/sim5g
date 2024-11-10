function [] = SingleFrequency_GUI(app)
        app.FrequencyGrid.RowHeight={'1x','1x','1x','2x','1x','0x','1x'};
        app.SingleFrequencyPanel.Enable=true;app.SingleFrequencyPanel.Visible=true;
        app.MFPanel.Enable=false;app.MFPanel.Visible=false;
        app.StartStopPanel.Visible=false;app.StartStopPanel.Enable=false;
        app.VectorPanel.Visible=false;app.VectorPanel.Enable=false;app.FrequencyMode=1;
end

