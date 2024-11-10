function [] = GUIActive(app,Message),app.Occupied=~app.Occupied;
    if(~app.Occupied),app.Lamp.Color=[0.47,0.67,0.19];app.ConditionLabel.Text=Message;app.ConditionLabel.FontColor=[0.47,0.67,0.19];
    else,app.Lamp.Color=[0.85,0.33,0.10];app.ConditionLabel.Text=Message;app.ConditionLabel.FontColor=[0.85,0.33,0.10];
    end
end

