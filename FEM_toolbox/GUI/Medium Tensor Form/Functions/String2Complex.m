function [num] = String2Complex(str)
    r=split(str,{'+','-'});
    if(isscalar(r)),num=str2double(r{1});
    elseif(numel(r)==2)
        if(isempty(r{1})),num=-str2double(r{2});
        else,num=str2double(r{1});im=r{2};im=im(1:end-1);imag_num=str2double(im);
            if(find(char(str)=='+')),num=num+imag_num*1i;else,num=num-imag_num*1i;end
        end
    elseif(numel(r)==3),num=-str2double(r{2});im=r{3};im=im(1:end-1);imag_num=str2double(im);
        if(find(char(str)=='+')),num=num+imag_num*1i;else,num=num-imag_num*1i;end
    end
end
