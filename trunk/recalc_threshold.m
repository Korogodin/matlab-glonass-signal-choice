function Ch = recalc_threshold( C, g, h, infC, supC )
%/**
% Функция приведения показателей к общей шкале
%@param C - значение показателя до приведения
%@param g - нижний новый предел
%@param h - верхний новый предел
%@param infC - нижний новый
%@param supC - верхний новый
%*/

Ch = g + (h-g)* (C - infC) ./ (supC - infC);

end

