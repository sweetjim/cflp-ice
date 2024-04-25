function [Fkhat] = fftshow(Fk)
Fkhat = fftshift(log10(abs(Fk).^2));
end

