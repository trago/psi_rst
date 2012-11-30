function [alpha] = AntiAliasing(alpha)
[m,k] = size(alpha);
for n=2:k
    if alpha(n)-alpha(n-1) > pi
        alpha(n) = alpha(n)-2*pi;
    end
    if alpha(n)-alpha(n-1) < -pi
        alpha(n) = alpha(n)+2*pi;
    end
end