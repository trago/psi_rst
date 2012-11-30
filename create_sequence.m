function [I dc]=create_sequence(N, phase, alpha, alpha_noise, frame_noise, dc_pow)

[m n] = size(phase);

dc = makeParabola(m,n,dc_pow);

makeI(dc,1,phase,alpha, alpha_noise,N,frame_noise);