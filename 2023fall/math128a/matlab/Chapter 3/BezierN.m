%
%
% draw the letter N and morph it into W 
% with Bezier curves
% 
% written by Ming Gu for Math 128A, Spring 2021
%
x = transpose([3 2 6 5 6.5]);
y = transpose([6 2 6 2 3]);
alphar = transpose([0 2.5 5.0 4.5 6.4]);
betar  = transpose([0 2.5 5.8 2.5 2.8]);
figure(1);
clf;
NP  = 10;
NPP = NP + 20;
scl = 0.75;
for ttt = 1:NPP
  tt = ttt - 1;
  alphar = scl*transpose([0 2.5-tt/NP/3 5.0 4.5 1.5-tt/NP]);
  betar  = scl*transpose([0 2.5+tt/NP/3 5.8 2.5 0.70+tt/NP]);
  x = transpose([3+tt/NP/10 2 6-tt/NP/3 5 5.05+tt/NP]);
  y = transpose([6 2 6-tt/NP/2 2 2.05+tt/NP]);
  alphal = scl*transpose([3.3 2.8 5.8 5.5 0]);
  betal  = scl*transpose([6.5 3.0 5.0 2.2 0]);
  [bx,by]= Bezier(x,y,alphal,betal,alphar,betar);
  bx = sqrt(1/(1+tt/NP))^3* bx;
  by = sqrt(1/(1+tt/NP))^3* by;
  n = length(x)-1;
  NT= 1000;
  t = transpose((0:NT))/NT;
  xp= [];
  yp= [];
  hold on
  axis([0, 6.52, 0, 6.52]);
  for k = 1:n
    xt = bx(k,1) + bx(k,2)*t + bx(k,3) *(t.^2) + bx(k,4) * (t.^2.*(1-t));
    yt = by(k,1) + by(k,2)*t + by(k,3) *(t.^2) + by(k,4) * (t.^2.*(1-t));
    xp = [xp;xt];
    yp = [yp;yt];
    plot(xt,yt,'b.-');
  end
  plot(xp,yp,'b.-');
  hold on
  title(['Morph letter N into W with bezier curves, t = ', num2str(tt)],'FontSize',15);
  if (ttt == 1)
    disp('hit return when ready');
    pause()
  end
  pause(0.05)
  if (ttt < NPP)
    clf;
  end
end
  
  
  
  
