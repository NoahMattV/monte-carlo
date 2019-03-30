% Generate the scattering arrays
% To be called at the start of the program

function out = generateScatt()
  global e;
  E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules

  global GTot_G;
  global GTot_X;
  global GTot_L;
  GTot_G = zeros(7, 1001);
  GTot_X = zeros(9, 1001);
  GTot_L = zeros(9, 1001);

  global GTot_G_Max;
  global GTot_X_Max;
  global GTot_L_Max;
  %GTot_G_Max = 0;
  %GTot_X_Max = 0;
  %GTot_L_Max = 0;


  global GG;
  global GX;
  global GL;

  global Gpop_abs; % gamma valley
  global Gpop_em;
  global Gpop_abs_X; % x valley
  global Gpop_em_X;
  global Gpop_abs_L; % l valley
  global Gpop_em_L;

  global Giv_abs_GX;
  global Giv_em_GX;
  global Giv_abs_GL;
  global Giv_em_GL;
  global Giv_abs_XG;
  global Giv_em_XG;
  global Giv_abs_XX;
  global Giv_em_XX;
  global Giv_abs_XL;
  global Giv_em_XL;
  global Giv_abs_LG;
  global Giv_em_LG;
  global Giv_abs_LX;
  global Giv_em_LX;
  global Giv_abs_LL;
  global Giv_em_LL;

  GTot_G(1,:) = GG;
  GTot_G(2,:) = GG + Gpop_abs;
  GTot_G(3,:) = GG + Gpop_abs + Gpop_em;
  GTot_G(4,:) = GG + Gpop_abs + Gpop_em + Giv_abs_GL;
  GTot_G(5,:) = GG + Gpop_abs + Gpop_em + Giv_abs_GL + Giv_em_GL;
  GTot_G(6,:) = GG + Gpop_abs + Gpop_em + Giv_abs_GL + Giv_em_GL + Giv_abs_GX;
  GTot_G(7,:) = GG + Gpop_abs + Gpop_em + Giv_abs_GL + Giv_em_GL + Giv_abs_GX + Giv_em_GX;
  GTot_G_Max = max(GTot_G(7,:));

%{
  GTot_G(1,:) = GG;
  GTot_G(2,:) = Gpop_abs;
  GTot_G(3,:) = Gpop_em;
  GTot_G(4,:) = Giv_abs_GL;
  GTot_G(5,:) = Giv_em_GL;
  GTot_G(6,:) = Giv_abs_GX;
  GTot_G(7,:) = Giv_em_GX;
  GTot_G_Max = max(GTot_G(7,:));
%}

%{
  figure();
  hold on;
  for i = 1:7
    plot(E(1:1001),GTot_G(i,:));
  end
  title('GTot_G old');
  legend('AC','POP_abs','POP_em','GL_abs','GL_em'...
      ,'GX_abs','GX_em');
  hold off;
%}
  
  GTot_X(1,:) = GX;
  GTot_X(2,:) = GX + Gpop_abs_X;
  GTot_X(3,:) = GX + Gpop_abs_X + Gpop_em_X;
  GTot_X(4,:) = GX + Gpop_abs_X + Gpop_em_X + Giv_abs_XG;
  GTot_X(5,:) = GX + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG;
  GTot_X(6,:) = GX + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL;
  GTot_X(7,:) = GX + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL + Giv_em_XL;
  GTot_X(8,:) = GX + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL + Giv_em_XL + Giv_abs_XX;
  GTot_X(9,:) = GX + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL + Giv_em_XL + Giv_abs_XX + Giv_em_XX;
  GTot_X_Max = max(GTot_X(9,:));

%{
  GTot_X(1,:) = GX;
  GTot_X(2,:) = Gpop_abs_X;
  GTot_X(3,:) = Gpop_em_X;
  GTot_X(4,:) = Giv_abs_XG;
  GTot_X(5,:) = Giv_em_XG;
  GTot_X(6,:) = Giv_abs_XL;
  GTot_X(7,:) = Giv_em_XL;
  GTot_X(8,:) = Giv_abs_XX;
  GTot_X(9,:) = Giv_em_XX;
  GTot_X_Max = max(GTot_X(9,:));
%}
%{
  figure();
  hold on;
  for i = 1:9
    plot(E(1:1001),GTot_X(i,:));
  end
  title('GTot_X');
  legend('AC','POP_abs','POP_em','XG_abs','XG_em'...
      ,'XL_abs','XL_em', 'XX_abs', 'XX_em');
  hold off;
%}
  GTot_L(1,:) = GL;
  GTot_L(2,:) = GL + Gpop_abs_L;
  GTot_L(3,:) = GL + Gpop_abs_L + Gpop_em_L;
  GTot_L(4,:) = GL + Gpop_abs_L + Gpop_em_L + Giv_abs_LG;
  GTot_L(5,:) = GL + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG;
  GTot_L(6,:) = GL + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL;
  GTot_L(7,:) = GL + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL + Giv_em_LL;
  GTot_L(8,:) = GL + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL + Giv_em_LL + Giv_abs_LX;
  GTot_L(9,:) = GL + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL + Giv_em_LL + Giv_abs_LX + Giv_em_LX;
  GTot_L_Max = max(GTot_L(9,:));

%{
  GTot_L(1,:) = GL;
  GTot_L(2,:) = Gpop_abs_L;
  GTot_L(3,:) = Gpop_em_L;
  GTot_L(4,:) = Giv_abs_LG;
  GTot_L(5,:) = Giv_em_LG;
  GTot_L(6,:) = Giv_abs_LL;
  GTot_L(7,:) = Giv_em_LL;
  GTot_L(8,:) = Giv_abs_LX;
  GTot_L(9,:) = Giv_em_LX;
  GTot_L_Max = max(GTot_L(9,:));
%}
%{
  figure();
  hold on;
  for i = 1:9
    plot(E(1:1001),GTot_L(i,:));
  end
  title('GTot_L');
  legend('AC','POP_abs','POP_em','LG_abs','LG_em'...
      ,'LL_abs','LL_em', 'LX_abs', 'LX_em');
  hold off;
%}
  out = 1;
end
