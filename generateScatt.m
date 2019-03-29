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
  GTot_G_Max = 0;
  GTot_X_Max = 0;
  GTot_L_Max = 0;


  global G;
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

  GTot_G(1,:) = G;
  GTot_G(2,:) = G + Gpop_abs;
  GTot_G(3,:) = G + Gpop_abs + Gpop_em;
  GTot_G(4,:) = G + Gpop_abs + Gpop_em + Giv_abs_GL;
  GTot_G(5,:) = G + Gpop_abs + Gpop_em + Giv_abs_GL + Giv_em_GL;
  GTot_G(6,:) = G + Gpop_abs + Gpop_em + Giv_abs_GL + Giv_em_GL + Giv_abs_GX;
  GTot_G(7,:) = G + Gpop_abs + Gpop_em + Giv_abs_GL + Giv_em_GL + Giv_abs_GX + Giv_em_GX;
  GTot_G_Max = max(GTot_G(7,:));

  %{
  figure();
  hold on;
  for i = 1:7
    plot(E(1:1001),GTot_G(i,:));
  end
  title('GTot_G');
  legend('AC','AC+POP_abs','AC+POP_abs+POP_em','AC+POP+GL_abs','AC+POP+GL_abs+GL_em'...
      ,'AC+POP+GL+GX_abs','AC+POP+GL+GX_abs+GL_em');
  %}

  GTot_X(1,:) = G;
  GTot_X(2,:) = G + Gpop_abs_X;
  GTot_X(3,:) = G + Gpop_abs_X + Gpop_em_X;
  GTot_X(4,:) = G + Gpop_abs_X + Gpop_em_X + Giv_abs_XG;
  GTot_X(5,:) = G + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG;
  GTot_X(6,:) = G + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL;
  GTot_X(7,:) = G + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL + Giv_em_XL;
  GTot_X(8,:) = G + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL + Giv_em_XL + Giv_abs_XX;
  GTot_X(9,:) = G + Gpop_abs_X + Gpop_em_X + Giv_abs_XG + Giv_em_XG + Giv_abs_XL + Giv_em_XL + Giv_abs_XX + Giv_em_XX;
  GTot_X_Max = max(GTot_X(9,:));

  %{
  figure();
  hold on;
  for i = 1:9
    plot(E(1:1001),GTot_X(i,:));
  end
  title('GTot_X');
  %}

  GTot_L(1,:) = G;
  GTot_L(2,:) = G + Gpop_abs_L;
  GTot_L(3,:) = G + Gpop_abs_L + Gpop_em_L;
  GTot_L(4,:) = G + Gpop_abs_L + Gpop_em_L + Giv_abs_LG;
  GTot_L(5,:) = G + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG;
  GTot_L(6,:) = G + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL;
  GTot_L(7,:) = G + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL + Giv_em_LL;
  GTot_L(8,:) = G + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL + Giv_em_LL + Giv_abs_LX;
  GTot_L(9,:) = G + Gpop_abs_L + Gpop_em_L + Giv_abs_LG + Giv_em_LG + Giv_abs_LL + Giv_em_LL + Giv_abs_LX + Giv_em_LX;
  GTot_L_Max = max(GTot_L(9,:));

  %{
  figure();
  hold on;
  for i = 1:9
    plot(E(1:1001),GTot_L(i,:));
  end
  title('GTot_L');
  %}

  out = 1;
end
