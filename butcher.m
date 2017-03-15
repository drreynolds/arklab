function B = butcher(method_name)
% Usage: B = butcher(method_name)
% 
% Returns the butcher table associated with the method given by the input
% argument 'method_name'.  We compute the internal stage times (root nodes)
% c, the Butcher matrix A, the method coefficients b, and the order
% of accuracy for the method q.   The output table has the block structure
%     B = [c, A; q, b]
% for a standard Runge-Kutta method.  
%
% If the method has an embedded error indicator, we also compute the embedded
% method's coefficients b2 and order of accuracy p, and provide the output
%     B = [c, A; q, b; p, b2]
% Methods with this form are indicated below with a (*). 
%
% Method types are specified by the abbreviations:
%       ERK - explicit Runge Kutta (strictly lower-triangular A)
%      DIRK - diagonally-implicit Runge Kutta (lower-triangular A)
%     SDIRK - singly-diagonally implicit Runge Kutta (DIRK with fixed diagonal)
%    ESDIRK - SDIRK method with an initial explicit stage
%
% Allowed methods are listed below.  They are grouped by category
% (ERK, DIRK).  Each method lists the order of accuracy (q)
% and number of stages (s).  Methods with embeddings are marked
% with (*), and their embedding order of accuracy is listed (p):
%
% Explicit:
%    ERK-1-1                            q=1, s=1
%    Heun-Euler-ERK (*)                 q=2, s=2, p=1
%    Ascher(2,3,2)-ERK                  q=2, s=3
%    Ascher(2,2,2)-ERK                  q=2, s=3
%    ARK3(2)4L[2]SA-ERK (*)             q=3, s=4, p=2
%    Bogacki-Shampine-ERK (*)           q=3, s=4, p=2
%    Cash-Karp-ERK (*)                  q=3, s=6, p=3
%    ERK-2-2                            q=3, s=2
%    ERK-3-3                            q=3, s=3
%    Ascher(2,3,3)-ERK                  q=3, s=3
%    Cooper4-ERK                        q=3, s=4
%    Ascher(3,4,3)-ERK                  q=3, s=4
%    Ascher(4,4,3)-ERK                  q=3, s=5
%    Zonneveld-4-3-ERK (*)              q=4, s=5, p=3
%    Fehlberg-ERK (*)                   q=4, s=6, p=3
%    ARK4(3)6L[2]SA-ERK (*)             q=4, s=6, p=3
%    Sayfy-Aburub-4-3-ERK (*)           q=4, s=6, p=3
%    Dormand-Prince-ERK (*)             q=4, s=7, p=3
%    ERK-4-4                            q=4, s=4
%    Merson-5-4-ERK (*)                 q=5, s=5, p=4
%    ARK5(4)8L[2]SA-ERK (*)             q=5, s=8, p=4
%    Cooper6-ERK                        q=5, s=6
%    Verner-6-5-ERK (*)                 q=6, s=8, p=5
%    Fehlberg-8-7-ERK (*)               q=8, s=13, p=7
%
% Diagonally implicit:
%    SDIRK-2-2                          q=2, s=2
%    Ascher(2,3,2)-SDIRK                q=2, s=2
%    Ascher(2,2,2)-SDIRK                q=2, s=2
%    TRBDF2-ESDIRK (*)                  q=3, s=3, p=2
%    Billington-SDIRK (*)               q=3, s=3, p=2
%    Kvaerno(4,2,3)-ESDIRK (*)          q=3, s=4, p=2
%    ARK3(2)4L[2]SA-ESDIRK (*)          q=3, s=4, p=2
%    SDIRK-5-4 (*)                      q=3, s=5, p=3
%    Ascher(2,3,3)-SDIRK                q=3, s=2
%    Ascher(3,4,3)-SDIRK                q=3, s=3
%    Ascher(4,4,3)-SDIRK                q=3, s=4
%    Cooper4-ESDIRK                     q=3, s=4
%    TRX2-ESDIRK (*)                    q=4, s=3, p=2
%    Kvaerno(5,3,4)-ESDIRK (*)          q=4, s=4, p=3
%    Cash(5,3,4)-SDIRK (*)              q=4, s=5, p=3
%    Cash(5,2,4)-SDIRK (*)              q=4, s=5, p=2
%    Sayfy-Aburub-4-3-DIRK (*)          q=4, s=6, p=3
%    ARK4(3)6L[2]SA-ESDIRK (*)          q=4, s=6, p=3
%    Kvaerno(7,4,5)-ESDIRK (*)          q=5, s=7, p=4
%    Ismail(7,4,5)-ESDIRK (*)           q=5, s=7, p=4 (seems broken)
%    ARK5(4)8L[2]SA-ESDIRK (*)          q=5, s=8, p=4
%    Cooper6-ESDIRK                     q=5, s=6
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% set the butcher table
if (strcmp(method_name,'ARK3(2)4L[2]SA-ERK'))

   c = [0; 1767732205903/2027836641118; 3/5; 1];
   b = [1471266399579/7840856788654, -4482444167858/7529755066697, ...
      11266239266428/11593286722821, 1767732205903/4055673282236];
   b2 = [2756255671327/12835298489170, -10771552573575/22201958757719, ...
      9247589265047/10645013368117, 2193209047091/5459859503100];
   A = [0, 0, 0, 0;
      1767732205903/2027836641118, 0, 0, 0;
      5535828885825/10492691773637, 788022342437/10882634858940, 0, 0;
      6485989280629/16251701735622, -4246266847089/9704473918619, ...
	  10755448449292/10357097424841, 0];
   q = 3;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'ARK3(2)4L[2]SA-ESDIRK'))

   c = [0; 1767732205903/2027836641118; 3/5; 1];
   b = [1471266399579/7840856788654, -4482444167858/7529755066697, ...
      11266239266428/11593286722821, 1767732205903/4055673282236];
   b2 = [2756255671327/12835298489170, -10771552573575/22201958757719,  ...
      9247589265047/10645013368117, 2193209047091/5459859503100];
   gamma = 1767732205903/4055673282236;
   A = [0, 0, 0, 0;
      1767732205903/4055673282236, gamma, 0, 0;
      2746238789719/10658868560708, -640167445237/6845629431997, gamma, 0;
      1471266399579/7840856788654, -4482444167858/7529755066697, ...
	  11266239266428/11593286722821, gamma];
   q = 3;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'ARK4(3)6L[2]SA-ERK'))

   c = [0; 1/2; 83/250; 31/50; 17/20; 1];
   b = [82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4];
   b2 = [4586570599/29645900160, 0, 178811875/945068544, ...
      814220225/1159782912, -3700637/11593932, 61727/225920];
   A = [0, 0, 0, 0, 0, 0;
      1/2, 0, 0, 0, 0, 0;
      13861/62500, 6889/62500, 0 0 0 0;
      -116923316275/2393684061468, -2731218467317/15368042101831, ...
	  9408046702089/11113171139209, 0, 0, 0;
      -451086348788/2902428689909, -2682348792572/7519795681897, ...
	  12662868775082/11960479115383, 3355817975965/11060851509271, 0, 0;
      647845179188/3216320057751, 73281519250/8382639484533, ...
	  552539513391/3454668386233, 3354512671639/8306763924573, ...
	  4040/17871, 0];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];
   
elseif (strcmp(method_name,'ARK4(3)6L[2]SA-ESDIRK'))

   c = [0; 1/2; 83/250; 31/50; 17/20; 1];
   b = [82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4];
   b2 = [4586570599/29645900160, 0, 178811875/945068544, ...
      814220225/1159782912, -3700637/11593932, 61727/225920];
   A = [0, 0, 0, 0, 0, 0;
      1/4, 1/4, 0, 0, 0, 0;
      8611/62500, -1743/31250, 1/4, 0, 0, 0;
      5012029/34652500, -654441/2922500, 174375/388108, 1/4, 0, 0;
      15267082809/155376265600, -71443401/120774400, 730878875/902184768, ...
         2285395/8070912, 1/4, 0;
      82889/524892, 0, 15625/83664, 69875/102672, -2260/8211, 1/4];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];
   
elseif (strcmp(method_name,'ARK5(4)8L[2]SA-ERK'))

   c = [0; 41/100; 2935347310677/11292855782101; ...
      1426016391358/7196633302097; 92/100; 24/100; 3/5; 1];
   b = [-872700587467/9133579230613, 0, 0, 22348218063261/9555858737531, ...
	  -1143369518992/8141816002931, -39379526789629/19018526304540, ...
	  32727382324388/42900044865799, 41/200];
   b2 = [-975461918565/9796059967033, 0, 0, 78070527104295/32432590147079, ...
	  -548382580838/3424219808633, -33438840321285/15594753105479, ...
	  3629800801594/4656183773603, 4035322873751/18575991585200];
   A = [0, 0, 0, 0, 0, 0, 0, 0;
        41/100, 0, 0, 0, 0, 0, 0, 0;
        367902744464/2072280473677, 677623207551/8224143866563, 0, 0, 0, 0, 0, 0;
	1268023523408/10340822734521, 0, 1029933939417/13636558850479, 0,0,0,0,0;
	14463281900351/6315353703477, 0, 66114435211212/5879490589093, ...
	    -54053170152839/4284798021562, 0,0,0,0;
	14090043504691/34967701212078, 0, 15191511035443/11219624916014, ...
	    -18461159152457/12425892160975, -281667163811/9011619295870, 0,0,0;
	19230459214898/13134317526959, 0, 21275331358303/2942455364971, ...
	    -38145345988419/4862620318723, -1/8, -1/8, 0, 0;
	-19977161125411/11928030595625, 0, -40795976796054/6384907823539, ...
	    177454434618887/12078138498510, 782672205425/8267701900261, ...
	    -69563011059811/9646580694205, 7356628210526/4942186776405, 0];
   q = 5;
   p = 4;
   B = [c, A; q, b; p, b2];
   
elseif (strcmp(method_name,'ARK5(4)8L[2]SA-ESDIRK'))

   c = [0; 41/100; 2935347310677/11292855782101; ...
      1426016391358/7196633302097; 92/100; 24/100; 3/5; 1];
   b = [-872700587467/9133579230613, 0, 0, 22348218063261/9555858737531, ...
	  -1143369518992/8141816002931, -39379526789629/19018526304540, ...
	  32727382324388/42900044865799, 41/200];
   b2 = [-975461918565/9796059967033, 0, 0, 78070527104295/32432590147079, ...
	  -548382580838/3424219808633, -33438840321285/15594753105479, ...
	  3629800801594/4656183773603, 4035322873751/18575991585200];
   A = [0, 0, 0, 0, 0, 0, 0, 0;
      41/200, 41/200, 0, 0, 0, 0, 0, 0;
      41/400, -567603406766/11931857230679, 41/200, 0, 0, 0, 0, 0;
      683785636431/9252920307686, 0, -110385047103/1367015193373, 41/200, 0,0,0,0;
      3016520224154/10081342136671, 0, 30586259806659/12414158314087, ...
	  -22760509404356/11113319521817, 41/200, 0, 0, 0;
      218866479029/1489978393911, 0, 638256894668/5436446318841, ...
	  -1179710474555/5321154724896, -60928119172/8023461067671, 41/200, 0,0;
      1020004230633/5715676835656, 0, 25762820946817/25263940353407, ...
	  -2161375909145/9755907335909, -211217309593/5846859502534, ...
	  -4269925059573/7827059040749, 41/200, 0;
      -872700587467/9133579230613, 0, 0, 22348218063261/9555858737531, ...
	  -1143369518992/8141816002931, -39379526789629/19018526304540, ...
	  32727382324388/42900044865799, 41/200];
   q = 5;
   p = 4;
   B = [c, A; q, b; p, b2];
   
elseif (strcmp(method_name,'Sayfy-Aburub-4-3-ERK'))
   
   A = [0, 0, 0, 0, 0, 0; ...
        1/2, 0, 0, 0, 0, 0; 
        -1, 2, 0, 0, 0, 0; 
        1/6, 2/3, 1/6, 0, 0, 0;
        0.137, 0.226, 0.137, 0, 0, 0;
        0.452, -0.904, -0.548, 0, 2, 0];
   b = [1/6, 1/3, 1/12, 0, 1/3, 1/12];
   b2 = [1/6, 2/3, 1/6, 0, 0, 0];
   c = [0; 1/2; 1; 1; 1/2; 1];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Sayfy-Aburub-4-3-DIRK'))
   
   A = [0, 0, 0, 0, 0, 0; ...
        1, 0.788675134594813, 0, 0, 0, 0; ...
        2.943375672974064, -2.732050807568877, 0.788675134594813, 0, 0, 0; ...
        1/6, 2/3, 1/6, 0, 0, 0; ...
        -0.423883252702594, 0.20464164, 0.719241612702594, 0, 1.9318, 0; ...
        2.695533010810374, -5.391066021620748, 1.695533010810374, 0, 2, 1.9318];
   b = [1/6, 1/3, 1/12, 0, 1/3, 1/12];
   b2 = [1/6, 2/3, 1/6, 0, 0, 0];
   c = [0; 1/2; 1; 1; 1/2; 1];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Ascher(2,3,3)-ERK'))

   gamma = (3 + sqrt(3))/6;
   c = [0; gamma; 1-gamma];
   b = [0, 1/2, 1/2];
   A = [0, 0, 0;
      gamma, 0, 0;
      gamma-1, 2*(1-gamma), 0];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(2,3,3)-SDIRK'))

   gamma = (3 + sqrt(3))/6;
   c = [gamma; 1-gamma];
   b = [1/2, 1/2];
   A = [gamma, 0;
        1-2*gamma, gamma];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(2,3,2)-ERK'))

   gamma = (2-sqrt(2))/2;
   delta = -2*sqrt(2)/3;
   c = [0; gamma; 1];
   b = [0, 1-gamma, gamma];
   A = [0, 0, 0;
        gamma, 0, 0;
        delta, 1-delta, 0];
   q = 2;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(2,3,2)-SDIRK'))

   gamma = (2-sqrt(2))/2;
   c = [gamma; 1];
   b = [1-gamma, gamma];
   A = [gamma, 0;
        1-gamma, gamma];
   q = 2;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(2,2,2)-ERK'))

   gamma = (2-sqrt(2))/2;
   delta = 1-1/(2*gamma);
   c = [0; gamma; 1];
   b = [delta, 1-delta, 0];
   A = [0, 0, 0;
        gamma, 0, 0;
        delta, 1-delta, 0];
   q = 2;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(2,2,2)-SDIRK'))

   gamma = (2-sqrt(2))/2;
   c = [gamma; 1];
   b = [1-gamma, gamma];
   A = [gamma, 0;
        1-gamma, gamma];
   q = 2;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(3,4,3)-ERK'))

   gamma = 0.4358665215;
   c = [0; gamma; (1+gamma)/2; 1];
   b = [0, -3/2*gamma^2+4*gamma-1/4, 3/2*gamma^2-5*gamma+5/4, gamma];
   A = [0, 0, 0, 0;
        gamma, 0, 0, 0;
	0.3212788860, 0.3966543747, 0, 0;
        -0.105858296, 0.5529291479, 0.5529291479, 0];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(3,4,3)-SDIRK'))

   gamma = 0.4358665215;
   c = [gamma; (1+gamma)/2; 1];
   b = [-3/2*gamma^2+4*gamma-1/4, 3/2*gamma^2-5*gamma+5/4, gamma];
   A = [gamma, 0, 0;
        (1-gamma)/2, gamma, 0;
        -3/2*gamma^2 + 4*gamma - 1/4, 3/2*gamma^2-5*gamma+5/4, gamma];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(4,4,3)-ERK'))

   c = [0; 1/2; 2/3; 1/2; 1];
   b = [1/4, 7/4, 3/4, -7/4, 0];
   A = [0, 0, 0, 0, 0;
        1/2, 0, 0, 0, 0;
	11/18, 1/18, 0, 0, 0;
	5/6, -5/6, 1/2, 0, 0;
	1/4, 7/4, 3/4, -7/4, 0];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Ascher(4,4,3)-SDIRK'))

   c = [1/2; 2/3; 1/2; 1];
   b = [3/2, -3/2, 1/2, 1/2];
   A = [1/2, 0, 0, 0;
        1/6, 1/2, 0, 0;
	-1/2, 1/2, 1/2, 0;
	3/2, -3/2, 1/2, 1/2];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Cooper4-ERK'))

   c = [0; 2/3; 2/3; 1];
   b = [1/4, 1/4, 1/2, 0];
   A = [0, 0, 0, 0;
        2/3, 0, 0, 0;
	1/6, 1/2, 0, 0;
	1/4, 1/4, 1/2, 0];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Cooper4-ESDIRK'))

   c = [0; 2/3; 2/3; 1];
   b = [1/4, 1/4, 1/2, 0];
   A = [0, 0, 0, 0;
        (1-sqrt(3))/6, (3+sqrt(3))/6, 0, 0;
	(5+sqrt(3))/12, -(1+sqrt(3))/4, (3+sqrt(3))/6, 0;
	1/4, 1/4, 1/2, 0];
   q = 3;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Cooper6-ERK'))

   c = [0; 1/2; 1/2; 1/2; 1; 1];
   b = [1/6, 0, 0, 2/3, 1/6, 0];
   A = [0, 0, 0, 0, 0, 0;
        1/2, 0, 0, 0, 0, 0;
	1/4, 1/4, 0, 0, 0, 0;
	1/4, 1/4, 0, 0, 0, 0;
	0, -1, 0, 2, 0, 0;
	1/6, 0, 0, 2/3, 1/6, 0];
   q = 5;
   B = [c, A; q, b];
   
elseif (strcmp(method_name,'Cooper6-ESDIRK'))

   beta = 1.0685790213;
   c = [0; 1/2; 1/2; 1/2; 1; 1];
   b = [1/6, 0, 0, 2/3, 1/6, 0];
   A = [0, 0, 0, 0, 0, 0;
        (1-2*beta)/2, beta, 0, 0, 0, 0;
	1/4, (1-4*beta)/4, beta, 0, 0, 0;
	1/4, beta/2, (1-6*beta)/4, beta, 0, 0;
	0, -2*beta, (1-6*beta-8*beta^2)/(1-4*beta), 4*beta/(1-4*beta), 0, 0;
	1/6, 0, 0, 2/3, 1/6, 0];
   q = 5;
   B = [c, A; q, b];

elseif (strcmp(method_name,'Heun-Euler-ERK'))

   A = [ 0, 0; 1, 0];
   b = [ 0.5, 0.5];
   b2 = [ 1, 0];
   c = [ 0; 1];
   q = 2;
   p = 1;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Bogacki-Shampine-ERK'))

   A = [0, 0, 0, 0; 1/2, 0, 0, 0; 0, 3/4, 0, 0; 2/9, 1/3, 4/9, 0];
   b = [2/9, 1/3, 4/9, 0];
   b2 = [7/24, 1/4, 1/3, 1/8];
   c = [0; 1/2; 3/4; 1];
   q = 3;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Fehlberg-ERK'))

   A = [0, 0, 0, 0, 0, 0; ...
        1/4, 0, 0, 0, 0, 0; ...
	3/32, 9/32, 0, 0, 0, 0; ...
	1932/2197, -7200/2197, 7296/2197, 0, 0, 0; ...
	439/216, -8, 3680/513, -845/4104, 0, 0; ...
	-8/27, 2, -3544/2565, 1859/4104, -11/40, 0];
   b = [ 25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
   b2 = [ 16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55];
   c = [ 0; 1/4; 3/8; 12/13; 1; 1/2];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Cash-Karp-ERK'))

   A = [ 0, 0, 0, 0, 0, 0; ...
         1/5, 0, 0, 0, 0, 0; ...
	 3/40, 9/40, 0, 0, 0, 0; ...
	 3/10, -9/10, 6/5, 0, 0, 0; ...
	 -11/54, 5/2, -70/27, 35/27, 0, 0; ...
	 1631/55296, 175/512, 575/13824, 44275/110592, 253/4096, 0];
   b = [ 37/378, 0, 250/621, 125/594, 0, 512/1771];
   b2 = [ 2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 1/4];
   c = [ 0; 1/5; 3/10; 3/5; 1; 7/8];
   q = 3;
   p = 3;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Dormand-Prince-ERK'))

   A = [ 0, 0, 0, 0, 0, 0, 0; ...
         1/5, 0, 0, 0, 0, 0, 0; ...
	 3/40, 9/40, 0, 0, 0, 0, 0; ...
	 44/45, -56/15, 32/9, 0, 0, 0, 0; ...
	 19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0; ...
	 9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0; ...
	 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
   b = [ 35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
   b2 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
   c = [ 0; 1/5; 3/10; 4/5; 8/9; 1; 1];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'TRBDF2-ESDIRK'))

   A = [0, 0, 0; ...
      (2-sqrt(2))/2, (2-sqrt(2))/2, 0; ...
      sqrt(2)/4, sqrt(2)/4, (2-sqrt(2))/2];
   b = [(1-sqrt(2)/4)/3, (3*sqrt(2)/4+1)/3, (2-sqrt(2))/6];
   b2 = [sqrt(2)/4, sqrt(2)/4, (2-sqrt(2))/2];
   c = [0; 2-sqrt(2); 1];
   q = 3;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'TRX2-ESDIRK'))

   A = [ 0, 0, 0; 0.25, 0.25, 0; 0.25, 0.5, 0.25];
   b = [ 1/6, 2/3, 1/6];
   b2 = [ 0.25, 0.5, 0.25];
   c = [ 0; 0.5; 1];
   q = 4;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Billington-SDIRK'))

   A = [0.292893218813, 0, 0; 0.798989873223, 0.292893218813, 0; ...
        0.740789228841, 0.259210771159, 0.292893218813];
   b = [ 0.691665115992, 0.503597029883, -0.195262145876];
   b2 = [ 0.740789228840, 0.259210771159, 0];
   c = [ 0.292893218813; 1.091883092037; 1.292893218813];
   q = 3;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Cash(5,2,4)-SDIRK'))

   A = [0.435866521508, 0, 0, 0, 0; ...
       -1.13586652150, 0.435866521508, 0, 0, 0; ...
        1.08543330679, -0.721299828287, 0.435866521508, 0, 0; ...
	0.416349501547, 0.190984004184, -0.118643265417, 0.435866521508, 0; ...
	0.896869652944, 0.0182725272734, -0.0845900310706, ...
	    -0.266418670647, 0.435866521508];
   b = [0.896869652944, 0.0182725272734, -0.0845900310706, ...
	  -0.266418670647, 0.435866521508];
   b2 = [(-0.7-0.5)/(-0.7-0.435866521508), ...
	  (0.5-0.435866521508)/(-0.7-0.435866521508), 0, 0, 0];
   c = [0.435866521508; -0.7; 0.8; 0.924556761814; 1];
   q = 4;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Cash(5,3,4)-SDIRK'))

   A = [0.435866521508, 0, 0, 0, 0; ...
       -1.13586652150, 0.435866521508, 0, 0, 0; ...
        1.08543330679, -0.721299828287, 0.435866521508, 0, 0; ...
        0.416349501547, 0.190984004184, -0.118643265417, 0.435866521508, 0;...
        0.896869652944, 0.0182725272734, -0.0845900310706, ...
	  -0.266418670647, 0.435866521508];
   b = [0.896869652944, 0.0182725272734, -0.0845900310706, ...
	  -0.266418670647, 0.435866521508];
   b2 = [0.776691932910, 0.0297472791484, -0.0267440239074, 0.220304811849, 0];
   c = [0.435866521508; -0.7; 0.8; 0.924556761814; 1];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Kvaerno(4,2,3)-ESDIRK'))

   A = [0, 0, 0, 0; ...
        0.4358665215, 0.4358665215, 0, 0; ...                      
        0.490563388419108, 0.073570090080892, 0.4358665215, 0; ...
        0.308809969973036, 1.490563388254106, -1.235239879727145, 0.4358665215];
   b = [0.308809969973036, 1.490563388254106, -1.235239879727145, 0.4358665215];
   b2 = [0.490563388419108, 0.073570090080892, 0.4358665215, 0];
   c = [0; 0.871733043; 1; 1];
   q = 3;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Kvaerno(5,3,4)-ESDIRK'))

   A = [0, 0, 0, 0, 0; ...
        0.4358665215, 0.4358665215, 0, 0, 0; ...
        0.140737774731968, -0.108365551378832, 0.4358665215, 0, 0;...
        0.102399400616089, -0.376878452267324, 0.838612530151233, 0.4358665215, 0;...
        0.157024897860995, 0.117330441357768, 0.61667803039168, -0.326899891110444, 0.4358665215];
   b = [0.157024897860995, 0.117330441357768, 0.61667803039168, -0.326899891110444, 0.4358665215];
   b2 = [0.102399400616089, -0.376878452267324, 0.838612530151233, 0.4358665215, 0];
   c = [0; 0.871733043; 0.468238744853136; 1; 1];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Kvaerno(7,4,5)-ESDIRK'))
   
   A = [0, 0, 0, 0, 0, 0, 0; ...
        0.26, 0.26, 0, 0, 0, 0, 0; ...
        0.13, 0.84033320996790809, 0.26, 0, 0, 0, 0; ...
        0.22371961478320505, 0.47675532319799699, -0.06470895363112615, 0.26, 0, 0, 0; ...
        0.16648564323248321, 0.10450018841591720, 0.03631482272098715, -0.13090704451073998, 0.26, 0, 0; ...
        0.13855640231268224, 0, -0.04245337201752043, 0.02446657898003141, 0.61943039072480676, 0.26, 0; ...
        0.13659751177640291, 0, -0.05496908796538376, -0.04118626728321046, 0.62993304899016403, 0.06962479448202728, 0.26];

   b = [0.13659751177640291, 0, -0.05496908796538376, -0.04118626728321046, 0.62993304899016403, 0.06962479448202728, 0.26];
   b2 = [0.13855640231268224, 0, -0.04245337201752043, 0.02446657898003141, 0.61943039072480676, 0.26, 0];
   c = [0; 0.52; 1.230333209967908; 0.895765984350076; 0.436393609858648; 1; 1];
   q = 5;
   p = 4;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Ismail(7,4,5)-ESDIRK'))

   A = [0, 0, 0, 0, 0, 0, 0; ...
        0.28589, 0.28589, 0, 0, 0, 0, 0; ...
        0.142945000375866, 0.924011005, 0.28589, 0, 0, 0, 0; ...
        0.168035986, -0.04941651, -0.004509476, 0.28589, 0, 0, 0; ...
        0.182315003, -0.112951603, -0.027793233, 0.422539833, 0.28589, 0, 0; ...
        0.247563917, -0.425378071, -0.107036282, 0.395700134, 0.503260302, 0.28589, 0; ...
        0.130018035, 0, -0.019290177, 0.535386266, 0.234313169, -0.166317293, 0.28589];
   b = [0.130018035, 0, -0.019290177, 0.535386266, 0.234313169, -0.166317293, 0.28589];
   c = [0; 0.57178; 1.352846005375866; 0.4; 0.75; 0.9; 1];

   % embedding is broken (though this is exactly what's in their paper)
   b2 = [-0.094388662, 0, -0.039782614, 0.745608552, -0.505129807, 0.704915206, 0.28589];
   q = 5;
   p = 4;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'ERK-1-1'))
   
   A = [0];
   b = [1];
   c = [0];
   q = 1;
   B = [c, A; q, b];

elseif (strcmp(method_name,'ERK-2-2'))
   
   A = [ 0, 0; 2/3, 0];
   b = [ 1/4, 3/4];
   c = [ 0; 2/3];
   q = 3;
   B = [c, A; q, b];

elseif (strcmp(method_name,'ERK-3-3'))
   
   A = [ 0, 0, 0; 1/2, 0, 0; -1, 2, 0];
   b = [ 1/6, 2/3, 1/6];
   b2 = [0, 1, 0];
   c = [ 0; 1/2; 1];
   q = 4;
   p = 2;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'ERK-4-4'))
   
   A = [0, 0, 0, 0; 1/2, 0, 0, 0; 0, 1/2, 0, 0; 0, 0, 1, 0];
   b = [1/6, 1/3, 1/3, 1/6];
   c = [0; 1/2; 1/2; 1];
   q = 4;
   B = [c, A; q, b];

elseif (strcmp(method_name,'Merson-5-4-ERK'))
   
   A = [0, 0, 0, 0, 0; 1/3, 0, 0, 0, 0; 1/6, 1/6, 0, 0, 0; ...
        1/8, 0, 3/8, 0, 0; 1/2, 0, -3/2, 2, 0];
   b = [1/6, 0, 0, 2/3, 1/6];
   b2 = [1/10, 0, 3/10, 2/5, 1/5];
   c = [0; 1/3; 1/3; 1/2; 1];
   q = 5;
   p = 4;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Zonneveld-4-3-ERK'))
   
   A = [0, 0, 0, 0, 0; 1/2, 0, 0, 0, 0; 0, 1/2, 0, 0, 0; ...
        0, 0, 1, 0, 0; 5/32, 7/32, 13/32, -1/32, 0];
   b = [1/6, 1/3, 1/3, 1/6, 0];
   b2 = [-1/2, 7/3, 7/3, 13/6, -16/3];
   c = [0; 1/2; 1/2; 1; 3/4];
   q = 4;
   p = 3;
   B = [c, A; q, b; p, b2];


elseif (strcmp(method_name,'Verner-6-5-ERK'))

   A = [0, 0, 0, 0, 0, 0, 0, 0; ...
        1/6, 0, 0, 0, 0, 0, 0, 0; ...
        4/75, 16/75, 0, 0, 0, 0, 0, 0; ...
        5/6, -8/3, 5/2, 0, 0, 0, 0, 0; ...
        -165/64, 55/6, -425/64, 85/96, 0, 0, 0, 0; ...
        12/5, -8, 4015/612, -11/36, 88/255, 0, 0, 0; ...
        -8263/15000, 124/75, -643/680, -81/250, 2484/10625, 0, 0, 0; ...
        3501/1720, -300/43, 297275/52632, -319/2322, 24068/84065, 0, 3850/26703, 0];
   b = [3/40, 0, 875/2244, 23/72, 264/1955, 0, 125/11592, 43/616];
   b2 = [13/160, 0, 2375/5984, 5/16, 12/85, 3/44, 0, 0];

   c = [0; 1/6; 4/15; 2/3; 5/6; 1; 1/15; 1];
   q = 6;
   p = 5;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'Fehlberg-8-7-ERK'))

   A = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
        1/20, 0, 0, 1/4, 1/5, 0, 0, 0, 0, 0, 0, 0, 0; ...
        -25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0, 0; ...
        31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0, 0; ...
        2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0, 0; ...
        -91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0, 0; ...
        2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0, 0; ...
        3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0, 0; ...
        -1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1, 0];
   b = [0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840];
   b2 = [41/840, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 41/840, 0, 0];

   c = [0; 2/27; 1/9; 1/6; 5/12; 1/2; 5/6; 1/6; 2/3; 1/3; 1; 0; 1];
   q = 8;
   p = 7;
   B = [c, A; q, b; p, b2];

elseif (strcmp(method_name,'SDIRK-2-2'))
   
   A = [1-1/sqrt(2), 0; 1/sqrt(2), 1-1/sqrt(2)];
   b = [ 1/sqrt(2), 1 - 1/sqrt(2)];
   c = [ 1-1/sqrt(2); 1];
   q = 2;
   B = [c, A; q, b];

else
   
   B = 0;
   fprintf('Butcher error, method %s not defined\n',method_name);
   
end

% end of function
