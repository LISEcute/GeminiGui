#include "CYrast.h"


#include <QTextStream>
#include <QFile>
#include <cmath>


#include <QDebug>

CYrast* CYrast::fInstance = 0;  // singleton

double CYrast::addBar = 0.;
double const CYrast::pi=acos(-1.);
double const CYrast::deltaJ= 3.;
bool CYrast::first = 1;
bool CYrast::bForceSierk = 0;
double const CYrast::kRotate=41.563;

//RLDM constants
double const CYrast::x1h[11][6]={
{0.0,0.0,0.0,0.0,0.0,0.0},
{-0.0057,-0.0058,-0.006,-0.0061,-0.0062,-0.0063},
{-0.0193,-0.0203,-0.0211,-0.022,-0.023,-0.0245},
{-0.0402, -0.0427,-0.0456,-0.0497,-0.054,-0.0616},
{-0.0755,-0.0812,-0.0899,-0.0988,-0.109,-0.12},
{-0.1273,-0.1356,-0.147,-0.1592,-0.1745,-0.1897},
{-0.1755,-0.1986,-0.2128,-0.2296,-0.251,-0.26},
{-0.255,-0.271,-0.291,-0.301,-0.327,-0.335},
{-0.354,-0.36,-0.365,-0.372, -0.403,-0.42},
{0.0,0.0,0.0,0.0,0.0,0.0},
{0.0,0.0,0.0,0.0,0.0,0.0}};

double const CYrast::x2h[11][6]={
{0.0,0.0,0.0,0.0,0.0,0.0},
{-0.0018,-0.0019,-0.00215,-0.0024,-0.0025,-0.003},
{-0.0063,-0.00705,-0.0076,-0.0083,-0.0091,-0.0095},
{-0.015,-0.0158,-0.0166,-0.0192,-0.0217,-0.025},
{-0.0245,-0.0254,-0.029,-0.0351,-0.0478,-0.0613},
{-0.0387,-0.0438,-0.0532,-0.0622,-0.0845,-0.0962},
{-0.0616,-0.0717,-0.0821,-0.0972,-0.1123,-0.1274},
{-0.0793, -0.1014,-0.1138,-0.1262,-0.1394,-0.1526},
{-0.12,-0.134,-0.1503,-0.1666,-0.1829,-0.1992},
{-0.1528,-0.171,-0.1907,-0.2104,-0.2301,-0.2498}
,{0.0,0.0,0.0,0.0,0.0,0.0}};

double const CYrast::x3h[20][10]={
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
{-0.00012,-0.00014,-0.00016,-0.00018,-0.0002,-0.00024,-0.00029,-0.00036,
-0.00065,-0.00089},
{-0.00047,-0.0005,-.00058,-.00065,-.00074,-.00085,-.00101,-.00124,
-.00138,-.00178},
{-0.001,-0.00105,-0.00124,-0.00138,-0.00156,-0.00179,-0.00275,-0.00292,
-0.003,-0.003},
{-0.00176,-0.0019,-0.00211,-0.00235,-0.00263,-0.00298,-0.00449,-0.0053,
-0.0053,-0.0053},
{-0.003,-0.00308,-0.00318,-0.00352,-0.00392,-0.00417,-0.0062,-0.0062,
-0.0062,-0.0062},
{-0.00374,-0.0041,-0.00444,-0.00488,-0.00521,-0.00545,-0.0066,-0.0066,
-0.0066,-0.0066},
{-0.0053,-0.0055,-0.00585,-0.0064,-0.00695,-0.007,-0.007,-0.007,-0.007,-0.007},
{-0.00632,-0.007,-0.00742,-0.00792,-0.00856,-0.009,-0.009,-0.009,
-0.009,-0.009},
{-0.0079,-0.0085,-0.01022,-0.0119,-0.012,-0.012,-0.012,-0.012,-0.012,-0.012},
{-0.00944,-0.0102,-0.0142,-0.0182,-0.019,-0.019,-0.019,-0.019,-0.019,-0.019},
{-0.0112,-0.0133,-0.0182,-0.0238,-0.024,-0.024,-0.024,-0.024,-0.024,-0.024},
{-0.01303,-0.0178,-0.0226,-0.0274,-0.028,-0.028,-0.028,-0.028,-0.028,-0.028},
{-0.0165,-0.0254,-0.0343,-0.0343,-0.034,-0.034,-0.034,-0.034,-0.034,-0.034},
{-0.0203,-0.033,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04,-0.04},
{-0.025,-0.0406,-0.046,-0.047,-0.047,-0.047,-0.047,-0.047,-0.047,-0.047},
{-0.03036,-0.0482,-0.048,-0.048,-0.048,-0.048,-0.048,-0.048,-0.048,-0.048},
{-0.0363,-0.0558,-0.056,-0.056,-0.056,-0.056,-0.056,-0.056,-0.056,-0.056},
{-0.04234,-0.0634,-0.064,-0.064,-0.064,-0.064,-0.064,-0.064-0.064,-0.064},
{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};


double const CYrast::x1b[11][6]={
{0.28,.243,.221,.208,.195,.18},
{.211,.186,.17,.1506,.136,.12},
{0.152,.131,.1155,.096,.0795,.0625},
{.09725,.0795,.065,.0506,.0375,0.0253},
{.05771,.0455,.03414,.0235,.014,.0065},
{.03325,.0235,.0153,.0081,0.001,.0},
{.01625,.009,.0032,.0,.0,.0},
{.0071,.0,.0,.0,.0,.0},
{.0,.0,0.0,.0,.0,.0},
{.0,.0,.0,.0,.0,.0},
{0.,0.,0.,0.,0.,0.}};

double const CYrast::x2b[11][6]={
{.18,.1695,.1515,.133,.1155,.0949},
{.1495,.1363,.1165,.099,.0815,.0594},
{.12,.1032,.0864,.0678,.0469,.028},
{.09,.0725,.0556,.037,.019,.0057},
{.0625,.045,.0304,.016,.005,0.},
{.0406,.0264,.0151,.0052,0.,0.},
{.0253,.0144,.0027,0.,0.,0.},
{.0141,.006,0.,0.,0.,0.},
{.0065,.0008, 0.,0.,0.,0.},
{.002,0.,0.,0.,0.,0.},
{0.,0.,0.,0.,0.,0.}};


double const CYrast::x3b[20][10]={
{.0949,.0755,.0564,.0382,.0223,.0121,.00588,.00242,.00069,.0001},
{.0873,.0684,.049,.0306,.0162,.0074,.00267,.00055,0.,0.},
{.0801,.061,.0418,.0235,.0108,.00373,.00071,0.,0.,0.},
{.073,.054,.035,.0178,.0062,.00125,0.,0.,0.,0.},
{.0661,.047,.0284,.012,.0025,0.,0.,0., 0.,0.},
{.0594,.0404,.022,.0065,0.,0.,0.,0.,0.,0.},
{.0528,.034,.0159,.002,0.,0.,0.,0.,0.,0.},
{.0465,.0277,.01,0.,0.,0.,0.,0.,0.,0.},
{.0401,.0217,.0044,0.,0.,0.,0.,0.,0.,0.},
{.0339,.0158,.00024,0.,0.,0.,0., 0.,0.,0.},
{.028,.0106,0.,0.,0.,0.,0.,0.,0.,0.},
{.0219,.0064,0.,0.,0., 0.,0.,0.,0.,0.},
{.0164,.0025,0.,0.,0.,0.,0.,0.,0.,0.},
{.0122,0.,0.,0., 0.,0.,0.,0.,0.,0.},
{.0085,0.,0.,0.,0.,0.,0.,0.,0.,0.},
{.0057,0.,0., 0.,0.,0.,0.,0.,0.,0.},
{.0035,0.,0.,0.,0.,0.,0.,0.,0.,0.},
{.0016,0.,0.,0.,0.,0.,0.,0.,0.,0.},
{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
{0.,0.,0., 0.,0.,0.,0.,0.,0.,0.}};


//sierk constants
double const CYrast:: emncof[4][5]={
{-9.01100e2,-1.40818e3, 2.77000e3,-7.06695e2, 8.89867e2},
{1.35355e4,-2.03847e4, 1.09384e4,-4.86297e3,-6.18603e2},
{-3.26367e3, 1.62447e3, 1.36856e3, 1.31731e3, 1.53372e2},
{7.48863e3,-1.21581e4, 5.50281e3,-1.33630e3, 5.05367e-2}};

double const CYrast:: elmcof[4][5]={
{ 1.84542e3,-5.64002e3, 5.66730e3,-3.15150e3, 9.54160e2},
{-2.24577e3, 8.56133e3,-9.67348e3, 5.81744e3,-1.86997e3},
{ 2.79772e3,-8.73073e3, 9.19706e3,-4.91900e3, 1.37283e3},
{-3.01866e1, 1.41161e3,-2.85919e3, 2.13016e3,-6.49072e2}};

double const CYrast::emxcof[5][7]={
{-4.10652732e6, 1.00064947e7,-1.09533751e7, 7.84797252e6,-3.78574926e6,
 1.12237945e6,-1.77561170e5},
{ 1.08763330e7,-2.63758245e7, 2.85472400e7,-2.01107467e7, 9.48373641e6,
-2.73438528e6, 4.13247256e5},
{-8.76530903e6, 2.14250513e7,-2.35799595e7, 1.70161347e7,-8.23738190e6,
 2.42447957e6,-3.65427239e5},
{ 6.30258954e6,-1.52999004e7, 1.65640200e7,-1.16695776e7, 5.47369153e6,
-1.54986342e6, 2.15409246e5},
{-1.45539891e6, 3.64961835e6,-4.21267423e6, 3.24312555e6,-1.67927904e6, 
   5.23795062e5,-7.66576599e4}};


double const CYrast::elzcof[7][7]={
{ 5.11819909e5,-1.30303186e6, 1.90119870e6,-1.20628242e6, 5.68208488e5,
 5.48346483e4,-2.45883052e4},
{-1.13269453e6,2.97764590e6,-4.54326326e6, 3.00464870e6,-1.44989274e6,
-1.02026610e5, 6.27959815e4},
{ 1.37543304e6,-3.65808988e6,5.47798999e6,-3.78109283e6, 1.84131765e6,
 1.53669695e4,-6.96817834e4},
{-8.56559835e5, 2.48872266e6,-4.07349128e6,3.12835899e6,-1.62394090e6,
 1.19797378e5, 4.25737058e4},
{3.28723311e5,-1.09892175e6, 2.03997269e6,-1.77185718e6,9.96051545e5,
-1.53305699e5,-1.12982954e4},
{ 4.15850238e4,7.29653408e4,-4.93776346e5, 6.01254680e5,-4.01308292e5,
 9.65968391e4,-3.49596027e3},
{-1.82751044e5, 3.91386300e5,-3.03639248e5, 1.15782417e5,-4.24399280e3,
-6.11477247e3, 3.66982647e2}};

double const CYrast::egscof[5][7][5]={
 {{-1.781665232e6,-2.849020290e6, 9.546305856e5, 2.453904278e5, 3.656148926e5},
{ 4.358113622e6, 6.960182192e6,-2.381941132e6,-6.262569370e5,-9.026606463e5},
{-4.804291019e6,-7.666333374e6, 2.699742775e6, 7.415602390e5, 1.006008724e6},
{ 3.505397297e6, 5.586825123e6,-2.024820713e6,-5.818008462e5,-7.353683218e5},
{-1.740990985e6,-2.759325148e6, 1.036253535e6, 3.035749715e5, 3.606919356e5},
{ 5.492532874e5, 8.598827288e5,-3.399809581e5,-9.852362945e4,-1.108872347e5},
{-9.229576432e4,-1.431344258e5, 5.896521547e4, 1.772385043e4, 1.845424227e4}},
 {{ 4.679351387e6, 7.707630513e6,-2.718115276e6,-9.845252314e5,-1.107173456e6},
{-1.137635233e7,-1.870617878e7, 6.669154225e6, 2.413451470e6, 2.691480439e6},
{ 1.237627138e7, 2.030222826e7,-7.334289876e6,-2.656357635e6,-2.912593917e6},
{-8.854155353e6,-1.446966194e7, 5.295832834e6, 1.909275233e6, 2.048899787e6},
{ 4.290642787e6, 6.951223648e6,-2.601557110e6,-9.129731614e5,-9.627344865e5},
{-1.314924218e6,-2.095971932e6, 8.193066795e5, 2.716279969e5, 2.823297853e5},
{ 2.131536582e5, 3.342907992e5,-1.365390745e5,-4.417841315e4,-4.427025540e4}},
{{-3.600471364e6,-5.805932202e6, 1.773029253e6, 4.064280430e5, 7.419581557e5},
{ 8.829126250e6, 1.422377198e7,-4.473342834e6,-1.073350611e6,-1.845960521e6},
{-9.781712604e6,-1.575666314e7, 5.161226883e6, 1.341287330e6,2.083994843e6},
{ 7.182555931e6, 1.156915972e7,-3.941330542e6,-1.108259560e6,-1.543982755e6},
{-3.579820035e6,-5.740079339e6, 2.041827680e6, 5.981648181e5, 7.629263278e5},
{ 1.122573403e6, 1.777161418e6,-6.714631146e5,-1.952833263e5,-2.328129775e5},
 {-1.839672155e5,-2.871137706e5, 1.153532734e5, 3.423868607e4, 3.738902942e4}},
 {{ 2.421750735e6, 4.107929841e6,-1.302310290e6,-5.267906237e5,-6.197966854e5},
{-5.883394376e6,-9.964568970e6, 3.198405768e6, 1.293156541e6, 1.506909314e6},
{ 6.387411818e6, 1.079547152e7,-3.517981421e6,-1.424705631e6,-1.629099740e6},
{-4.550695232e6,-7.665548805e6, 2.530844204e6, 1.021187317e6, 1.141553709e6},
{ 2.182540324e6, 3.646532772e6,-1.228378318e6,-4.813626449e5,-5.299974544e5},
{-6.518758807e5,-1.070414288e6, 3.772592079e5, 1.372024952e5, 1.505359294e5},
 { 9.952777968e4, 1.594230613e5,-6.029082719e4,-2.023689807e4,-2.176008230e4}},
 {{-4.902668827e5,-8.089034293e5, 1.282510910e5,-1.704435174e4, 8.876109934e4},
{ 1.231673941e6, 2.035989814e6,-3.727491110e5, 4.071377327e3,-2.375344759e5},
{-1.429330809e6,-2.376692769e6, 5.216954243e5, 7.268703575e4, 3.008350125e5},
{ 1.114306796e6, 1.868800148e6,-4.718718351e5,-1.215904582e5,-2.510379590e5},
{-5.873353309e5,-9.903614817e5, 2.742543392e5, 9.055579135e4, 1.364869036e5},
{ 1.895325584e5, 3.184776808e5,-9.500485442e4,-3.406036086e4,-4.380685984e4},
{-2.969272274e4,-4.916872669e4, 1.596305804e4, 5.741228836e3, 6.669912421e3}}};


  double const CYrast::aizroc[5][6]={
{ 2.34441624e4,-5.88023986e4, 6.37939552e4,-4.79085272e4, 
 2.27517867e4,-5.35372280e3},
{-4.19782127e4, 1.09187735e5,-1.24597673e5, 9.93997182e4,
-4.95141312e4, 1.19847414e4},
{ 4.18237803e4,-1.05557152e5, 1.16142947e5,-9.00443421e4,
 4.48976290e4,-1.10161792e4},
{-8.27172333e3, 2.49194412e4,-3.39090117e4, 3.33727886e4,
-1.98040399e4, 5.37766241e3},  
{ 5.79695749e2,-1.61762346e3, 2.14044262e3,-3.55379785e3,
    3.25502799e3,-1.15583400e3}};

double const CYrast::ai70c[5][6]={
{ 3.11420101e4,-7.54335155e4, 7.74456473e4,-4.79993065e4,
 2.23439118e4,-4.81961155e3},
{-7.24025043e4, 1.72276697e5,-1.72027101e5, 1.03891065e5,
-4.83180786e4, 1.08040504e4},
{ 7.14932917e4,-1.72792523e5, 1.75814382e5,-1.07245918e5,
 4.86163223e4,-1.10623761e4},
{-2.87206866e4, 6.76667976e4,-6.50167483e4, 3.67161268e4,
-1.74755753e4, 4.67495427e3},
{1.67914908e4,-3.97304542e4, 3.81446552e4,-2.04628156e4,
   7.20091899e3,-1.49978283e3}};

double const CYrast::ai95c[5][6]={
{-6.17201449e5, 1.45561724e6,-1.47514522e6, 9.37798508e5,
-3.74435017e5, 7.81254880e4},
{ 1.24304280e6,-2.94179116e6, 3.00170753e6,-1.92737183e6,
 7.79238772e5,-1.64803784e5},
{-1.49648799e6, 3.52658199e6,-3.56784327e6, 2.26413602e6,
-9.02243251e5, 1.88619658e5},
{ 7.27293223e5,-1.72140677e6,1.75634889e6,-1.12885888e6,
 4.57150814e5,-9.74833991e4},
{-3.75965723e5, 8.83032946e5,-8.87134867e5, 5.58350462e5,
 -2.20433857e5, 4.62178756e4}};


double const CYrast::aimaxc[5][6]={
  {-1.07989556e6, 2.54617598e6,-2.56762409e6, 1.62814115e6,
   -6.39575059e5, 1.34017942e5},
  { 2.17095357e6,-5.13081589e6, 5.19610055e6,-3.31651644e6,
    1.31229476e6,-2.77511450e5},
  {-2.66020302e6, 6.26593165e6,-6.31060776e6, 3.99082969e6, 
   -1.56447660e6, 3.25613262e5},
  { 1.29464191e6,-3.05746938e6, 3.09487138e6,-1.97160118e6,
    7.79696064e5,-1.63704652e5},
  {-7.13073644e5, 1.67482279e6,-1.67984330e6, 1.05446783e6,
   -4.10928559e5, 8.43774143e4}};

double const CYrast::ai952c[5][6]={
  {-7.37473153e5, 1.73682827e6,-1.75850175e6, 1.11320647e6, 
   -4.41842735e5, 9.02463457e4},
  { 1.49541980e6,-3.53222507e6, 3.59762757e6,-2.29652257e6,
    9.21077757e5,-1.90079527e5},
  {-1.80243593e6, 4.24319661e6,-4.29072662e6, 2.71416936e6,
   -1.07624953e6, 2.20863711e5},
  { 8.86920591e5,-2.09589683e6,2.13507675e6,-1.36546686e6,
    5.48868536e5,-1.14532906e5},
  {-4.62131503e5, 1.08555722e6,-1.09187524e6, 6.87308217e5,
   -2.70986162e5, 5.61637883e4}};

double const CYrast::aimax2c[5][6]={
  {-1.16343311e6, 2.74470544e6,-2.77664273e6, 1.76933559e6,
   -7.02900226e5, 1.49345081e5},
  { 2.36929777e6,-5.60655122e6,5.70413177e6,-3.66528765e6,
    1.47006527e6,-3.15794626e5},
  {-2.82646077e6, 6.66086824e6,-6.72677653e6, 4.27484625e6,
   -1.69427298e6, 3.58429081e5},
  { 1.39112772e6,-3.29007553e6, 3.34544584e6,-2.14723142e6,
    8.61118401e5,-1.84500129e5},
  {-7.21329917e5, 1.69371794e6,-1.69979786e6, 1.07037781e6,
   -4.20662028e5, 8.80728361e4}};

double const CYrast::aimax3c[4][4]={
  {-2.88270282e3, 5.30111305e3,-3.07626751e3,6.56709396e2},
  {5.84303930e3,-1.07450449e4,6.24110631e3,-1.33480875e3},
  {-4.20629939e3,7.74058373e3,-4.50256063e3, 9.65788439e2},
  {1.23820134e3,-2.28228958e3, 1.33181316e3,-2.87363568e2}};

double const CYrast::aimax4c[4][4]={
  {-3.34060345e3, 6.26384099e3,-3.77635848e3, 8.57180868e2},
  { 6.76377873e3,-1.26776571e4, 7.64206952e3,-1.73406840e3},
  {-4.74821371e3, 8.89857519e3,-5.36266252e3, 1.21614216e3},
  { 1.46369384e3,-2.74251101e3, 1.65205435e3,-3.74262365e2}};

double const CYrast::bizroc[4][6]={
  { 5.88982505e2,-1.35630904e3, 1.32932125e3,-7.78518395e2,
    2.73122883e2,-3.49600841e1},
  {-9.67701343e2, 2.24594418e3,-2.24303790e3, 1.35440047e3,
   -4.96538939e2, 6.66791793e1},
  { 1.17090267e3,-2.71181535e3, 2.67008958e3,-1.58801770e3,
    5.66896359e2,-8.21530057e1},
  {-3.83031864e2, 9.05191483e2,-9.30560410e2, 5.96618532e2,
   -2.34403480e2, 3.97909172e1}};

double const CYrast::bi70c[4][6]={
  { 2.32414810e3,-5.42381778e3, 5.40202710e3,-3.26923144e3,
    1.18318943e3,-1.93186467e2},
  {-4.38084778e3, 1.03523570e4,-1.05573803e4, 6.59901160e3,
   -2.47601209e3, 4.19497260e2},
  { 4.35377377e3,-1.01728647e4, 1.01311246e4,-6.14038462e3,
    2.21957562e3,-3.62854365e2},
  {-1.84533539e3, 4.41613298e3,-4.59403284e3, 2.95951225e3,
   -1.14630148e3, 2.02702459e2}};

double const CYrast::bi95c[4][6]={
  { 1.55359266e3,-3.58209715e3, 3.50693744e3,-2.03992913e3,
    7.05498010e2,-1.49075519e2},
  {-2.86876240e3, 6.77107086e3,-6.90300614e3, 4.20246063e3,
   -1.50290693e3, 3.13662258e2},
  { 2.60138185e3,-5.95414919e3, 5.70261588e3,-3.17188958e3,
    9.89207911e2,-1.76320647e2},
  {-1.75198402e3, 4.16635208e3,-4.25212424e3, 2.59953301e3,
   -9.09813362e2, 1.51070448e2}};
double const CYrast::bimaxc[4][6]={
  { 4.17708254e3,-8.59358778e3, 6.46392215e3,-8.84972189e2,
    -1.59735594e3, 1.39662071e3},
  {-1.56318394e4, 3.54574417e4,-3.35945173e4, 1.65495998e4,
   -3.32021998e3,-1.46150905e3},
  { 1.41292811e4,-3.11818487e4, 2.77454429e4,-1.19628827e4,
    1.28008968e3, 1.66111636e3},
  {-1.92878152e4, 4.56505796e4,-4.66413277e4, 2.89229633e4,
   -1.07284346e4, 1.50513815e3}};

double const CYrast::b[8][5][5]={
 { {-17605.486,    7149.518,   23883.053,   21651.272,    -348.402}, 
{ 67617.841,  -43739.724,  -72255.975,  -79130.794,  -44677.571},
{-38284.438,    1629.405,   49928.158,   41804.076,   -4744.661},
{ 32124.997,  -30258.577,  -45416.684,  -40347.682,  -20718.220},
{ -7009.141,  -12274.997,    5703.797,    4375.711,   -2519.378}},
{{ 44633.464,  -19664.988,  -59871.761,  -55426.914,    -834.542},
{-168222.448,  111291.480,  179939.092,  198092.569,  111347.893},
{ 97361.913,   -9404.674, -125757.922, -108003.071,    7218.634},
{-79483.915,   77325.943,  113433.935,  100793.473,   51134.715},
{ 17974.287,   28084.202,  -15187.024,  -12244.758,    4662.638}},
{{-52775.944,   26372.246,   69061.418,   66412.802,    4456.150},
{191743.968, -131672.398, -204896.966, -227867.388, -126609.907},
{-115704.145,   22053.324,  146315.778,  131364.763,    1176.787},
{  89552.374,  -92216.043, -129901.381, -115369.145,  -56853.433},
{ -21565.916,  -26595.626,   19490.177,   16758.541,   -2102.079}},
{{  43089.492,  -24205.277,  -54145.651,  -54916.885,   -6982.880},
{-148893.417,  105958.013,  157835.761,  178691.634,   97246.541},
{  94988.919,  -28100.912, -115999.016, -110379.031,  -10165.655},
{ -68277.650,   75226.302,  100957.904,   89635.423,   41887.718},
{   17720.511,   14837.634,  -17476.245,  -15762.519,   -1297.194}},
{{-25673.002,   15584.030,   30291.575,   33008.500,    6161.347},
{  83671.721,  -60540.161,  -86493.251, -101211.404,  -53525.729},
{ -56872.290,   21837.845,   65632.326,   67209.922,   11367.245},
{  37417.011,  -43832.496,  -56023.110,  -49995.410,  -21521.767},
{ -10406.752,   -4483.654,   11314.537,   10435.311,    2206.541}},
{{  11115.697,   -6672.659,  -11696.751,  -14183.580,   -3586.280},
{ -33914.223,   23453.203,   32630.486,   40847.174,   21223.389},
{  24704.794,  -10226.568,  -25497.313,  -29082.833,   -7190.865},
{ -14734.982,   17457.735,   21457.946,   19706.266,    7745.219},
{   4283.242,     214.103,   -5058.390,   -4715.987,   -1303.557}},
{{  -3196.271,    1579.967,    2718.777,    3997.275,    1361.541},
{   9152.527,   -5182.345,   -7412.465,  -10822.583,   -5804.143},
{  -7177.709,    2317.015,    5874.736,    8208.479,    2858.079},
{   3919.329,   -4068.562,   -4982.296,   -5065.564,   -1948.636},
{  -1147.488,     265.063,    1372.992,    1333.320,     416.369}},
{{    442.331,    -219.051,    -337.663,    -575.255,    -232.634},
{  -1217.351,     590.169,     854.438,    1468.791,     811.481},
{   1031.157,    -227.306,    -651.183,   -1174.558,    -528.960},
{   -536.534,     413.429,     533.977,     663.064,     278.203},
{      154.548,     -54.623,    -169.817,    -183.006,     -63.999}}};



double const CYrast::hbarc=197.32858;
double const CYrast::alfinv=137.035982;
double const CYrast::srznw=1.16;
double const CYrast::aknw=2.3;
double const CYrast::bb=0.99;
double const CYrast::um=931.5016;
double const CYrast::elm=0.51104;
double const CYrast::spdlt=2.99792458e23;
double const CYrast::asnw=21.13;
double const CYrast::kx[8]={0.0,0.08299,0.1632,0.2435,0.3217,0.402
			    ,0.51182,0.617422};
double const CYrast::ky[6]={0.0,0.02,0.05,0.10,0.15,0.20};
double const CYrast::ka[11]={0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.};
double const CYrast::r0=1.225;
double const CYrast::sep=2.;



//******************************************************************
/**
 * Constructor
 */
CYrast::CYrast()
{

  mass = CMass::instance(); // mass singleton
  QString fileName(":/tbl/sad.tbl");

  QFile iff(fileName);
  if(!iff.open(QIODevice::ReadOnly | QIODevice::Text))
    {
     cout << "could not open " << fileName.toStdString() << " in CYrast" << endl;
     abort();

  }

  QTextStream ifFile(&iff);

  QString line;
  line = ifFile.readLine();
  //getline(ifFile,line);
 // qDebug() << "line = " << line << endl;
  for (int i5=0;i5<2;i5++)
    for (int i4=0;i4<11;i4++)
      for (int i3=0;i3<2;i3++)
	{
            //getline(ifFile,line);
          ifFile.readLine();
            //cout << line << endl;
            //cout << i3+1 << " " << i4+1 << " " << i5+1 << endl;
	    for (int i2=0;i2<8;i2++)
	      {
              for(int i1=0;i1<6;i1++)ifFile >> c[i1][i2][i3][i4][i5];
              ifFile.readLine();
             // getline(ifFile,line);
	      }
	}
  ifFile.flush();
  iff.close();
}

CYrast* CYrast::instance()
{
    if (fInstance == 0) {
        fInstance = new CYrast;
    }
    return fInstance;
}

//*********************************************************************
/**
 * Returns the yrast energy in MeV from the Rotating Liquid Drop Model
 * \param iZ is the proton number
 * \param iA is the mass number
 * \param fL is the angular momentum in hbar
 */
double CYrast::getYrastRLDM(int iZ, int iA, double fL)
{
    const double fZ = (double)iZ;
    const double fA = (double)iA;
    const double fN = fA - fZ;

    const double asym = (fN - fZ) / fA;
    double paren = 1.0f - 1.7826f * (asym * asym);
    if (paren <= 0.0f) paren = 0.1f;

    // Replace pow(fA, 2/3), pow(fA, 5/3), pow(fA, 7/3) using cbrt
    const double a13 = cbrt(fA);          // A^(1/3)
    const double a23 = a13 * a13;          // A^(2/3)
    const double a53 = fA * a23;           // A^(5/3) = A * A^(2/3)
    const double a73 = (fA * fA) * a13;    // A^(7/3) = A^2 * A^(1/3)

    const double eso = 17.9439f * paren * a23;
    const double x   = 0.019655f * fZ * (fZ / fA) / paren;

    const double fL2 = fL * fL;
    const double ero = 34.548f * fL2 / a53;
    const double y   = 1.9254f * fL2 / (paren * a73);

    int ix = (int)(20.0f * x + 1.0f);
    const double bx = 20.0f * x + 1.0f;
    const double dx = bx - (double)ix;

    double hf = 0.0f;

    if (x <= 0.25f)
    {
        double by = 10.0f * y + 1.0f;
        if (by > 9.0f) by = 9.0f;
        if (by < 1.0f) by = 1.0f;

        const int iy = (int)by;
        const double dy = by - (double)iy;

        const double h1 = (x1h[iy - 1][ix] - x1h[iy - 1][ix - 1]) * dx + x1h[iy - 1][ix - 1];

        // NOTE: kept your original expression (it has x1h[iy][ix]-x1h[iy][ix] which is 0)
        const double h2 = (x1h[iy][ix] - x1h[iy][ix]) * dx + x1h[iy][ix];

        hf = (h2 - h1) * dy + h1;
    }
    else if (x < 0.5f)
    {
        double by = 20.0f * y + 1.0f;
        if (by > 11.0f) by = 10.0f;
        if (by < 1.0f)  by = 1.0f;

        ix = ix - 5;

        const int iy = (int)by;
        const double dy = by - (double)iy;

        const double h1 = (x2h[iy - 1][ix] - x2h[iy - 1][ix - 1]) * dx + x2h[iy - 1][ix - 1];
        const double h2 = (x2h[iy][ix]     - x2h[iy][ix - 1])     * dx + x2h[iy][ix - 1];

        hf = (h2 - h1) * dy + h1;
    }
    else
    {
        double xx = x;
        if (xx > 0.94999999f) xx = 0.949999999f;

        ix = (int)(20.0f * xx + 1.0f);
        ix = ix - 10;

        double by = 100.0f * y + 1.0f;
        if (by > 19.0f) by = 19.0f;
        if (by < 1.0f)  by = 1.0f;

        const int iy = (int)by;
        const double dy = by - (double)iy;

        const double h1 = (x3h[iy - 1][ix] - x3h[iy - 1][ix - 1]) * dx + x3h[iy - 1][ix - 1];
        const double h2 = (x3h[iy][ix]     - x3h[iy][ix - 1])     * dx + x3h[iy][ix - 1];

        hf = (h2 - h1) * dy + h1;
    }

    return ero + hf * eso;
}

//*************************************************
/**
 * calculates a array of Legendre Polynomials
 * \param x is the independent variable
 * \param n is the maximum order of the array
 * \param pl is the pointer to the array
 */
void CYrast::lpoly(double x, int n, double* pl)
{
    pl[0] = 1.0;
    pl[1] = x;

    for (int i = 2; i < n; i++)
    {
        pl[i] = ((double)(2 * i - 1) * x * pl[i - 1] - (double)(i - 1) * pl[i - 2]) / (double)i;
    }
}

//************************************************
/**
 * Returns the relative surface area of the saddle-point shape
 * from the model of Sierk. The surface is relative to the
 * value for a spherical nucleus
 * \param J is the angular momentum
 */
double CYrast::getBsSierk(double J)
{
    const double xa = A / 320.0;
    const double JJ = (double)J / Jmax;

    double pl[9];
    lpoly(JJ, 9, pl);

    double paa[5];
    lpoly(xa, 5, paa);

    double pzz[8];
    lpoly(zz, 8, pzz);

    double bs = 0.0;

    for (int izz = 0; izz < 8; izz++)
        for (int iaa = 0; iaa < 5; iaa++)
            for (int ill = 0; ill < 5; ill++)
                bs += b[izz][iaa][ill] * pzz[izz] * paa[iaa] * pl[2 * ill];

    if (bs < 1.0) bs = 1.0;
    return (double)bs;
}

//**************************************
/**
 * Returns the Maximum angular momentum where the nucleus is stable
 * against fission in the model of Sierk.
 * Units are in hbar.
 * \param iZ is proton number
 * \param iA is the mass number
 */
double CYrast::getJmaxSierk(int iZ, int iA)
{
    if (iZ < 19 || iZ > 111)
    {
        cout << " Z out of range for CYrast::Jmax" << endl;
        abort();
    }

    Z = (double)iZ;
    A = (double)iA;

    // replace pow(Z,2) with Z*Z
    const double Z2 = Z * Z;
    amin = 1.2 * Z + 0.01 * Z2;
    amax = 5.8 * Z - 0.024 * Z2;

    if (A < amin || A > amax)
    {
        cout << "A out of limits for CYrast::Jmax" << endl;
        abort();
    }

    const double aa = 2.5e-3 * A;
    zz = 1.0e-2 * Z;

    Jmax = 0.0;
    lpoly(zz, 7, pz);
    lpoly(aa, 7, pa);

    for (int i = 0; i < 5; i++)
        for (int k = 0; k < 7; k++)
            Jmax += emxcof[i][k] * pz[k] * pa[i];

    return (double)Jmax;
}

//******************************************
/**
 * Returns the fission barrier from the model of Sierk. Units are in MeV.
 * The function getJmaxSierk() must be called beforehand
 * \param J is the angular momentum in hbar
 */
double CYrast::getBarrierFissionSierk(double J)
{
    const double dJ = (double)J;

    if (Z > 102.0 && J > 0.0f)
    {
        cout << "z and J out of limits" << endl;
        return 999.0f;
    }

    double bfis = 0.0;
    for (int i = 0; i < 7; i++)
        for (int k = 0; k < 7; k++)
            bfis += elzcof[i][k] * pz[k] * pa[i];

    if (dJ < 1.0) return (double)bfis;

    double J80 = 0.0;
    double J20 = 0.0;

    for (int i = 0; i < 4; i++)
        for (int k = 0; k < 5; k++)
        {
            J80 += elmcof[i][k] * pz[k] * pa[i];
            J20 += emncof[i][k] * pz[k] * pa[i];
        }

    if (dJ <= J20)
    {
        const double J20_2 = J20 * J20;
        const double J80_2 = J80 * J80;
        const double J80_3 = J80_2 * J80;
        const double J20_3 = J20_2 * J20;

        const double q  = 0.2 / (J20_2 * J80_2 * (J20 - J80));
        const double qa = q * (4.0 * J80_3 - J20_3);
        const double qb = -q * (4.0 * J80_2 - J20_2);

        const double dJ2 = dJ * dJ;
        const double dJ3 = dJ2 * dJ;

        bfis *= (1.0 + qa * dJ2 + qb * dJ3);
    }
    else
    {
        const double x = (double)(J20 / Jmax);
        const double y = (double)(J80 / Jmax);

        const double x2 = x * x;
        const double x4 = x2 * x2;
        const double x5 = x4 * x;

        const double y2 = y * y;
        const double y4 = y2 * y2;
        const double y5 = y4 * y;

        const double ym1 = (y - 1.0f);
        const double xm1 = (x - 1.0f);

        const double aj = (-20.0f * x5 + 25.0f * x4 - 4.0f) * (ym1 * ym1) * y2;
        const double ak = (-20.0f * y5 + 25.0f * y4 - 1.0f) * (xm1 * xm1) * x2;

        const double oneMinusX = (1.0f - x);
        const double oneMinusY = (1.0f - y);
        const double t = (oneMinusX * oneMinusY * x * y);
        const double q = 0.2f / ((y - x) * (t * t));

        const double qa = q * (aj * y - ak * x);
        const double qb = -q * (aj * (2.0f * y + 1.0f) - ak * (2.0f * x + 1.0f));

        const double zzz = (double)(dJ / Jmax);
        const double z2 = zzz * zzz;
        const double z4 = z2 * z2;
        const double z5 = z4 * zzz;

        const double a1 = 4.0f * z5 - 5.0f * z4 + 1.0f;
        const double a2 = qa * (2.0f * zzz + 1.0f);

        const double zm1 = (zzz - 1.0f);
        bfis *= (double)(a1 + zm1 * (a2 + qb * zzz) * z2 * zm1);
    }

    if (bfis <= 0.0) bfis = 0.0;
    if (J > (double)Jmax) bfis = 0.0;

    return (double)bfis;
}

//*****************************************
/**
 * Returns the yrast energy in the model of Sierk.
 * Units are in MeV. The function getJmaxSierk must be called beforehand.
 * \param J is the angular momentum
 */
double CYrast::getYrastSierk(double J)
{
    if (Z > 102.0 && J > 0.0f)
    {
        cout << "z and J out of limits" << endl;
        return 999.0f;
    }

    if (A < amin || A > amax)
    {
        cout << "A is out of range" << endl;
        return 999.0f;
    }

    double Erot = 0.0;
    const double JJ = (double)J / Jmax;

    double pl[9];
    lpoly(JJ, 9, pl);

    for (int k = 0; k < 5; k++)
        for (int l = 0; l < 7; l++)
            for (int m = 0; m < 5; m++)
                Erot += egscof[k][l][m] * pz[l] * pa[k] * pl[2 * m];

    if (Erot < 0.0) Erot = 0.0;
    return (double)Erot;
}

//*****************************************************
/**
 * Calculates the three principal moments of inertia
 * associated with the saddle-point configuration in the model
 * of Sierk. The function getJmaxSierk must be called beforehand.
 * \param J is the angular momentum
 */
double CYrast::getMomentOfInertiaSierk(double J)
{
    const double JJ = (double)J / Jmax;

    double aizro  = 0.0;
    double ai70   = 0.0;
    double aimax  = 0.0;
    double ai95   = 0.0;
    double bizro  = 0.0;
    double bi70   = 0.0;
    double bimax  = 0.0;
    double bi95   = 0.0;
    double aimax2 = 0.0;
    double ai952  = 0.0;

    for (int l = 0; l < 6; l++)
    {
        for (int k = 0; k < 5; k++)
        {
            aizro  += aizroc[k][l]  * pz[l] * pa[k];
            ai70   += ai70c[k][l]   * pz[l] * pa[k];
            ai95   += ai95c[k][l]   * pz[l] * pa[k];
            aimax  += aimaxc[k][l]  * pz[l] * pa[k];
            ai952  += ai952c[k][l]  * pz[l] * pa[k];
            aimax2 += aimax2c[k][l] * pz[l] * pa[k];
        }

        for (int k = 0; k < 4; k++)
        {
            bizro += bizroc[k][l] * pz[l] * pa[k];
            bi70  += bi70c[k][l]  * pz[l] * pa[k];
            bi95  += bi95c[k][l]  * pz[l] * pa[k];
            bimax += bimaxc[k][l] * pz[l] * pa[k];
        }
    }

    double ff1 = 1.0, ff2 = 0.0;
    double fg1 = 1.0, fg2 = 0.0;
    double aimidh = 0.0;

    if (Z > 70.0)
    {
        double aimaxh = 0.0;
        aimidh = 0.0;

        for (int l = 0; l < 4; l++)
            for (int k = 0; k < 4; k++)
            {
                aimaxh += aimax3c[k][l] * pz[l] * pa[k];
                aimidh += aimax4c[k][l] * pz[l] * pa[k];
            }

        if (Z > 80.0) ff1 = 0.0;
        if (Z >= 80.0) fg1 = 0.0;
        if (bimax > 0.95) fg1 = 0.0;
        if (aimaxh > aimax) ff1 = 0.0;

        ff2 = 1.0 - ff1;
        fg2 = 1.0 - fg1;

        aimax  = aimax  * ff1 + ff2 * aimaxh;
        aimax2 = aimax2 * ff1 + ff2 * aimidh;
    }

    bimax = bimax * fg1 + aimidh * fg2;

    const double saizro  = max(aizro,  0.0);
    const double sai70   = max(ai70,   0.0);
    const double sai95   = max(ai95,   0.0);
    const double saimax  = max(aimax,  0.0);
    const double sai952  = max(ai952,  0.0);
    const double simax2  = max(aimax2, 0.0);
    const double sbimax  = max(bimax,  0.0);
    const double sbi70   = max(bi70,   0.0);
    const double sbi95   = max(bi95,   0.0);
    const double sbizro  = max(bizro,  0.0);

    double q1 = -3.148849569;
    double q2 =  4.465058752;
    double q3 = -1.316209183;
    const double q4 =  2.26129233;
    const double q5 = -4.94743352;
    const double q6 =  2.68614119;

    const double gam  = -20.0 * log(fabs(saizro - sai95)  / fabs(saizro - saimax));
    const double aa   = q1 * saizro + q2 * sai70 + q3 * sai95;
    const double bb   = q4 * saizro + q5 * sai70 + q6 * sai95;

    const double gam2 = -20.0 * log(fabs(saizro - sai952) / fabs(saizro - simax2));
    const double aa2  = q1 * saizro + q2 * sai70 + q3 * sai952;
    const double bb2  = q4 * saizro + q5 * sai70 + q6 * sai952;

    const double aa3  = q1 * sbizro + q2 * sbi70 + q3 * sbi95;
    const double bb3  = q4 * sbizro + q5 * sbi70 + q6 * sbi95;

    const double gam3 = 60.0;

    const double alpha = pi * (JJ - 0.7);
    const double beta  = 5.0 * pi * (JJ - 0.9);

    const double JJ2 = JJ * JJ;
    const double JJ4 = JJ2 * JJ2;

    const double silt  = saizro + aa  * JJ2 + bb  * JJ4;
    const double sjlt  = sbizro + aa3 * JJ2 + bb3 * JJ4;
    const double silt2 = saizro + aa2 * JJ2 + bb2 * JJ4;

    if (JJ <= 0.70)
    {
        momInertiaMin = sjlt;
        momInertiaMax = silt;
        momInertiaMid = silt2;
    }
    else
    {
        const double sigt  = saizro + (saimax - saizro) * exp(gam  * (JJ - 1.0));
        const double sjgt  = sbi95  + (sbimax - sbi95)  * exp(gam3 * (JJ - 1.0));
        const double sigt2 = saizro + (simax2 - saizro) * exp(gam2 * (JJ - 1.0));

        if (JJ <= 0.95)
        {
            const double ca = cos(alpha), sa = sin(alpha);
            const double ca2 = ca * ca,  sa2 = sa * sa;

            momInertiaMax = silt  * ca2 + sigt  * sa2;
            momInertiaMin = sjlt  * ca2 + sjgt  * sa2;
            momInertiaMid = silt2 * ca2 + sigt2 * sa2;
        }
        else
        {
            const double cb = cos(beta), sb = sin(beta);
            const double cb2 = cb * cb, sb2 = sb * sb;

            momInertiaMax = silt  * cb2 + sigt  * sb2;
            momInertiaMin = sjlt  * cb2 + sjgt  * sb2;
            momInertiaMid = silt2 * cb2 + sigt2 * sb2;
        }
    }

    if (ff2 > 0.01 && fg2 > 0.01)
    {
        q1 = 4.001600640;
        q2 = 0.960784314;
        q3 = 2.040816327;

        const double aa3n =  q1 * sai70 - q2 * saimax - (1.0 + q3) * saizro;
        const double bb3n = -q1 * sai70 + (1.0 + q2) * saimax + q3 * saizro;

        const double aa4n =  q1 * sai70 - q2 * simax2 - (1.0 + q3) * saizro;
        const double bb4n = -q1 * sai70 + (1.0 + q2) * simax2 + q3 * saizro;

        momInertiaMax = saizro + aa3n * JJ2 + bb3n * JJ4;
        momInertiaMid = saizro + aa4n * JJ2 + bb4n * JJ4;
    }

    if (momInertiaMid > momInertiaMax) momInertiaMid = momInertiaMax;
    momInertiaMid = max(momInertiaMid, (double)0.0);

    // momInertiaSphere = 0.4 * r0^2 * A^(5/3)
    const double A13 = cbrt(A);
    const double A23 = A13 * A13;
    const double A53 = A * A23;
    const double momInertiaSphere = (double)(0.4 * (r0 * r0) * A53);

    momInertiaMax *= momInertiaSphere;
    momInertiaMid *= momInertiaSphere;
    momInertiaMin *= momInertiaSphere;

    return momInertiaMax;
}

//*******************************************************************
/**
 * Returns the fission barrier in MeV from the Rotating Liquid Drop Model.
 * \param iZ is the proton number.
 * \param iA is the mass number
 * \param fJ is the angular momentum
 */
double CYrast::getBarrierFissionRLDM(int iZ, int iA, double fJ)
{
    const double A = (double)iA;
    const double Z = (double)iZ;
    const double N = A - Z;

    const double asym = (N - Z) / A;
    double paren = 1.0f - 1.7826f * (asym * asym);

    const double a13 = cbrt(A);
    const double a23 = a13 * a13;
    const double a73 = (A * A) * a13; // A^(7/3)

    const double eso = 17.9439f * paren * a23;
    const double x   = 0.019655f * Z * (Z / A) / paren;

    const double fJ2 = fJ * fJ;
    const double y   = 1.9254f * fJ2 / (paren * a73);

    int ix = (int)(20.0f * x + 0.999f);
    const double bx = 20.0f * x + 0.999f;
    const double dx = bx - (double)ix;

    double bf = 0.0f;

    if (x <= 0.25f)
    {
        double by = 10.0f * y + 0.999f;
        if (by > 9.0f) by = 9.0f;
        if (by < 1.0f) by = 1.0f;

        const int iy = (int)by;
        const double dy = by - (double)iy;

        const double b2 = (x1b[iy][ix]     - x1b[iy][ix - 1])     * dx + x1b[iy][ix - 1];
        const double b1 = (x1b[iy - 1][ix] - x1b[iy - 1][ix - 1]) * dx + x1b[iy - 1][ix - 1];

        bf = (b2 - b1) * dy + b1;
    }
    else if (x <= 0.5f)
    {
        double by = 20.0f * y + 0.999f;
        if (by > 11.0f) by = 11.0f;
        if (by < 1.0f)  by = 1.0f;

        ix = ix - 5;

        const int iy = (int)by;
        const double dy = by - (double)iy;

        const double b1 = (x2b[iy - 1][ix] - x2b[iy - 1][ix - 1]) * dx + x2b[iy - 1][ix - 1];
        const double b2 = (x2b[iy][ix]     - x2b[iy][ix - 1])     * dx + x2b[iy][ix - 1];

        bf = (b2 - b1) * dy + b1;
    }
    else
    {
        double xx = x;
        if (xx > 0.95f) xx = 0.95f;

        ix = (int)(20.0f * xx + 0.999f);
        ix = ix - 10;

        double by = 100.0f * y + 0.999f;
        if (by > 19.0f) by = 19.0f;
        if (by < 1.0f)  by = 1.0f;

        const int iy = (int)by;
        const double dy = by - (double)iy;

        const double b1 = (x3b[iy - 1][ix] - x3b[iy - 1][ix - 1]) * dx + x3b[iy - 1][ix - 1];
        const double b2 = (x3b[iy][ix]     - x3b[iy][ix - 1])     * dx + x3b[iy][ix - 1];

        bf = (b2 - b1) * dy + b1;
    }

    return bf * eso;
}

//*****************************************************
/**
 * Cubic polynomial used in 2D cubic spline interpolation
 */
double CYrast::cubic(double a, double b, double c, double d, double e, double f)
{
    // no pow here already
    return a + f * (e * c + f * (3.0f * (b - a) - (d + 2.0f * c) * e
                                 + f * (2.0f * (a - b) + (c + d) * e)));
}

//*****************************************************
/**
 * Prepares for the calculation of conditional saddle-point
 * energies in MeV
 * \param iZ0 is the proton number
 * \param iA0 is the mass number
 * \param fJ0 is the angular momentum in hbar
 */
void CYrast::prepareAsyBarrier(int iZ0, int iA0, double fJ0)
{
    iZ = iZ0;
    iA = iA0;
    fJ = fJ0;

    const double A = (double)iA;
    const double Z = (double)iZ;

    Narray = (int)(A / 2.0f);

    int ia = -1;

    const double ac = 0.6f * hbarc / alfinv / srznw;

    // pow(A,1/3) and pow(A,2/3)
    const double A13 = cbrt(A);
    const double A23 = A13 * A13;

    const double rz = srznw * A13;

    const double asym = (A - 2.0f * Z) / A;
    double cs = asnw * (1.0f - aknw * (asym * asym));

    // (bb/(1.4142*srznw))^3  and (Z/A)^2
    const double tbb = bb / (1.4142f * srznw);
    const double tbb3 = tbb * tbb * tbb;

    const double ZoverA = Z / A;
    const double ZoverA2 = ZoverA * ZoverA;

    const double delcs =
        45.0f * hbarc / (8.0f * srznw * alfinv) * tbb3 * ZoverA2;

    cs = cs + delcs;

    const double esz = cs * A23;

    const double emz = um * A - elm * Z;

    // zsqoa = Z^2 / A (avoid pow)
    const double zsqoa = (Z * Z) / A;

    double x = ac * zsqoa / (2.0f * cs);

    if (x > 0.61743f)
    {
        if (first)
        {
            //cout << "No barriers available for this nucleus" << endl;
            //cout << "Z= " << Z << " A= " << A << endl;
            //cout << "using scaled 194Hg barriers" << endl;
            first = 0;
        }
        x = 0.61743f;
    }

    //----find y fissility parameter-------------------------
    if (esz == 0.0f)
    {
        cout << "esz=0 in sad, Z= " << Z << " A= " << A << endl;
        abort();
    }

    if (emz / esz <= 0.0f)
    {
        cout << "in yrast.PrepareAsyBar square root of negative for A= "
             << A << " Z=" << Z << endl;
        cout << " emz= " << emz << " esz= " << esz << endl;
        abort();
    }

    const double tz = rz * sqrt(emz / esz) / spdlt;
    const double elz = esz * tz;
    const double elzohb = elz / 6.582173e-22f;

    // y = 1.25 * (fJ/elzohb)^2
    const double ratio = fJ / elzohb;
    double y = 1.25f * (ratio * ratio);

    //----------find normalization factor-------------------
    double sadf;

    if (Z >= 19.0f)
    {
        // polynomial in x: replace pow(x,n) with multiplies
        const double x2 = x * x;
        const double x3 = x2 * x;
        const double x4 = x2 * x2;
        const double x5 = x4 * x;

        const double sad0 =
            -5.5286e-2f + 190.03f * x + 235.189f * x2
            - 2471.249f * x3 + 4266.905f * x4 - 2576.048f * x5;

        double bfis;
        if (x >= 0.6174f)
            bfis = 13.86f; // more fissile than 194Hg
        else
        {
            getJmaxSierk(iZ, iA); // initialize sierk sub
            bfis = getBarrierFissionSierk(0.0f);
        }

        sadf = bfis / sad0;
    }
    else
    {
        const double x2 = x * x;
        const double x3 = x2 * x;

        const double ess = 2376.0f * x - 6062.4f * x2 + 8418.3f * x3;

        // sadf = (esz/ess)^2.12  -> keep pow (non-integer exponent)
        sadf = pow(esz / ess, 2.12f);
    }

    //-----select subroutines for cubic-spline interpolation--

    // find nearest knots for interpolation
    int iy = 1;
    bool yExtrapolation = false;
    double yy = 0.0f;

    if (y > ky[5])
    {
        iy = 5;
        yExtrapolation = true;
        yy = y;
        y = ky[5];
    }
    else
    {
        for (;;)
        {
            if (y <= ky[iy]) break;
            iy++;
            if (iy > 5) break;
        }
        iy = iy - 1;
    }

    int ix = 1;
    for (;;)
    {
        if (x <= kx[ix]) break;
        ix++;
        if (ix > 6) break;
    }
    ix = ix - 1;

    // now y is between knots ky[iy] & ky[iy+1] and x is
    // between knots kx[ix] & kx[ix+1]

    const double Dkx = kx[ix + 1] - kx[ix];
    const double Rx  = (x - kx[ix]) / Dkx;

    const double Dky = ky[iy + 1] - ky[iy];
    const double Ry  = (y - ky[iy]) / Dky;

    double Dka = 0.0f;
    double C0 = 0.0f;
    double C1 = 0.0f;
    double C2 = 0.0f;
    double C3 = 0.0f;

    //--- fill array with saddle point energies-------
    for (int ii = Narray; ii > 0; ii--)
    {
        const double alpha = (A - 2.0f * (double)ii) / A;

        if (alpha >= ka[ia + 1])
        {
            for (;;)
            {
                if (alpha < ka[ia + 1]) break;
                ia++;
            }

            double k1 = cubic(c[iy][ix][0][ia][0],     c[iy][ix + 1][0][ia][0],
                             c[iy][ix][1][ia][0],     c[iy][ix + 1][1][ia][0], Dkx, Rx);

            double k2 = cubic(c[iy][ix][0][ia + 1][0], c[iy][ix + 1][0][ia + 1][0],
                             c[iy][ix][1][ia + 1][0], c[iy][ix + 1][1][ia + 1][0], Dkx, Rx);

            double k3 = cubic(c[iy][ix][0][ia][1],     c[iy][ix + 1][0][ia][1],
                             c[iy][ix][1][ia][1],     c[iy][ix + 1][1][ia][1], Dkx, Rx);

            double k4 = cubic(c[iy][ix][0][ia + 1][1], c[iy][ix + 1][0][ia + 1][1],
                             c[iy][ix][1][ia + 1][1], c[iy][ix + 1][1][ia + 1][1], Dkx, Rx);

            Dka = ka[ia + 1] - ka[ia];

            C0 = k1;
            C1 = Dka * k3;
            C2 = 3.0f * (k2 - k1) - (k4 + 2.0f * k3) * Dka;
            C3 = 2.0f * (k1 - k2) + (k3 + k4) * Dka;

            if (fJ > 0.1f)
            {
                k1 = cubic(c[iy + 1][ix][0][ia][0],     c[iy + 1][ix + 1][0][ia][0],
                           c[iy + 1][ix][1][ia][0],     c[iy + 1][ix + 1][1][ia][0], Dkx, Rx);

                k2 = cubic(c[iy + 1][ix][0][ia + 1][0], c[iy + 1][ix + 1][0][ia + 1][0],
                           c[iy + 1][ix][1][ia + 1][0], c[iy + 1][ix + 1][1][ia + 1][0], Dkx, Rx);

                k3 = cubic(c[iy + 1][ix][0][ia][1],     c[iy + 1][ix + 1][0][ia][1],
                           c[iy + 1][ix][1][ia][1],     c[iy + 1][ix + 1][1][ia][1], Dkx, Rx);

                k4 = cubic(c[iy + 1][ix][0][ia + 1][1], c[iy + 1][ix + 1][0][ia + 1][1],
                           c[iy + 1][ix][1][ia + 1][1], c[iy + 1][ix + 1][1][ia + 1][1], Dkx, Rx);

                C0 = (k1 - C0) * Ry + C0;
                C1 = (Dka * k3 - C1) * Ry + C1;
                C2 = (3.0f * (k2 - k1) - (k4 + 2.0f * k3) * Dka - C2) * Ry + C2;
                C3 = (2.0f * (k1 - k2) + (k3 + k4) * Dka - C3) * Ry + C3;
            }
        }

        const double Ra = (alpha - ka[ia]) / Dka;

        // replace pow(Ra,2/3) with multiplies
        const double Ra2 = Ra * Ra;
        const double Ra3 = Ra2 * Ra;

        sadArray[ii] = (C0 + C1 * Ra + C2 * Ra2 + C3 * Ra3) * sadf;

        // extrapolation beyond last y knot
        if (yExtrapolation)
        {
            const double A1 = (double)ii;
            const double A2 = A - A1;

            const double A1_13 = cbrt(A1);

            const double r1 = A1_13 * r0;

            // NOTE: your original code had pow(A2,1/2) here (likely a bug),
            // but I keep it exactly as-is: A2^(1/2)
            const double r2 = sqrt(A2) * r0;

            const double Areduced = A1 * A2 / A;

            const double r1_2 = r1 * r1;
            const double r2_2 = r2 * r2;

            const double sum = (r1 + r2 + sep);
            const double sum2 = sum * sum;

            const double MomInertia =
                0.4f * A1 * r1_2 + 0.4f * A2 * r2_2 + Areduced * sum2;

            const double fJ2 = fJ * fJ;
            sadArray[ii] += fJ2 * (1.0f - y / yy) / (2.0f * MomInertia) * kRotate;
        }
    }

    for (int ii = 5; ii <= Narray; ii++)
    {
        const double A1 = (double)ii;
        const double A2 = A - A1;

        const double Z1 = (A1 / A) * Z;
        const double Z2 = Z - Z1;

        const int ia1 = (int)A1;
        const int ia2 = (int)A2;

        const int iz1 = (int)Z1;
        if (iz1 < 2) continue;

        const int iz2 = (int)Z2;

        double Mass1 = mass->getLDM(iz1, ia1);

        // turn off decays outside bounds
        if (Mass1 == -1000.0f)
        {
            sadArrayZA[ii] = 1000.0f;
            continue;
        }

        double tmp = mass->getLDM(iz1 + 1, ia1);
        Mass1 += (tmp - Mass1) * (Z1 - (double)iz1);

        double Mass2 = mass->getLDM(iz2, ia2);

        // turn off decays outside bounds
        if (Mass2 == -1000.0f)
        {
            sadArrayZA[ii] = 1000.0f;
            continue;
        }

        tmp = mass->getLDM(iz2 + 1, ia2);
        Mass2 += (tmp - Mass2) * (Z2 - (double)iz2);

        const double r1 = r0 * cbrt(A1);
        const double r2 = r0 * cbrt(A2);

        const double Ecoul = Z1 * Z2 * 1.44f / (r1 + r2 + sep);

        sadArrayZA[ii] = sadArray[ii] - Mass1 - Mass2 - Ecoul;
    }
}

//*****************************************************************************
/**
 * Prints out the array of conditional barriers
 */
void CYrast::printAsyBarrier()
{
    for (int i = 0; i < Narray; i++)
        cout << i << " " << sadArray[i] << endl;
}

//*****************************************************
/**
 * Returns the conditional saddle-point energy when both nascient fragments
 * have the same Z/A as the parent nucleus
 * \param A1 is the mass number of one of the nascient fragments
 */
double CYrast::getSaddlePointEnergy(double A1)
{
    if (A1 > (double)iA / 2.0f) A1 = (double)iA - A1;

    const int iA1 = (int)A1;

    const double Esaddle1 = sadArray[iA1];

    if (2 * iA1 + 1 >= iA) return Esaddle1;

    const double Esaddle2 = sadArray[iA1 + 1];

    return Esaddle1 + (Esaddle2 - Esaddle1) * (A1 - (double)iA1) + addBar;
}

//***********************************************
/**
 * Returns the constrained saddle-point energy for channel
 * where one of the fragments is iZ1,iA1
 * \param iZ1 is the proton number
 * \param iA1 is the mass number
 */
double CYrast::getSaddlePointEnergy(int iZ1, int iA1)
{
    const int iZ2 = iZ - iZ1;
    const int iA2 = iA - iA1;

    int iAmin = iA1;
    if (iA2 < iAmin) iAmin = iA2;

    const double mass1 = mass->getLDM(iZ1, iA1);
    const double mass2 = mass->getLDM(iZ2, iA2);

    const double r1 = r0 * cbrt((double)iA1);
    const double r2 = r0 * cbrt((double)iA2);

    const double Ecoul = 1.44f * (double)iZ1 * (double)iZ2 / (r1 + r2 + sep);

    return sadArrayZA[iAmin] + mass1 + mass2 + Ecoul + addBar;
}

//********************************************************
/**
 * General call to get Yrast energy according to Nuclear models
 * This handles problems if iA is too small or the angular momentum
 * is too large in extrapolating the Yrast line.
 * \param iZ is the proton number
 * \param iA is the mass number
 * \param fJ is the angular momentum
 */
double CYrast::getYrastModel(int iZ, int iA, double fJ)
{
    if (iZ < 19 || iZ > 102) return getYrastRLDM(iZ, iA, fJ);

    getJmaxSierk(iZ, iA);

    if (fJ < (double)Jmax - deltaJ) return getYrastSierk(fJ);

    double MInertia = getMomentOfInertiaSierk((double)Jmax - deltaJ);

    // Sierk routine has given negative inertias for very exotic nuclei
    if (MInertia <= 0.0f)
    {
        const double R = r0 * cbrt((double)iA);
        MInertia = 0.4f * (double)iA * (R * R);
    }

    const double fJ2 = fJ * fJ;
    const double J0  = (double)Jmax - deltaJ;
    const double J02 = J0 * J0;

    return getYrastSierk(fJ - deltaJ)
           + 0.5f * (fJ2 - J02) / MInertia * kRotate;
}

//********************************************************
/**
 * General call to get Yrast energy, includes a fix for light nuclei so that
 * the Yrast line does not become too steep. At a spin JChange, the Yrast line
 * continues with a constant slope.
 * \param iZ is the proton number
 * \param iA is the mass number
 * \param fJ is the angular momentum
 */
double CYrast::getYrast(int iZ, int iA, double fJ)
{
    if (bForceSierk) return getYrastModel(iZ, iA, fJ);

    const double fJChange = 0.319f * (double)iA;

    if (fJ < fJChange) return getYrastModel(iZ, iA, fJ);

    const double r1 = getYrastModel(iZ, iA, fJChange);
    const double r2 = getYrastModel(iZ, iA, fJChange + 1.0f);

    return r1 + (fJ - fJChange) * (r2 - r1);
}

//*******************************************************
/**
 * Returns the saddle-point energy for symmetric fission in MeV
 * \param iZ is the proton number
 * \param iA is the mass number
 * \param fJ is the angular momentum
 */
double CYrast::getSymmetricSaddleEnergy(int iZ, int iA, double fJ)
{
    if (iZ > 102)
        return getYrast(iZ, iA, fJ) + getBarrierFissionRLDM(iZ, iA, fJ);
    else
        return getYrast(iZ, iA, fJ) + getBarrierFissionSierk(fJ);
}

//**************************************************************
/**
 * Forces the use of Sierk yrast energies.
 */
void CYrast::forceSierk(bool b /*=1*/)
{
    bForceSierk = b;
}

//******************************************************************
/**
 * Prints out the parameters used in the Yrast class
 */
void CYrast::printParameters()
{
    cout << "forceSierk " << bForceSierk << endl;
}

//**********************************************
/**
 * calculates the Wigner energy in the mass formula used by Sierk
 * \param iZ is the proton number
 * \param iA is the mass number
 */
double CYrast::WignerEnergy(int iZ, int iA)
{
    const double azero = 2.693f;
    const double absI = fabs((double)(iA - 2 * iZ) / (double)iA);

    if (absI > 0.35f) return 6.716f + azero;

    // keep as-is; no pow here
    return 38.38f * absI * (1.0f - 0.5f * absI / 0.35f) + azero;
}
