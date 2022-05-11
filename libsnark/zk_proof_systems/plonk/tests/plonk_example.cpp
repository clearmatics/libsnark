#ifndef __PLONK_EXAMPLE_CPP__
#define __PLONK_EXAMPLE_CPP__

namespace libsnark
{

template<typename ppT>
void plonk_example<ppT>::initialize()
{
  using Field = libff::Fr<ppT>;
  using BaseField = libff::Fq<ppT>;

  // Hashes of transcript (Fiat-Shamir heuristic)
  this->beta = Field("3710899868510394644410941212967766116886736137326022751891187938298987182388");
  this->gamma = Field("11037930384083194587907709665332116843267274045828802249545114995763715746939");
  this->alpha = Field("37979978999274723893071781986484838492111162341880360022719385400306128734648");
  this->zeta = Field("43271972289218399355833643945502350270719103959803126415018065799136107272465");
  this->nu = Field("27515859833869775242150726508092341429478280783192379165155175653098691426347");
  this->u = Field("17817511439546966846324492112120565778288553881098836505706049265393896966778");

  // Prover Round 3
  this->t_poly_at_secret_g1 = 
    {
     {BaseField("3633475304544039580937168033891821181031270028948315156966430357290637750912918602224358395819959043217498580613188"), BaseField("1428090154951261810016759192966903623360639220861161704510358440208878251190328471919089961503194904379492282570328")}, // (x, y)   
     {BaseField("763634090322259543766607669979108502605520397912172619611323329140740033948682915660599655604319492439350037062593"), BaseField("2813678383705930006472398012708516812631766189864357429304341222779755096333176883586053913173384834727806732577514")}, // (x, y)   
     {BaseField("1133773332119974571006388114320487134122128432292374613610471191239740936855771046194807037399513728603857921779020"), BaseField("2371743385249340433047174208075481672774011018845240422821241326403735375534578397825283190840736410689009347296342")} // (x, y)   
    };
    // Prover Round 4
    this->a_zeta = Field("8901875463326310313456570652869873776746767429001289506712732994487869455294");
    this->b_zeta = Field("17059370482702287697833061796226204248201565415155528923232473993993212080397");
    this->c_zeta = Field("2409756965930008556696371654169913125397449986372522636184003898699439708220");
    this->S_0_zeta = Field("46143626155803287918330279428390848286076645428477353060129573054942492588828");
    this->S_1_zeta = Field("24392704635891252343143269633563768345145008520140360299402842967762646340846");
    this->t_zeta = Field("17704211079697158667898451781925539666888780633357685549668669638883218786797");
    this->z_poly_xomega_zeta = Field("28842520748921749858267479462666161290723351257502457358354355079408206613634");
    // Prover Round 5
    this->r_zeta = Field("9840183355354075764860129139740187852136872731945621853688663309524905254695");
    this->W_zeta =
      {
       Field("24252080758613024303887343640606716251261063588992752319177912308896977829530"), // x^0 +
       Field("10607607898724698632439123030466964271136356066368003670799965361976876241434"), // x^1 +
       Field("49455700908352041558678567416056380331270776757332153063055971634526190286005"), // x^2 +
       Field("33320271003589563855271608550901117652643277365972280669278308626702032462766"), // x^3 +
       Field("15852366479262249685982764572680658811329078253545174235941892300400403939930"), // x^4 +
       Field("27682430312398456895231605727369329939276517036904219760974903184897420916090"), // x^5 +
       Field("3408971485363305489368119320301696605838015578921377320829539984918748131464"), // x^6 +
       Field("6170013496654113309991805083384693399451071882671661442964650981706324964718"), // x^7 +
       Field("40466512078532355653582750156211155136663365180596769975433193406332496460952"), // x^8 +
       Field("27797551006624739724390744293569978212779205094709314274857289151379644029144") // x^9
      };
    this->W_zeta_omega =
      {
       Field("40531564636934977450521645157318380827315092224535055562348406035588261107649"), // x^0 +
       Field("40582772068296697445766687367246470550165702965822503698035189013689327598242"), // x^1 +
       Field("34692297336207296183471203952558746590281061299363350685811458233238031608744"), // x^2 +
       Field("36729769970326679258606619996857871496729667470839387898797395607660186144417"), // x^3 +
       Field("39200377238152535353153251658193113728678840254914604189132813789688481911843"), // x^4 +
       Field("37806771595110701245785306926880932563588614469577320079389456056927918503907"), // x^5 +
       Field("30089056814844115209371497598364690230296651200316391472693490173328399559224"), // x^6 +
       Field("50192235966668377039640655041398654654205952586673513120204602159086894893773"), // x^7 +
       Field("50464552726074655840942710135408889439064687734699646733572317469624082573673"), // x^8 +
       Field("2131891516651518828698089707163191982683101187444208330829153527689737950718") // x^9
      };
    this->W_zeta_at_secret =
      {
       BaseField("870060484392950318936407436339281693869104302033433437217724352106918879203859031480272360445121364037832041952039"), BaseField("861299420568581680683022801956519699850698057581115592449001068617883576936792838745784184017332399410405735556637")
      };
    this->W_zeta_omega_at_secret =
      {
       BaseField("655495718420813899433022178744291341055707361726913326276228363256669650441492160760400814686050651569669287015855"), BaseField("2294814434048221032840441892192852574031250430216011378262324574002174077316300644184560707652291246729455056384738")
      };
    // Verifier precomputation
    this->Q_polys_at_secret_g1 = 
      {
       // Q_polys_at_secret_G1[0](x, y)
       {BaseField("378726381176462718147300358135739414750401881865179496688045184488721729111812849671928978444183306485908524284480"), BaseField("791833214624823925392353480024899450306823972545233663950102448253555712343824442312037763251986444707559468789167")},
       // Q_polys_at_secret_G1[1](x, y)
       {BaseField("3463121064896078753387572937828548293123602506447011638964695206791719047814189335281941216929475385990343236354767"), BaseField("3756005527484874575873122885769433128395562570281048182601469211783316382910994832150136370603197748628978251617855")},
       // Q_polys_at_secret_G1[2](x, y)
       {BaseField("3892945593632124257345815822805968966792071575911821618411850268364131846825179296155052328167376852966904585186447"), BaseField("1438765187404699002937782524976196699749348173404933606055994266712582628388073773207681042085983472030360245506053")},
       // Q_polys_at_secret_G1[3](x, y)
       {BaseField("2139520368859124231883423692514970528859812258299930976367881008205428239812556102432321082924215659317599549354780"), BaseField("3160801670110892612311124955553991582667509305479627995819254406979123526107121808691993608981196593224342875940517")},
       // Q_polys_at_secret_G1[4](x, y)
       {BaseField("1494706197488050548473647189700590163018710083619049894201340042906428954701970245565963571442231938653450961032544"), BaseField("2824602300232981549589068758658165009614852853949575921795599675003541734658829364117937951492488210759581945412529")}
      };
    this->S_polys_at_secret_g1 = 
      {
       // S_polys_at_secret_G1[0](x, y)
       {BaseField("363680786286958401425193141820272230973918593225003485230375359988783103495068430872040401356145672343538149573956"), BaseField("3808611715527906761995774613532804514541116806047370905815923618967460517724494838454007652229237779324980791026678")},
       // S_polys_at_secret_G1[1](x, y)
       {BaseField("2269290281792079624043921770487152030658877587414533354808821851134353908789851247370093633774252241702882912887085"), BaseField("1202991994896142771673590759394828512458115068505041721692100348489247274098109255301291346677656642600399089495236")},
       // S_polys_at_secret_G1[2](x, y)
       {BaseField("1409279644592282847911075677223250748431662422519753258586136046888366320164223401617201095933133577003578716884726"), BaseField("1739252675334665913183857419645758811114836143883564006233172106501577251543276680000254892836167925142291524692240")}       
      };
    // Verifier Step 5: vanishing polynomial evaluation at zeta
    this->zh_zeta = Field("3231291541068425176363304614853753627837086743398580603659947674812058489341");
    // Verifier Step 6: compute Lagrange polynomial evaluation L1(zeta)
    this->L_0_zeta = Field("24902350288015759783891507297853223584346979652007592864411578966192407598153");
    // Verifier Step 7: evaluate public input polynomial at zeta
    this->PI_zeta = Field("52131913828646241587214622232821400526669239976784838364162926303506898992750");
    // Verifier Step 8: compute quotient polynomial evaluation r'(zeta) = r(zeta) - r0, where r0 is a constant term
    this->r_prime_zeta = Field("17704211079697158667898451781925539666888780633357685549668669638883218786797");
    // Verifier Step 9
    this->D1 =
      {
       // (x,y)
       BaseField("2950551289935145714650514276285309855008046336198238299122943701060700368614623093700969937886225712332445422779232"), 
       BaseField("1773574577433686795642287130007384519506152230425561199358320811214796879879223471447420612658036440896892253660321")
      };
      
}

} // namespace libsnark
  
#endif // __PLONK_EXAMPLE_CPP__
