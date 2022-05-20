#ifndef __PLONK_EXAMPLE_CPP__
#define __PLONK_EXAMPLE_CPP__

namespace libsnark
{

template<typename ppT>
//void plonk_example<ppT>::initialize()
plonk_example<ppT>::plonk_example()
{
  using Field = libff::Fr<ppT>;
  using BaseField = libff::Fq<ppT>;
  
  // Circuit data
  
  // number of gates / constraints. we have 6 gates for the example
  // circuit + 2 dummy gates to make it a power of 2 (for the fft)
  this->num_gates = 8;
    
  // number of q-polynomials
  this->num_qpolys = 5;
  
  // hard-coded gates matrix for the example circuit
  // P(x) = x**3 + x + 5 = 3
  // Each column is a q-vector
  this->gates_matrix =
    {
     // q_L     q_R        q_M         q_O       q_C
     {Field(0), Field(0), Field(1), -Field("1"), Field(0)},   // mul
     {Field(0), Field(0), Field(1), -Field("1"), Field(0)},   // mul
     {Field(1), Field(1), Field(0), -Field("1"), Field(0)},   // add
     {Field(0), Field(1), Field(0),  Field(0),  -Field("5")}, // con5
     {Field(0), Field(1), Field(0),  Field(0),   Field(0)},   // PI
     {Field(1), Field(1), Field(0), -Field("1"), Field(0)},   // add
     {Field(0), Field(0), Field(0),  Field(0),   Field(0)},   // dummy
     {Field(0), Field(0), Field(0),  Field(0),   Field(0)},   // dummy
    };
  
    this->gates_matrix_transpose = 
      {
       //  mul          mul          add          con5       PI           add        dum        dum
       {   Field(0),    Field(0),    Field(1),    Field(0),  Field(0),    Field(1),  Field(0),  Field(0)}, // q_L
       {   Field(0),    Field(0),    Field(1),    Field(1),  Field(1),    Field(1),  Field(0),  Field(0)}, // q_R
       {   Field(1),    Field(1),    Field(0),    Field(0),  Field(0),    Field(0),  Field(0),  Field(0)}, // q_M
       {-Field("1"), -Field("1"), -Field("1"),    Field(0),  Field(0), -Field("1"),  Field(0),  Field(0)}, // q_O
       {   Field(0),    Field(0),    Field(0), -Field("5"),  Field(0),    Field(0),  Field(0),  Field(0)}, // q_C
      };

    // witness values
    // w_L = a = [ 3,  9, 27,  1,  1, 30,  0,  0]
    // w_R = b = [ 3,  3,  3,  5, 35,  5,  0,  0]
    // w_O = c = [ 9, 27, 30,  5, 35, 35,  0,  0]
    // W = w_L + w_R + w_O
    this->witness = 
      {
       3,  9, 27,  1,  1, 30,  0,  0, // w_L 
       3,  3,  3,  5, 35,  5,  0,  0, // w_R
       9, 27, 30,  5, 35, 35,  0,  0  // w_O
      };
    
  // wire permutation (TODO: add function plonk_compute_permutation())
  this->wire_permutation =
    {9, 17, 18, 5, 4, 19, 7, 8, 10, 11, 1, 14, 21, 20, 15, 16, 2, 3, 6, 12, 22, 13, 23, 24};
  
  // public input (PI)
  this->public_input = Field(35);;
  // index of the row of the PI in the non-transposed gates_matrix 
  this->public_input_index = 4;
    
  // n-th root of unity omega in Fq (n=8 is the number of constraints
  // in the example). omega is a generator of the multiplicative
  // subgroup H.  Example (2**32)-th primitive root of unity in the
  // base field Fq of bls12-381 i.e. such that omega_base**(2**32) =
  // 1. The bls12-381 prime q is such that any power of 2 divides
  // (q-1). In particular 2**32|(q-1) and so the 2**32-th root of
  // unity exists.
  this->omega_base = Field("23674694431658770659612952115660802947967373701506253797663184111817857449850");

  // Constants k1,k2 to generate domains on which to evaluate the witness
  // polynomials. k can be random, but we fix it for debug to match
  // against the test vector values
  this->k1 = Field("7069874114745813936829552608791213902061117400356596714713673571023200548519"); 
  // Similarly, k2 can be random, but we fix it to match the test
  // vectors
  this->k2 = libff::power(k1, libff::bigint<1>(2));

    // H_gen contains the generators of H, k1 H and K2 H in one place
    // ie. omega, omega_k1 and omega_k2
    this->H_gen =
      {
       Field("1"),
       Field("23674694431658770659612952115660802947967373701506253797663184111817857449850"),
       Field("3465144826073652318776269530687742778270252468765361963008"),
       Field("8685283084174350996472453922654922162880456818468779543064782192722679779374"),
       Field("52435875175126190479447740508185965837690552500527637822603658699938581184512"),
       Field("28761180743467419819834788392525162889723178799021384024940474588120723734663"),
       Field("52435875175126190475982595682112313518914282969839895044333406231173219221505"),
       Field("43750592090951839482975286585531043674810095682058858279538876507215901405139"),
       Field("7069874114745813936829552608791213902061117400356596714713673571023200548519"),
       Field("629699347073785541163432187025543638744038878967631498033813889126570754409"),
       Field("21342932302184336920276494696368276078375881055107899047538799226396168459088"),
       Field("1267882965986156059121750719885063580789301664257248207931748295050157094180"),
       Field("45366001060380376542618187899394751935629435100171041107889985128915380635994"),
       Field("51806175828052404938284308321160422198946513621560006324569844810812010430104"),
       Field("31092942872941853559171245811817689759314671445419738775064859473542412725425"),
       Field("51167992209140034420325989788300902256901250836270389614671910404888424090333"),
       Field("18385921047832120226678773501713941136113926472106079933751452914767234401351"),
       Field("37876498801780702110808781477671785000889887091792406160116561883724775274403"),
       Field("20138953753137482450556344452362715831859795325354223322427824019686639989507"),
       Field("23457009893195201626861832175274146615046554789987966068815327283610815221310"),
       Field("34049954127294070252768967006472024701576626028421557888852205785171346783162"),
       Field("14559376373345488368638959030514180836800665408735231662487096816213805910110"),
       Field("32296921421988708028891396055823250005830757175173414500175834680251941195006"),
       Field("28978865281930988852585908332911819222643997710539671753788331416327765963203")
      };

    // H_gen permuted according to the wire permutation
    this->H_gen_permute =
      {
       Field("7069874114745813936829552608791213902061117400356596714713673571023200548519"),
       Field("18385921047832120226678773501713941136113926472106079933751452914767234401351"),
       Field("37876498801780702110808781477671785000889887091792406160116561883724775274403"),
       Field("52435875175126190479447740508185965837690552500527637822603658699938581184512"),
       Field("8685283084174350996472453922654922162880456818468779543064782192722679779374"),
       Field("20138953753137482450556344452362715831859795325354223322427824019686639989507"),
       Field("52435875175126190475982595682112313518914282969839895044333406231173219221505"),
       Field("43750592090951839482975286585531043674810095682058858279538876507215901405139"),
       Field("629699347073785541163432187025543638744038878967631498033813889126570754409"),
       Field("21342932302184336920276494696368276078375881055107899047538799226396168459088"),
       Field("1"),
       Field("51806175828052404938284308321160422198946513621560006324569844810812010430104"),
       Field("34049954127294070252768967006472024701576626028421557888852205785171346783162"),
       Field("23457009893195201626861832175274146615046554789987966068815327283610815221310"),
       Field("31092942872941853559171245811817689759314671445419738775064859473542412725425"),
       Field("51167992209140034420325989788300902256901250836270389614671910404888424090333"),
       Field("23674694431658770659612952115660802947967373701506253797663184111817857449850"),
       Field("3465144826073652318776269530687742778270252468765361963008"),
       Field("28761180743467419819834788392525162889723178799021384024940474588120723734663"),
       Field("1267882965986156059121750719885063580789301664257248207931748295050157094180"),
       Field("14559376373345488368638959030514180836800665408735231662487096816213805910110"),
       Field("45366001060380376542618187899394751935629435100171041107889985128915380635994"),
       Field("32296921421988708028891396055823250005830757175173414500175834680251941195006"),
       Field("28978865281930988852585908332911819222643997710539671753788331416327765963203")
      };

  // random hidden element secret (toxic waste). we fix it to a
  // constant in order to match against the test vectors
  this->secret = Field("13778279493383315901513166932749987230291710199728570152123261818328463629146");

  this->secret_powers_g1 =
    {
     {BaseField("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507"), BaseField("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569")},
     {BaseField("2074996566222649294561509797819818922173986132921308411358140896113079536003850694012630431408673580733705096617791"), BaseField("2077616938866677859533162977900874845236985734146041659617870201535410847474163664062006824876754586584291291275540")},
     {BaseField("2990806689640503496118106849118488768095304836122884937458915361027713460409591689021671430152609090653721640806775"), BaseField("1506910563119006480123164248508818566228973077756104948165605526387811452704182227509354555265743592891256699177447")},
     {BaseField("1162399422975338009465809685004642177085116952014234240402611534721365805950787073304245916047514747507357062412420"), BaseField("3243991324184282430666077430939015801339058147116419627536528847445203118923037149866561100103780253903700544760551")},
     {BaseField("665955504776459554223176911875912280457926009416658631176550711534449750755438228569302251439190379809780541605589"), BaseField("1563975451886550958452054831990339170002818201816232376486289252647203460244570019694049681767056597843416540080200")},
     {BaseField("1399715533473894434645348443452040407038595286119213991524009101772100073848392561767527055222979107182734281429257"), BaseField("2480196718611802082753242427207067840895581757500620454655673428668098605801043836809939785138234969928618554184939")},
     {BaseField("721989966033562332082331424380238800020350342699499474810940866633976625729399263772964415588210693897376315032210"), BaseField("2390971308801833445511417997165267564276479228294981140991537842761274654400117722746265926005170825975983654738210")},
     {BaseField("2494645203339367172765688786928769556042933531286628315780459840863253114375261836717471927002510598950105128289874"), BaseField("223089474026824320876534240998641224932692629684941862228540161286780230155004903395236049749437555611923506858957")},
     {BaseField("1369209115618410018904052993028425029992413842182015145017724608932835902017735132442291326610507957006806305596987"), BaseField("1669799775261473018291838448563339237460560373867009551270385685110329728891027501463328761827882531243313994834792")},
     {BaseField("3232461715720141976784119984008916789236798226172259100991524828275951110441373038297570085682311463031127657083499"), BaseField("2174858594096909044747114643824136389829031071589294112891897914222472713836544371512541030476560224224399906695672")},
     {BaseField("3294519054154733393573452083589753739284068002251231468836206591341493107913789146528048543749767535356614865447074"), BaseField("1078423904143836995685260825957768951259396575353368887368751656219548634950458729521798081562106307219496322329732")}
    };
  
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
    // Verifier Step 10: compute full batched polynomial commitment
    this->F1 =
      {
       // (x,y)
       BaseField("1826742492115868612582288009547441561233030549024396566895587306861582709499257146203756040438377200227370973259986"),
       BaseField("3174698003204604305970390381492147393289612967828928986277720440012247826580102584515518061223894220403192017952478")
      };
    // Verifier Step 11: compute group-encoded batch evaluation [E]_1
    this->E1 =
      {
       // (x,y)
       BaseField("896211057343888615494168233046980277322670342780717794143414541657336321268995118530971981505762603044924183423639"),
       BaseField("2430669190228433775014477375800217003654423165984629086408919423216148134888225770566675462506053600277852687652642")
      };
    // Verifier Step 12: batch validate all evaluations via pairing
    this->pairing_first_lhs =
      {
       // (x,y)
       BaseField("2594415027929961255375100283893836484494149195677830805487527382596454063206138678026974672787696177127251762793677"),
       BaseField("2352509828946078049326116392855469074202342093596555171015421209941442208769351281087086741535791209871063622971343")
      };
    this->pairing_first_rhs =
      {
       // (x,y)
#if 0 // pairing check e() = e()
       BaseField("2786086768974850272772873952573647917824800608941506526890251802946860616723974417792544909526285497672211326145554"),
       BaseField("3322299130191197851913727020573383402853736418546119609247136782819928908204641650829975455690084120791575793280234")
#endif  // #if 0 // pairing check e() = e()
       // (x,y)^-1
#if 1 // pairing check e() * e()^-1 = 1
       BaseField("2786086768974850272772873952573647917824800608941506526890251802946860616723974417792544909526285497672211326145554"), 
       BaseField("680110425030469541504062805162520753703146401392888276084921353304102742286196213612712173438931543246318479279553")
#endif // #if 1 // pairing check e() * e()^-1 = 1       
      };
}

} // namespace libsnark
  
#endif // __PLONK_EXAMPLE_CPP__
