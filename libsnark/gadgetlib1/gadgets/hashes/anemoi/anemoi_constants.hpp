/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_CONSTANTS_HPP_
#define LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_CONSTANTS_HPP_

namespace libsnark
{
// TODO: specialize by the field type + cast to the field
// see     setup_sha3_constants();
#if 0
// l = 1
FieldT C1[1][19] = {
    {39},
    {41362478282768062297187132445775312675360473883834860695283235286481594490621},
    {9548818195234740988996233204400874453525674173109474205108603996010297049928},
    {25365440569177822667580105183435418073995888230868180942004497015015045856900},
    {34023498397393406644117994167986720327178154686105264833093891093045919619309},
    {38816051319719761886041858113129205506758421478656182868737326994635468402951},
    {35167418087531820804128377095512663922179887277669504047069913414630376083753},
    {25885868839756469722325652387535232478219821850603640827385444642154834700231},
    {8867588811641202981080659274007552529205713737251862066053445622305818871963},
    {36439756010140137556111047750162544185710881404522379792044818039722752946048},
    {7788624504122357216765350546787885309160020166693449889975992574536033007374},
    {3134147137704626983201116226440762775442116005053282329971088789984415999550},
    {50252287380741824818995733304361249016282047978221591906573165442023106203143},
    {48434698978712278012409706205559577163572452744833134361195687109159129985373},
    {32960510617530186159512413633821386297955642598241661044178889571655571939473},
    {12850897859166761094422335671106280470381427571695744605265713866647560628356},
    {14578036872634298798382048587794204613583128573535557156943783762854124345644},
    {21588109842058901916690548710649523388049643745013696896704903154857389904594},
    {35731638686520516424752846654442973203189295883541072759390882351699754104989}};
FieldT D1[1][19] = {
    {14981678621464625851270783002338847382197300714436467949315331057125308909900},
    {28253420209785428420233456008091632509255652343634529984400816700490470131093},
    {51511939407083344002778208487678590135577660247075600880835916725469990319313},
    {46291121544435738125248657675097664742296276807186696922340332893747842754587},
    {3650460179273129580093806058710273018999560093475503119057680216309578390988},
    {45802223370746268123059159806400152299867771061127345631244786118574025749328},
    {11798621276624967315721748990709309216351696098813162382053396097866233042733},
    {42372918959432199162670834641599336326433006968669415662488070504036922966492},
    {52181371244193189669553521955614617990714056725501643636576377752669773323445},
    {23791984554824031672195249524658580601428376029501889159059009332107176394097},
    {33342520831620303764059548442834699069640109058400548818586964467754352720368},
    {16791548253207744974576845515705461794133799104808996134617754018912057476556},
    {11087343419860825311828133337767238110556416596687749174422888171911517001265},
    {11931207770538477937808955037363240956790374856666237106403111503668796872571},
    {3296943608590459582451043049934874894049468383833500962645016062634514172805},
    {7080580976521357573320018355401935489220216583936865937104131954142364033647},
    {25990144965911478244481527888046366474489820502460615136523859419965697796405},
    {33907313384235729375566529911940467295099705980234607934575786561097199483218},
    {25996950265608465541351207283024962044374873682152889814392533334239395044136}};

// l=2
FieldT C2[2][12] = {
    {39,
     17756515227822460609684409997111995494590448775258437999344446424780281143353},
    {41362478282768062297187132445775312675360473883834860695283235286481594490621,
     3384073892082712848969991795331397937188893616190315628722966662742467187281},
    {9548818195234740988996233204400874453525674173109474205108603996010297049928,
     51311880822158488881090781617710146800056386303122657365679608608648067582435},
    {25365440569177822667580105183435418073995888230868180942004497015015045856900,
     29347609441914902330741511702270026847909178228078752565372729158237774700914},
    {34023498397393406644117994167986720327178154686105264833093891093045919619309,
     2339620320400167830454536231899316133967303509954474267430948538955691907104},
    {38816051319719761886041858113129205506758421478656182868737326994635468402951,
     27338042530319738113354246208426108832239651080023276643867223794985578055610},
    {35167418087531820804128377095512663922179887277669504047069913414630376083753,
     42192983528513372869128514327443204912824559545179630597589572656156258515752},
    {25885868839756469722325652387535232478219821850603640827385444642154834700231,
     42721818980548514490325424436763032046927347769153393863616095871384405840432},
    {8867588811641202981080659274007552529205713737251862066053445622305818871963,
     23473499332437056484066006746048591864129988909190267521144125882222313735740},
    {36439756010140137556111047750162544185710881404522379792044818039722752946048,
     16497366583607480604161417644040292299204496829635795525393416854929276060989},
    {7788624504122357216765350546787885309160020166693449889975992574536033007374,
     16727395967350522643500778393489915391834352737211416857240725807058479128000},
    {3134147137704626983201116226440762775442116005053282329971088789984415999550,
     46525506418681456193255596516104416743523037046982280449529426136392814992763}};

FieldT D2[2][12] = {
    {14981678621464625851270783002338847382197300714436467949315331057125308909900,
     48720959343719104324739338388885839802998711550637402773896395605948383052052},
    {28253420209785428420233456008091632509255652343634529984400816700490470131093,
     6257781313532096835800460747082714697295034136932481743077166200794135826591},
    {51511939407083344002778208487678590135577660247075600880835916725469990319313,
     4386017178186728799761421274050927732938229436976005221436222062273391481632},
    {46291121544435738125248657675097664742296276807186696922340332893747842754587,
     13820180736478645172746469075181304604729976364812127548341524461074783412926},
    {3650460179273129580093806058710273018999560093475503119057680216309578390988,
     40385222771838099109662234020243831589690223478794847201235014486200724862134},
    {45802223370746268123059159806400152299867771061127345631244786118574025749328,
     50306980075778262214155693291132052551559962723436936231611301042966928400825},
    {11798621276624967315721748990709309216351696098813162382053396097866233042733,
     34806952212038537244506031612074847133207330427265785757809673463434908473570},
    {42372918959432199162670834641599336326433006968669415662488070504036922966492,
     22755759419530071315007011572076166983660942447634027701351681157370705921018},
    {52181371244193189669553521955614617990714056725501643636576377752669773323445,
     30334172084294870556875274308904688414158741457854908094300017436690480001547},
    {23791984554824031672195249524658580601428376029501889159059009332107176394097,
     19832360622723392584029764807971325641132953515557801717644226271356492507876},
    {33342520831620303764059548442834699069640109058400548818586964467754352720368,
     5828182614154296575131381170785760240834851189333374788484657124381010655319},
    {16791548253207744974576845515705461794133799104808996134617754018912057476556,
     23729797853490401568967730686618146850735129707152853256809050789424668284094}};

// l=3
FieldT C3[3][10] = {
    {39,
     17756515227822460609684409997111995494590448775258437999344446424780281143353,
     10188916128123599964772546147951904500865009616764646948187915341627970346879},
    {41362478282768062297187132445775312675360473883834860695283235286481594490621,
     3384073892082712848969991795331397937188893616190315628722966662742467187281,
     38536464596998108028197905645250196649287447208374169339784649587982292038621},
    {9548818195234740988996233204400874453525674173109474205108603996010297049928,
     51311880822158488881090781617710146800056386303122657365679608608648067582435,
     24596965950552905296088269899880882549715354660832391374009234980535928382152},
    {25365440569177822667580105183435418073995888230868180942004497015015045856900,
     29347609441914902330741511702270026847909178228078752565372729158237774700914,
     14356478667385969079309349540394948109414829921001045845599553435706989367858},
    {34023498397393406644117994167986720327178154686105264833093891093045919619309,
     2339620320400167830454536231899316133967303509954474267430948538955691907104,
     12136748919666286297989154404429099226154686992028401568133058190732008277996},
    {38816051319719761886041858113129205506758421478656182868737326994635468402951,
     27338042530319738113354246208426108832239651080023276643867223794985578055610,
     15580674179713644540398409523441814073810768449493940562136422009899312699155},
    {35167418087531820804128377095512663922179887277669504047069913414630376083753,
     42192983528513372869128514327443204912824559545179630597589572656156258515752,
     47389212411441573266379092392931599970417884729397156841216318364858334633325},
    {25885868839756469722325652387535232478219821850603640827385444642154834700231,
     42721818980548514490325424436763032046927347769153393863616095871384405840432,
     5855288403637341107158034195599277569854359593529752399086836976954392351035},
    {8867588811641202981080659274007552529205713737251862066053445622305818871963,
     23473499332437056484066006746048591864129988909190267521144125882222313735740,
     5696063807157149622355481994320806474692190935543821893362808351446578125354},
    {36439756010140137556111047750162544185710881404522379792044818039722752946048,
     16497366583607480604161417644040292299204496829635795525393416854929276060989,
     31479323495970113713816467604460499675889579912370034974841212556442942086146}};

FiledT D3[3][10] = {
    {14981678621464625851270783002338847382197300714436467949315331057125308909900,
     48720959343719104324739338388885839802998711550637402773896395605948383052052,
     11709610427641952476226704950218052763560489079301307464225164120801969364960},
    {28253420209785428420233456008091632509255652343634529984400816700490470131093,
     6257781313532096835800460747082714697295034136932481743077166200794135826591,
     11966422202069200811427605007493817363680804416274031195624148724039857787313},
    {51511939407083344002778208487678590135577660247075600880835916725469990319313,
     4386017178186728799761421274050927732938229436976005221436222062273391481632,
     663227665329044490605880474899933274574966982371072793854806732105730575244},
    {46291121544435738125248657675097664742296276807186696922340332893747842754587,
     13820180736478645172746469075181304604729976364812127548341524461074783412926,
     21821175320697611197161277831984495658213397245419754392657307036488476373765},
    {3650460179273129580093806058710273018999560093475503119057680216309578390988,
     40385222771838099109662234020243831589690223478794847201235014486200724862134,
     20738601554725926373596082603265918636164823648026470243422423735982938342408},
    {45802223370746268123059159806400152299867771061127345631244786118574025749328,
     50306980075778262214155693291132052551559962723436936231611301042966928400825,
     9105861908793877437599087016640061747418296780065295891365798855886560153752},
    {11798621276624967315721748990709309216351696098813162382053396097866233042733,
     34806952212038537244506031612074847133207330427265785757809673463434908473570,
     10559431278588446438155840088055546145087872298641007742921718770142881700525},
    {42372918959432199162670834641599336326433006968669415662488070504036922966492,
     22755759419530071315007011572076166983660942447634027701351681157370705921018,
     8881354201366797207686592249590682298565723459695719800911380560885170725516},
    {52181371244193189669553521955614617990714056725501643636576377752669773323445,
     30334172084294870556875274308904688414158741457854908094300017436690480001547,
     35548861917762862971011720475855172816698712671893796030607658203859222685056},
    {23791984554824031672195249524658580601428376029501889159059009332107176394097,
     19832360622723392584029764807971325641132953515557801717644226271356492507876,
     5370567718707734490084045178883836972105253285449736908577321570876055642415}};

// l=4
FieldT C4[4][10] = {
    {39,
     17756515227822460609684409997111995494590448775258437999344446424780281143353,
     10188916128123599964772546147951904500865009616764646948187915341627970346879,
     3814237141406755457246679946340702245820791055503616462386588886553626328449},
    {41362478282768062297187132445775312675360473883834860695283235286481594490621,
     3384073892082712848969991795331397937188893616190315628722966662742467187281,
     38536464596998108028197905645250196649287447208374169339784649587982292038621,
     37592197675289757358471908199906415982484124338112374453435292524131427342810},
    {9548818195234740988996233204400874453525674173109474205108603996010297049928,
     51311880822158488881090781617710146800056386303122657365679608608648067582435,
     24596965950552905296088269899880882549715354660832391374009234980535928382152,
     34036826250287807194659359129722586818079652442547178531030410684351456041117},
    {25365440569177822667580105183435418073995888230868180942004497015015045856900,
     29347609441914902330741511702270026847909178228078752565372729158237774700914,
     14356478667385969079309349540394948109414829921001045845599553435706989367858,
     9488013611624811735432450930006811652991761655550510302915118428283918068143},
    {34023498397393406644117994167986720327178154686105264833093891093045919619309,
     2339620320400167830454536231899316133967303509954474267430948538955691907104,
     12136748919666286297989154404429099226154686992028401568133058190732008277996,
     19442569822772655270268482835742480365499256802520510905846953360427433130058},
    {38816051319719761886041858113129205506758421478656182868737326994635468402951,
     27338042530319738113354246208426108832239651080023276643867223794985578055610,
     15580674179713644540398409523441814073810768449493940562136422009899312699155,
     4362660876979205605782410963041525734654031488177761934879852229226211686053},
    {35167418087531820804128377095512663922179887277669504047069913414630376083753,
     42192983528513372869128514327443204912824559545179630597589572656156258515752,
     47389212411441573266379092392931599970417884729397156841216318364858334633325,
     41487656259632727393098274178738763934249662924287956242704596746920012242443},
    {25885868839756469722325652387535232478219821850603640827385444642154834700231,
     42721818980548514490325424436763032046927347769153393863616095871384405840432,
     5855288403637341107158034195599277569854359593529752399086836976954392351035,
     18845851722124019325834426094831743068408557621685658713002749358354699910772},
    {8867588811641202981080659274007552529205713737251862066053445622305818871963,
     23473499332437056484066006746048591864129988909190267521144125882222313735740,
     5696063807157149622355481994320806474692190935543821893362808351446578125354,
     48558031599255072862103809681060565464555437399403822458902024251997890071747},
    {36439756010140137556111047750162544185710881404522379792044818039722752946048,
     16497366583607480604161417644040292299204496829635795525393416854929276060989,
     31479323495970113713816467604460499675889579912370034974841212556442942086146,
     52327065242455117582590188333899352706031813782154293138553490341266149456684}};

FieldT D4[4][10] = {
    {14981678621464625851270783002338847382197300714436467949315331057125308909900,
     48720959343719104324739338388885839802998711550637402773896395605948383052052,
     11709610427641952476226704950218052763560489079301307464225164120801969364960,
     3188799073106888901912065951229864304299742047220134499402570163601813730969},
    {28253420209785428420233456008091632509255652343634529984400816700490470131093,
     6257781313532096835800460747082714697295034136932481743077166200794135826591,
     11966422202069200811427605007493817363680804416274031195624148724039857787313,
     8876022912542631074912834764773050492660953075192093830253524158063181475941},
    {51511939407083344002778208487678590135577660247075600880835916725469990319313,
     4386017178186728799761421274050927732938229436976005221436222062273391481632,
     663227665329044490605880474899933274574966982371072793854806732105730575244,
     7956955597245727322388196907364651338722736293265717471854714933795446618648},
    {46291121544435738125248657675097664742296276807186696922340332893747842754587,
     13820180736478645172746469075181304604729976364812127548341524461074783412926,
     21821175320697611197161277831984495658213397245419754392657307036488476373765,
     14806577897118234786495606424219372997573800509149076370951604526939593458489},
    {3650460179273129580093806058710273018999560093475503119057680216309578390988,
     40385222771838099109662234020243831589690223478794847201235014486200724862134,
     20738601554725926373596082603265918636164823648026470243422423735982938342408,
     25898290090014076279086638237202313571292864987698437102115051403552551578909},
    {45802223370746268123059159806400152299867771061127345631244786118574025749328,
     50306980075778262214155693291132052551559962723436936231611301042966928400825,
     9105861908793877437599087016640061747418296780065295891365798855886560153752,
     48177591413367409915642056167048753041735583848456612607691620273026228709602},
    {11798621276624967315721748990709309216351696098813162382053396097866233042733,
     34806952212038537244506031612074847133207330427265785757809673463434908473570,
     10559431278588446438155840088055546145087872298641007742921718770142881700525,
     2511742758961381498086249076485723904703122022711664665388729650078747694082},
    {42372918959432199162670834641599336326433006968669415662488070504036922966492,
     22755759419530071315007011572076166983660942447634027701351681157370705921018,
     8881354201366797207686592249590682298565723459695719800911380560885170725516,
     19725785152035256359574211351446161592903393017031483635806025440159666669692},
    {52181371244193189669553521955614617990714056725501643636576377752669773323445,
     30334172084294870556875274308904688414158741457854908094300017436690480001547,
     35548861917762862971011720475855172816698712671893796030607658203859222685056,
     23828822166916376664523534857031979764654878164406016294521947902346141831375},
    {23791984554824031672195249524658580601428376029501889159059009332107176394097,
     19832360622723392584029764807971325641132953515557801717644226271356492507876,
     5370567718707734490084045178883836972105253285449736908577321570876055642415,
     24072177097374519292068993110945703798030958684413852593268331853573451397392}};
#endif
} // namespace libsnark

#endif // LIBSNARK_GADGETLIB1_GADGETS_HASHES_ANEMOI_CONSTANTS_HPP_
