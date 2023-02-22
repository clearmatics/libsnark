/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by Clearmatics Ltd
 *             (originally developed by SCIPR Lab) and contributors
 *             (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include "libsnark/gadgetlib1/gadgets/hashes/anemoi/tests/anemoi_outputs.hpp"

namespace libsnark
{

std::vector<libff::Fr<libff::bls12_381_pp>> anemoi_expected_output_one_round(
    const size_t &NumStateColumns)
{
    std::vector<libff::Fr<libff::bls12_381_pp>> Y_expect_one_round;

    assert(
        ((NumStateColumns == 1) || (NumStateColumns == 2) ||
         (NumStateColumns == 3) || (NumStateColumns == 4)));

    // Expected output for 1 round, L=1: Y_left || Y_right
    if (NumStateColumns == 1) {
        Y_expect_one_round = {
            libff::Fr<libff::bls12_381_pp>(
                "38051718563229095456356396838757622428877349000988080406936"
                "541035058348383373"),
            libff::Fr<libff::bls12_381_pp>(
                "46541011259834287958249207092806566220478802569831738513953"
                "284817094880352697"),
        };
    }

    // Expected output for 1 round, L=2: Y_left || Y_right
    if (NumStateColumns == 2) {
        Y_expect_one_round = {
            libff::Fr<libff::bls12_381_pp>(
                "15150541060175709103777475248496599766370694616692747879011"
                "019662924685442224"),
            libff::Fr<libff::bls12_381_pp>(
                "29843552910061109491271060352447525363732495152033190876686"
                "584609246472017584"),
            libff::Fr<libff::bls12_381_pp>(
                "26146505138638275289195845765973260067149980640992713539811"
                "89299166526395914"),
            libff::Fr<libff::bls12_381_pp>(
                "49824839783019326099876978052724035902783619514814699665266"
                "556433444232424513"),
        };
    }

    // Expected output for 1 round, L=3: Y_left || Y_right
    if (NumStateColumns == 3) {
        Y_expect_one_round = {
            libff::Fr<libff::bls12_381_pp>(
                "10213223669833360114287009308428395240580814943870872556412"
                "118775684298316596"),
            libff::Fr<libff::bls12_381_pp>(
                "21664318220192052598342324987452326886678438558734363929475"
                "644562220502833005"),
            libff::Fr<libff::bls12_381_pp>(
                "12646567985368940694364168913172674258570854544904909064049"
                "69831326941799205"),
            libff::Fr<libff::bls12_381_pp>(
                "27292794672043705408598844612721784937283668235394047343788"
                "755455342406044808"),
            libff::Fr<libff::bls12_381_pp>(
                "38119908930143426720630804902252966609368611078523098634130"
                "250397080737556763"),
            libff::Fr<libff::bls12_381_pp>(
                "33144463221517343347312859079453261424067069247167408451667"
                "76226968891742442"),
        };
    }

    // Expected output for 1 round, L=4: Y_left || Y_right
    if (NumStateColumns == 4) {
        Y_expect_one_round = {
            libff::Fr<libff::bls12_381_pp>(
                "32728029339990442022355611963591129142873176406157617761037"
                "682452539219819088"),
            libff::Fr<libff::bls12_381_pp>(
                "49311030695492657479127064149718592806428520311648889321916"
                "443844338668678323"),
            libff::Fr<libff::bls12_381_pp>(
                "32674068623897120493809932137805441331335518262582687356675"
                "282427934524464097"),
            libff::Fr<libff::bls12_381_pp>(
                "47598191392555380432433599763649181065734671583813678477762"
                "106446112059712955"),
            libff::Fr<libff::bls12_381_pp>(
                "33816823933838499380471531208184598784697974497023773863342"
                "709156643799272839"),
            libff::Fr<libff::bls12_381_pp>(
                "39287313504077148049871180527977732866247140639315642469731"
                "817272046249721803"),
            libff::Fr<libff::bls12_381_pp>(
                "11801804404965184856947462551550819967745058632589947991520"
                "598952500171782048"),
            libff::Fr<libff::bls12_381_pp>(
                "13464045906969562610953775876862024564309085018098935959036"
                "17609367350545744"),
        };
    }

    return Y_expect_one_round;
}

// Output values automatically generated with SAGE script parameters.sage on
// 22/2/2023 at 15:42:52

std::vector<libff::Fr<libff::bls12_381_pp>> anemoi_expected_output_sec128(
    const size_t &NumStateColumns)
{
    std::vector<libff::Fr<libff::bls12_381_pp>> Y_expect;

    assert(
        ((NumStateColumns == 1) || (NumStateColumns == 2) ||
         (NumStateColumns == 3) || (NumStateColumns == 4)));

    // Expected output for X rounds, L=1: Y_left || Y_right
    if (NumStateColumns == 1) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_1_COL_128_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "19463313543534248726432829720355949992115481142527733903500993"
                "416436681359462"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "24097836352777748145579459445766630401317354855690625231253784"
                "292194604293571"),
        };
    }

    // Expected output for X rounds, L=2: Y_left || Y_right
    if (NumStateColumns == 2) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_2_COL_128_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "20636196687232276298438196241398690677938366239053904563010834"
                "541309527854991"),
            libff::Fr<libff::bls12_381_pp>(
                "39253916165022249533082065477479430579750918696028787690371014"
                "52687525354031"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "47061893317231315260263422053135889963204648975511761903522774"
                "86385893809395"),
            libff::Fr<libff::bls12_381_pp>(
                "45793334370241758458544071272207820701498478315901461355994185"
                "085013035430493"),
        };
    }

    // Expected output for X rounds, L=3: Y_left || Y_right
    if (NumStateColumns == 3) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_3_COL_128_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "28715769826978654676588546137117600748550476950408785857193565"
                "153875705420975"),
            libff::Fr<libff::bls12_381_pp>(
                "30425229241496163422803804141580461486929064146611232216062976"
                "254737479882403"),
            libff::Fr<libff::bls12_381_pp>(
                "22789248741438984532214659536852322247722236081022972364883893"
                "282168923848450"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "83340274330612138707941626621488069857370856538793932637490252"
                "41344734749219"),
            libff::Fr<libff::bls12_381_pp>(
                "38341449265920577314536214955239005955863226132621074752657782"
                "328514346130364"),
            libff::Fr<libff::bls12_381_pp>(
                "28073280982375140090919520522401329658310313884355685437076001"
                "811914604043364"),
        };
    }

    // Expected output for X rounds, L=4: Y_left || Y_right
    if (NumStateColumns == 4) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_4_COL_128_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "19356149758649626294450657712280508522363341709284507844575489"
                "610192191514962"),
            libff::Fr<libff::bls12_381_pp>(
                "43155080790752903177252471274668523810674926564664489782972593"
                "085074051296840"),
            libff::Fr<libff::bls12_381_pp>(
                "41787658699965881793824486198005321620242003410745831867294329"
                "180279986418075"),
            libff::Fr<libff::bls12_381_pp>(
                "45207360550431743775255889352626333429066098966478431547522957"
                "807501170661928"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "37717615359665988354289045467804278397912413188575884030395250"
                "092911699801716"),
            libff::Fr<libff::bls12_381_pp>(
                "36147903245738985095503661861599596609184057062914456862066826"
                "887006062632577"),
            libff::Fr<libff::bls12_381_pp>(
                "41936718433766085193462460233543805461189771221168643535359299"
                "387832399336324"),
            libff::Fr<libff::bls12_381_pp>(
                "35115766378970860468260290086983943422955413756722324468203198"
                "670156540279242"),
        };
    }

    return Y_expect;
}

std::vector<libff::Fr<libff::bls12_381_pp>> anemoi_expected_output_sec256(
    const size_t &NumStateColumns)
{
    std::vector<libff::Fr<libff::bls12_381_pp>> Y_expect;

    assert(
        ((NumStateColumns == 1) || (NumStateColumns == 2) ||
         (NumStateColumns == 3) || (NumStateColumns == 4)));

    // Expected output for X rounds, L=1: Y_left || Y_right
    if (NumStateColumns == 1) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_1_COL_256_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "23592301003500664995929917733266743854375415916040326652210568"
                "078795386647097"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "46415751505026927687934363401060779170128061056083999362978931"
                "09307338071265"),
        };
    }

    // Expected output for X rounds, L=2: Y_left || Y_right
    if (NumStateColumns == 2) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_2_COL_256_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "14658224666697525543665581557177920897018725487885393770688961"
                "291142066565838"),
            libff::Fr<libff::bls12_381_pp>(
                "22696721298764611334499842832788494984317317977616725617062055"
                "973718288462083"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "41254164318270696844698147640474784924833116221680739817900833"
                "75015002260367"),
            libff::Fr<libff::bls12_381_pp>(
                "85516268014959707663369712918470356171127697856604948332123300"
                "88434570655220"),
        };
    }

    // Expected output for X rounds, L=3: Y_left || Y_right
    if (NumStateColumns == 3) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_3_COL_256_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "39354654009222930922076215480960296260844083612260032711385075"
                "821419598992495"),
            libff::Fr<libff::bls12_381_pp>(
                "35857978005837223530941291838645726295300924278186042951565477"
                "753001911648704"),
            libff::Fr<libff::bls12_381_pp>(
                "73133080760004009022948041840133483135280922655558522590705858"
                "39828489100756"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "27964702120041308262791208736057137783851503845042628202338953"
                "806313118894432"),
            libff::Fr<libff::bls12_381_pp>(
                "32605073953207604732414082895795832873262647917998760850707996"
                "197372829415078"),
            libff::Fr<libff::bls12_381_pp>(
                "26789211140092257781038997601682984740469870292740801416378551"
                "890114989182776"),
        };
    }

    // Expected output for X rounds, L=4: Y_left || Y_right
    if (NumStateColumns == 4) {
        Y_expect = {
            // A_BLS_12_381_SCALARFIELD_4_COL_256_BITS
            // Left outputs
            libff::Fr<libff::bls12_381_pp>(
                "37935275437371366775666670594948501766272500250774338143461366"
                "789714115968944"),
            libff::Fr<libff::bls12_381_pp>(
                "46013056163891430736716124133522777555340005169288884686271429"
                "476146703061386"),
            libff::Fr<libff::bls12_381_pp>(
                "34292192884279029656147093024511389074726243647468508480039887"
                "915986427463842"),
            libff::Fr<libff::bls12_381_pp>(
                "19361591618166593449330002683801850115916835743219702190183265"
                "100898653128350"),
            // Right outputs
            libff::Fr<libff::bls12_381_pp>(
                "50276857710105470204512033208185084240001478327548370600191258"
                "668539356034348"),
            libff::Fr<libff::bls12_381_pp>(
                "30612042191114496995886018278834216278786456996263573491525930"
                "890999097880377"),
            libff::Fr<libff::bls12_381_pp>(
                "44865148646864607507025476306280271506300928383843934665366106"
                "390302816945546"),
            libff::Fr<libff::bls12_381_pp>(
                "12299595100938587761671298225363168979705772605172431805112941"
                "997639442458604"),
        };
    }

    return Y_expect;
}

} // namespace libsnark
