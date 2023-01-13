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
    const size_t &NumStateColumns_L)
{
    std::vector<libff::Fr<libff::bls12_381_pp>> Y_expect_one_round;

    assert(
        ((NumStateColumns_L == 1) || (NumStateColumns_L == 2) ||
         (NumStateColumns_L == 3) || (NumStateColumns_L == 4)));

    // Expected output for 1 round, L=1: Y_left || Y_right
    if (NumStateColumns_L == 1) {
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
    if (NumStateColumns_L == 2) {
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
    if (NumStateColumns_L == 3) {
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
    if (NumStateColumns_L == 4) {
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

} // namespace libsnark
