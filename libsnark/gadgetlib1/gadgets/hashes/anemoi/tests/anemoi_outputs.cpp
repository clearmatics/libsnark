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
                "16886130779841338454685815420547293349389636023659763220356242"
                "191375349870484"),
            libff::Fr<libff::bls12_381_pp>(
                "25110401710643410245000460073920576016937556620262171420761605"
                "822971640233467")};
    }

    // Expected output for 1 round, L=2: Y_left || Y_right
    if (NumStateColumns_L == 2) {
        Y_expect_one_round = {
            libff::Fr<libff::bls12_381_pp>(
                "1175380475393374013893416228835342241561414"
                "1018732075580217124721246447701179"),
            libff::Fr<libff::bls12_381_pp>(
                "5134979568569792732148078272800590308409837"
                "4830632868723399131618498422926132"),
            libff::Fr<libff::bls12_381_pp>(
                "3727673203363049293464225228142251472192664"
                "376552350730266112321515955286600"),
            libff::Fr<libff::bls12_381_pp>(
                "1961381747459264641165654006524664686879236"
                "9344942410316218713185877778202296")};
    }

    // Expected output for 1 round, L=3: Y_left || Y_right
    if (NumStateColumns_L == 3) {
        Y_expect_one_round = {
            libff::Fr<libff::bls12_381_pp>(
                "34028487665369117219732649464748473746959762451803899118620720"
                "068928912690818"),
            libff::Fr<libff::bls12_381_pp>(
                "23656628590950976492304898550690305252497146294682073762427630"
                "253215753757615"),
            libff::Fr<libff::bls12_381_pp>(
                "50489879152443944169393928755102146404491652050126328897483718"
                "597599270365689"),
            libff::Fr<libff::bls12_381_pp>(
                "48428584324500458296574300843220437759994411990569443572079173"
                "839484309978119"),
            libff::Fr<libff::bls12_381_pp>(
                "11634702739956834435596217222823133671340499123246628858964783"
                "053003807419954"),
            libff::Fr<libff::bls12_381_pp>(
                "21780864109569366263860246478259162683407782976390816039852496"
                "900891891832562")};
    }

    // Expected output for 1 round, L=4: Y_left || Y_right
    if (NumStateColumns_L == 4) {
        //        Y_expect_one_round = {0, 0, 0, 0, 0, 0, 0, 0};
        Y_expect_one_round = {
            libff::Fr<libff::bls12_381_pp>(
                "17041780326162393669606616053498620031347218348431899833333341"
                "812514930377576,"),
            libff::Fr<libff::bls12_381_pp>(
                "30554768545485035781675432957894191123666596590949762067679457"
                "455353618615899,"),
            libff::Fr<libff::bls12_381_pp>(
                "42087151902926409355779626374174416555426768306907549437379349"
                "373759646788220,"),
            libff::Fr<libff::bls12_381_pp>(
                "14414386843800103337001169291703636864672990050174868447277039"
                "577479472767841"),
            libff::Fr<libff::bls12_381_pp>(
                "31520834268440017158677087966007788520972054097794068312302569"
                "32773289540941,"),
            libff::Fr<libff::bls12_381_pp>(
                "38028057056256168614342506654290572930068443685959285901110523"
                "458386878293466,"),
            libff::Fr<libff::bls12_381_pp>(
                "29773101755555543082186573665458252067898132238258590143518897"
                "95128931607232,"),
            libff::Fr<libff::bls12_381_pp>(
                "51937274826670493627644353962426329420272866410811452792266413"
                "243735248502826")};
    }

    return Y_expect_one_round;
}

} // namespace libsnark
