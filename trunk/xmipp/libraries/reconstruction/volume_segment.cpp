/***************************************************************************
 *
 * Authors:    Carlos Oscar            coss@cnb.csic.es (2002)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include <data/args.h>
#include <data/morphology.h>
#include <data/filters.h>
#include "volume_segment.h"

// Read arguments ==========================================================
void Prog_segment_prm::read(int argc, char **argv)
{
    fn_vol = getParameter(argc, argv, "-i");
    voxel_mass = textToFloat(getParameter(argc, argv, "-voxel_mass", "-1"));
    dalton_mass = textToFloat(getParameter(argc, argv, "-dalton_mass", "-1"));
    aa_mass = textToFloat(getParameter(argc, argv, "-aa_mass", "-1"));
    sampling_rate = textToFloat(getParameter(argc, argv, "-sampling_rate", "-1"));
    fn_mask = getParameter(argc, argv, "-o", "");
    en_threshold = checkParameter(argc, argv, "-threshold");
    if (en_threshold)
        threshold = textToFloat(getParameter(argc, argv, "-threshold"));
    otsu=checkParameter(argc,argv,"-otsu");

    //Sjors
    wang_radius = textToInteger(getParameter(argc, argv, "-wang", "3."));
    do_prob = checkParameter(argc, argv, "-prob");

}

// Show ====================================================================
std::ostream & operator << (std::ostream &out, const Prog_segment_prm &prm)
{
    out << "Input file   : " << prm.fn_vol        << std::endl
        << "Voxel mass   : " << prm.voxel_mass    << std::endl
        << "Dalton mass  : " << prm.dalton_mass   << std::endl
        << "AA mass      : " << prm.aa_mass       << std::endl
        << "Sampling rate: " << prm.sampling_rate << std::endl
        << "Output mask  : " << prm.fn_mask       << std::endl
        << "Enable thres.: " << prm.en_threshold  << std::endl
        << "Threshold    : " << prm.threshold     << std::endl
        << "Otsu         : " << prm.otsu          << std::endl
        << "Wang radius  : " << prm.wang_radius   << std::endl
        << "Probabilistic: " << prm.do_prob       << std::endl
    ;
    return out;
}

// usage ===================================================================
void Prog_segment_prm::usage() const
{
    std::cerr << "   -i <input volume>       : Volume to segment\n"
    << "  [-voxel_mass  <mass>  |  : Mass in voxels\n"
    << "   [-dalton_mass <mass> |  : Mass in daltons\n"
    << "    -aa_mass     <mass>]   : Mass in aminoacids\n"
    << "   -sampling_rate <Tm>]    : Sampling rate (A/pix)\n"
    << "  [-o <output mask=\"\">]    : Output mask\n"
    << "  [-threshold <th=-1>]     : Thresholding method\n"
    << "  [-otsu]                  : Otsu's method segmentation\n"
    << "  [-wang <rad=3>]          : Radius [pix] for B.C. Wang cone\n"
    << "  [-prob]                  : Calculate probabilistic solvent mask\n"
    ;
}

// Produce side information ================================================
void Prog_segment_prm::produce_side_info()
{
    V.read(fn_vol);
    double sampling_rate3 = sampling_rate * sampling_rate * sampling_rate;
    if (voxel_mass == -1 && !en_threshold && !otsu)
    {
        if ((dalton_mass == -1 && aa_mass == -1) || sampling_rate == -1)
            REPORT_ERROR(1, "Prog_segment_prm: No way to compute voxel mass");
        if (dalton_mass != -1)
            voxel_mass = dalton_mass * 1.207 / sampling_rate3;
        else
            voxel_mass = aa_mass * 110 * 1.207 / sampling_rate3;
    }
    std::cout << std::endl << "Derived voxel_mass=" << voxel_mass << std::endl;
}

// Count voxels ============================================================
// Segment with a given threshold and compute the number of voxels in the
// biggest piece
//#define DEBUG
double segment_threshold(const Image<double> *V_in, Image<double> *V_out,
                         double threshold, bool do_prob)
{
    Image<double> aux;

    // Binarize input volume
    (*V_out)() = (*V_in)();
    (*V_out)().threshold("below", threshold, threshold);
    (*V_out)().binarize(threshold);

#ifdef DEBUG
    std::cout << threshold << std::endl;
    Image<double> save;
    save() = (*V_in)();
    save.write("PPP0.vol");
    save() = (*V_out)();
    save.write("PPP1.vol");
#endif

    if (!do_prob)
    {
        // Apply morphological opening to input volume
        aux().resize((*V_out)());
        opening3D((*V_out)(), aux(), 18, 0, 1);
        closing3D(aux(), (*V_out)(), 18, 0, 1);
#ifdef DEBUG
        save() = (*V_out)();
        save.write("PPP2.vol");
#endif

    }

    // Count the number of different objects
    int no_comp = label_image3D((*V_out)(), aux());
    Matrix1D<double> count(no_comp + 1);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(aux())
    count((int)aux(k, i, j))++;

#ifdef DEBUG
    std::cout << count << std::endl << std::endl;
    std::cout << "Press any key\n";
    char c;
    std::cin >> c;
#endif

    // Pick the maximum
    count(0) = 0; // We don't want to pick the background
    int imax;
    count.maxIndex(imax);

    // Select the mask with only that piece
    FOR_ALL_ELEMENTS_IN_ARRAY3D((*V_out)())
    (*V_out)(k, i, j) = aux(k, i, j) == imax;

    return count(imax);
}

void wang_smoothing(const Image<double> *V_in, Image<double> *V_out, int radius)
{

    int r2, radius2 = radius * radius;
    double sumw, weight;

    (*V_out)().resize((*V_in)());

    FOR_ALL_ELEMENTS_IN_ARRAY3D((*V_in)())
    {
        sumw = 0.;
        A3D_ELEM((*V_out)(), k, i, j) = 0.;
        for (int kp = k - radius; kp < k + radius; kp++)
        {
            if (kp > STARTINGZ((*V_in)()) && kp < FINISHINGZ((*V_in)()))
            {
                for (int ip = i - radius; ip < i + radius; ip++)
                {
                    if (ip > STARTINGY((*V_in)()) && ip < FINISHINGY((*V_in)()))
                    {
                        for (int jp = j - radius; jp < j + radius; jp++)
                        {
                            if (jp > STARTINGX((*V_in)()) && jp < FINISHINGX((*V_in)()))
                            {
                                r2 = (kp - k) * (kp - k) + (ip - i) * (ip - i) + (jp - j) * (jp - j);
                                if ((r2 < radius2) && (A3D_ELEM((*V_in)(), kp, ip, jp) > 0.))
                                {
                                    weight = 1. - sqrt((double)(r2 / radius2));
                                    A3D_ELEM((*V_out)(), k, i, j) += weight * XMIPP_MAX(0., A3D_ELEM((*V_in)(), kp, ip, jp));
                                    sumw += weight;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (sumw > 0.)
            A3D_ELEM((*V_out)(), k, i, j) /= sumw;
        else
            A3D_ELEM((*V_out)(), k, i, j) = 0.;
    }

}


void probabilistic_solvent(Image<double> *V_in, Image<double> *V_out)
{

    // Calculate mean and sigma for protein and solvent regions
    // according to the traditional segmentation
    double Np, sump, sum2p, Ns, sums, sum2s;
    double avgp, sigp, avgs, sigs, aux, solv_frac, prot_frac;
    double p_prot, p_solv;

    (*V_in)().setXmippOrigin();
    (*V_out)().setXmippOrigin();

    Np = sump = sum2p = Ns = sums = sum2s = 0.;
    FOR_ALL_ELEMENTS_IN_ARRAY3D((*V_in)())
    {
        aux = A3D_ELEM((*V_in)(), k, i, j);
        if (A3D_ELEM((*V_out)(), k, i, j) < 0.5)
        {
            sums += aux;
            sum2s += aux * aux;
            Ns += 1.;
        }
        else
        {
            sump += aux;
            sum2p += aux * aux;
            Np += 1.;
        }
    }
    if (Np > 0. && Ns > 0.)
    {
        avgs = sums / Ns;
        sigs = sum2s / Ns - avgs * avgs;
        avgp = sump / Np;
        sigp = sum2p / Np - avgp * avgp;
        prot_frac = Np / (Np + Ns);
        solv_frac = 1. - prot_frac;
    }
    else
    {
        REPORT_ERROR(1, "Prog_segment_prm: empty solvent or protein region");
    }

    // Terwilliger-like calculation of P(x|solv) & P(x|prot)
    // Bayes: P(prot|x)= P(x|prot)/{P(x|prot)+P(x|solv)}
    FOR_ALL_ELEMENTS_IN_ARRAY3D((*V_in)())
    {
        aux = A3D_ELEM((*V_in)(), k, i, j) - avgs;
        p_solv = solv_frac * exp(-aux * aux / (2 * sigs));
        aux = A3D_ELEM((*V_in)(), k, i, j) - avgp;
        p_prot = prot_frac * exp(-aux * aux / (2 * sigp));
        A3D_ELEM((*V_out)(), k, i, j) = p_prot / (p_prot + p_solv);
    }

}

// Really segment ==========================================================
void Prog_segment_prm::segment(Image<double> &mask)
{
    double th_min, th_max, val_min, val_max;
    V().computeDoubleMinMax(val_min, val_max);
    th_min = val_min;
    th_max = val_max;
    double mass_min = MULTIDIM_SIZE(V());
    double mass_max = 1;

    bool ok = false;
    if (!otsu)
    {
        if (!en_threshold)
        {
            // Perform a bracketing search until the mass is
            // within a 0.1% of the desired mass
            do
            {
                double th_med = (th_min + th_max) * 0.5;
                double mass_med = segment_threshold(&V, &mask, th_med, do_prob);
                std::cout << "Threshold= " << th_med
                << " mass of the main piece= " << mass_med << std::endl;
                if (ABS(mass_med - voxel_mass) / voxel_mass < 0.001)
                {
                    ok = true;
                    break;
                }
                if ((th_max - th_min) / (val_max - val_min) < 0.0001)
                    break;
                if (mass_med < voxel_mass)
                {
                    th_max = th_med;
                    mass_max = mass_med;
                }
                else
                {
                    th_min = th_med;
                    mass_min = mass_med;
                }
            }
            while (true);
        }
        else
        {
            // Perform a single thresholding
            double mass_med = segment_threshold(&V, &mask, threshold, do_prob);
            std::cout << "Threshold= " << threshold
            << " mass of the main piece= " << mass_med << std::endl;
            ok = true;
        }
    }

    if (otsu)
    {
        mask()=V();
        EntropyOtsuSegmentation(mask());
        ok=true;
    }

    if (do_prob)
    {
        // Wang-Leslie like modification of the input volume
        if (wang_radius >= 3)
        {
            Image<double> Vwang;
            wang_smoothing(&V, &Vwang, wang_radius);
            V = Vwang;
        }
        // Terwilliger-like calculation of P(solv|x) through P(x|solv) & P(x|prot)
        probabilistic_solvent(&V, &mask);
    }
    
    // Save mask if necessary
    if (fn_mask != "" && (ok || do_prob))
        mask.write(fn_mask);
    if (!ok && !do_prob)
        std::cout << "Segment: Cannot find an appropriate threshold\n";
}

