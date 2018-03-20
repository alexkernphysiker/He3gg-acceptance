#include <iostream>
#include <math_h/lorentzvector.h>
#include <gnuplot_wrap.h>
#include <Experiment/experiment_conv.h>
#include <Kinematics/particles.h>
#include "common.h"
using namespace std;
using namespace MathTemplates;
using namespace GnuplotWrap;
int main()
{
    	Plotter::Instance().SetOutput(".","simulation-geom");
	const RandomValueTableDistr<>pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
        const RandomUniform<> Pb_distr(p_beam_low,p_beam_hi); 
        const auto M_thr=Particle::he3().mass()+Particle::eta().mass();
        const size_t ev_count=100000;
        {//3He 2gamma
            Distribution1D<> registered_count(BinsByStep(-70.,2.5,30.));
            Distribution1D<> all_count(BinsByStep(-70.,2.5,30.));
            for(size_t event=0;event<ev_count;event++){
                const auto C=Compound(Pb_distr,pf3);
                const auto Q=((C.he3+C.eta_).M()-M_thr)*1000.;
                const auto gammas_cme=binaryDecay(C.eta_.M(),0.0,0.0,randomIsotropic<3>());
                const auto g1Plab=gammas_cme.first.Transform(-C.eta_.Beta());
                const auto g2Plab=gammas_cme.second.Transform(-C.eta_.Beta());
                const auto he3dir=direction(C.he3.P());
                all_count.Fill(Q);
                if((he3dir.th()>(3.0*PI()/180.))&&(he3dir.th()<(18.0*PI()/180.))){
                    const auto g1dir=direction(g1Plab.P());
                    const auto g2dir=direction(g2Plab.P());
                    if(
                        (g1dir.th()>(20.*PI()/180.))&&(g1dir.th()<(169.*PI()/180.))&&
                        (g2dir.th()>(20.*PI()/180.))&&(g2dir.th()<(169.*PI()/180.))
                    )
                        registered_count.Fill(Q);
                }
            }
            Plot().Hist(registered_count*100./all_count)<<"set yrange [0:]"
            <<"set ylabel 'Acceptance, percents'"<<"set xlabel 'Q, MeV'"
            <<"set title 'pd->3He 2gamma'";
        }
        {//3He 6gamma
            Distribution1D<> registered_count(BinsByStep(-70.,2.5,30.));
            Distribution1D<> all_count(BinsByStep(-70.,2.5,30.));
            for(size_t event=0;event<ev_count;event++){
                const auto C=Compound(Pb_distr,pf3,pow(3.0*Particle::pi0().mass(),2));//eta must have enough energy for 3 pi0's
                const auto Q=((C.he3+C.eta_).M()-M_thr)*1000.;
                auto gammas=ThreePi0Decay(C.eta_);
                const auto he3dir=direction(C.he3.P());
                all_count.Fill(Q);
                if((he3dir.th()>(3.0*PI()/180.))&&(he3dir.th()<(18.0*PI()/180.))){
                    bool accepted=true;
                    for(const auto&gamma:gammas){
                        const auto gdir=direction(gamma.P());
                        accepted&=((gdir.th()>(20.*PI()/180.))&&(gdir.th()<(169.*PI()/180.)));
                    }
                    if(accepted)registered_count.Fill(Q);
                }
            }
            Plot().Hist(registered_count*100./all_count)<<"set yrange [0:]"
            <<"set ylabel 'Acceptance, percents'"<<"set xlabel 'Q, MeV'"
            <<"set title 'pd->3He 6gamma'";
        }
}
