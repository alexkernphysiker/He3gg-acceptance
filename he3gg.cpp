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
    	Plotter::Instance().SetOutput(".","sim-2g");
	const RandomValueTableDistr<>pf3=ReadPfFromFile("distributions/he3eta-pf-90-20.txt");
        const RandomUniform<> Pb_distr(p_beam_low,p_beam_hi); 
        const size_t ev_count=100000;
        size_t count=0;
        for(size_t event=0;event<ev_count;event++){
            const auto C=Compound(Pb_distr,pf3);
            const auto gammas_cme=binaryDecay(C.eta_.M(),0.0,0.0,randomIsotropic<3>());
            const auto g1Plab=gammas_cme.first.Transform(-C.eta_.Beta());
            const auto g2Plab=gammas_cme.second.Transform(-C.eta_.Beta());
            const auto he3dir=direction(C.he3.P());
            
            if((he3dir.th()>(3.0*PI()/180.))&&(he3dir.th()<(18.0*PI()/180.))){
                const auto g1dir=direction(g1Plab.P());
                const auto g2dir=direction(g2Plab.P());
                if(
                    (g1dir.th()>(20.*PI()/180.))&&(g1dir.th()<(169.*PI()/180.))&&
                    (g2dir.th()>(20.*PI()/180.))&&(g2dir.th()<(169.*PI()/180.))
                )
                    count++;
            }
        }
        const auto acceptance=std_error(double(count))/ev_count;
        cout<<"Geometry acceptance is "<<acceptance<<endl;
}
