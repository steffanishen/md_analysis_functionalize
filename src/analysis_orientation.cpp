/*
    read_dcd : c++ class + main file example for reading a CHARMM dcd file
    Copyright (C) 2013  Florent Hedin
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//***************** Partially contributed by Meng Shen since 06/10/2018 ***************************
//***************** Email: steffanishen@gmail.com ***************************************


#include <cstdlib>
#include <stdio.h>
#include <limits>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include <math.h>
#include "analysis_orientation.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_ORIENTATION::ANALYSIS_ORIENTATION(PSF *system, GROUP *sel1, int vector1d, int vector2d, int voidf, string filename, string name0, string name1, string name2, string name3, int nbins)
{
    this->system = system;
    this->sel1 = sel1;
    this->sel2 = sel2;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->name0 = name0;
    this->name1 = name1;
    this->name2 = name2;
    this->name3 = name3;
    this->dist_crit = dist_crit;
    this->nbins = nbins;
    this->rdf_count.resize(nbins);
    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    this->iframe = 0;


}

void ANALYSIS_ORIENTATION::init() {
}

vector<float> ANALYSIS_ORIENTATION::compute_vector() {
    sel1->anglezs.clear();
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> r2(3,0.0);
    vector<float> r3(3,0.0);
    vector<float> disp1(3,0.0);
    vector<float> disp2(3,0.0);
    float dist2;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float dist_crit2 = this->dist_crit*this->dist_crit;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    this->iframe += 1;
    float costheta = 0.0;
    vector<double> box;
    box.resize(3);
    box[2] = system->box_first_frame[2];
    float shift = box[2] * 0.5;
    float dz = box[2]/float(this->nbins);

    //cout <<"dist_crit: " << this->dist_crit<< "  dr: " << dr << " rdf bin: " << nbins << endl;
    //cout << "rdf nbins: " << nbins << endl;

    this->rdf_count.clear();
    this->rdf_count.resize(nbins);
    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);

    vector<float> costheta2s(nbins,0.0);
    vector<float> order_parameters(nbins,0.0);

    if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    //cout << "sel1->NATOM: " << sel1->NATOM << endl; //for debug purpose
    //cout << "sel2->NATOM: " << sel2->NATOM << endl; //for debug purpose
//    cout << "M_PI" << M_PI << endl;// for debug purpose
 
    for (auto &residue:sel1->residues) {
	    for (int ind : residue) {
            if (system->atomname[ind] == name0) {
                r[0] = system->x[ind];
	            r[1] = system->y[ind];
	            r[2] = system->z[ind];
                //cout << system->atomname[ind] << endl;
            }

            if (system->atomname[ind] == name1) {
                r1[0] = system->x[ind];
	            r1[1] = system->y[ind];
	            r1[2] = system->z[ind];
             //   cout << system->atomname[ind] << endl;
            }
	    
            if (system->atomname[ind] == name2) {
                r2[0] = system->x[ind];
	            r2[1] = system->y[ind];
	            r2[2] = system->z[ind];
             //   cout << system->atomname[ind] << endl;
            }

            if (system->atomname[ind] == name3) {
                r3[0] = system->x[ind];
	            r3[1] = system->y[ind];
	            r3[2] = system->z[ind];
             //   cout << system->atomname[ind] << endl;
            }
 
	    }
        disp1 = getDistPoints(r1,r2);
        disp2 = getDistPoints(r3,r2);

        vector<float> plane_norm = cross_product(disp1,disp2);
        vector<float> nz = {0.0, 0.0, 1.0};
        costheta = dot_product(plane_norm,nz)/(norm(plane_norm)*norm(nz));

        float costheta2 = costheta * costheta; 

	//cout << "costheta2: " << costheta2 << endl;

        int izbin = int((r[2] + shift)/dz);
        if (izbin >=0 and izbin < this->nbins) {
            this->rdf_count[izbin] += 1.0;
            costheta2s[izbin] += costheta2;
        }
    }


    for (int izbin = 0; izbin <= nbins; izbin++) {
        if (this->rdf_count[izbin] > 0.000001) {
	    costheta2s[izbin] = costheta2s[izbin] / this->rdf_count[izbin];
            order_parameters[izbin] = (3.0 * costheta2s[izbin] - 1.0) * 0.5;
//            order_parameters[izbin] = (3.0 * costheta2s[izbin] - 1.0 ) * 0.5;
            //order_parameters[izbin] = this->rdf_count[izbin];
        }
 //       cout <<"rdf: " << rdf[ibin] << endl;
    }

/*
*/

    return order_parameters;
}


ANALYSIS_ORIENTATION::~ANALYSIS_ORIENTATION()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
    this->rdf_count.clear();
}

