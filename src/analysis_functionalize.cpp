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
#include "analysis_functionalize.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_FUNCTIONALIZE::ANALYSIS_FUNCTIONALIZE(PSF *system, GROUP *sel1, GROUP *sel2, int vector1d, int vector2d, int voidf, string filename, float dist_crit, float dr)
{
    this->system = system;
    this->sel1 = sel1;
    this->sel2 = sel2;
    this->vector1d = vector1d;
    this->vector2d = vector2d;
    this->voidf = voidf;
    this->filename = filename;
    this->dist_crit = dist_crit;
    this->dr = dr;
    nbins = int(this->dist_crit/dr);
    this->rdf_count.resize(nbins);
    fill(this->rdf_count.begin(), this->rdf_count.end(),0.0);
    this->iframe = 0;
}

void ANALYSIS_FUNCTIONALIZE::init() {
}

vector<vector<vector<int>>> ANALYSIS_FUNCTIONALIZE::head_cell(vector<vector<int>> segments_ind) {
    vector<vector<vector<int>>> head;
    head.resize(xcount,vector<vector<int>>(ycount,vector<int>(zcount,-1)));//,vector<int>(zcount)));

    int ixcell,iycell,izcell;
     for (auto &segment:segments_ind) {
	    for (int ind : segment) {
            ixcell = int ((system->x[ind] + system->pbc[0]*0.5)/cellsize);
            if (ixcell < 0) ixcell = 0;
            else if (ixcell > xcount-1) ixcell = xcount - 1;

            iycell = int ((system->y[ind] + system->pbc[1]*0.5)/cellsize);
            if (iycell < 0) iycell = 0;
            else if (iycell > ycount-1) iycell = ycount - 1;

            izcell = int ((system->z[ind] + system->pbc[2]*0.5)/cellsize);
            if (izcell < 0) izcell = 0;
            else if (izcell > zcount-1) izcell = zcount - 1;

            linkedlist[ind] = head[ixcell][iycell][izcell];
            head[ixcell][iycell][izcell] = ind;
        }
    }
    return head;
}

void ANALYSIS_FUNCTIONALIZE::compute_void() {
    sel1->anglezs.clear();
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> disp(3,0.0);
    float dist2;
    int izbin;
    float Avogadro = 6.02e23;
    float Atocm = 1.0e-8;
    float dist_crit2 = this->dist_crit*this->dist_crit;
    float kBT = 1.38e-23*300.0;
    float etoCoul = 1.6e-19;
    float permittivity = 8.854187e-12;
    float dist;
    float dr = this->dr;
    int nbins = int(this->dist_crit/dr);
    this->iframe += 1;

    //cout <<"dist_crit: " << this->dist_crit<< "  dr: " << dr << " rdf bin: " << nbins << endl;
    //cout << "rdf nbins: " << nbins << endl;

    vector<float> rdf(nbins,0.0);

    if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    if (sel2->NATOM == 0) error1.error_exit("ERROR: sel2 doesn't contain any atoms!");
    //cout << "sel1->NATOM: " << sel1->NATOM << endl; //for debug purpose
    //cout << "sel2->NATOM: " << sel2->NATOM << endl; //for debug purpose
//    cout << "M_PI" << M_PI << endl;// for debug purpose

    xcount = int(system->pbc[0]/cellsize);
    ycount = int(system->pbc[1]/cellsize);
    zcount = int(system->pbc[2]/cellsize);

    linkedlist.resize(system->NATOM);

    vector<vector<vector<int>>> head1 = head_cell(sel1->segments_ind);
    vector<vector<vector<int>>> head2 = head_cell(sel2->segments_ind);

    int xcount = head1.size();
    int ycount = head1[0].size();
    int zcount = head1[0][0].size();

    
 
    for (auto &segment:sel1->segments_ind) {
	for (int ind : segment) {
	    r[0] = system->x[ind];
	    r[1] = system->y[ind];
	    r[2] = system->z[ind];

    //cout << "sel1->atom_index: " << system->atom_index[ind] << endl; //for debug purpose
    //cout << "sel1->atomtype: " << system->atomtype[ind] << endl; //for debug purpose

            for (auto &segment:sel2->segments_ind) {
	        for (int ind1 : segment) {
	            r1[0] = system->x[ind1];
	            r1[1] = system->y[ind1];
	            r1[2] = system->z[ind1];
                    
	            disp = getDistPoints(r, r1);
                    dist2 = disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
                    if (dist2 < dist_crit2) {
                        dist = sqrt(dist2);
                        int ibin = int(dist/dr);
                        this->rdf_count[ibin] += 1.0;
                    }
 
                }
            }
	}
    }


    for (int ibin = 0; ibin <= nbins; ibin++) {
        float drp = dr * ibin;
        float drn = dr * (ibin+1);
        rdf[ibin] = this->rdf_count[ibin]/(this->iframe * 4.0/3.0 * PI * (drn*drn*drn - drp*drp*drp));
 //       cout <<"rdf: " << rdf[ibin] << endl;
    }

    for (int ibin = 0; ibin <= nbins; ibin++) {
        float nideal = rdf[nbins-1];
        rdf[ibin] /= nideal; 
    }
/*
*/


}


ANALYSIS_FUNCTIONALIZE::~ANALYSIS_FUNCTIONALIZE()
{
    system = NULL;
    sel1 = NULL;
    sel2 = NULL;
    this->rdf_count.clear();
}

