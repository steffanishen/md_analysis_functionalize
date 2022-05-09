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
#include <string.h>
#include <algorithm>

#include <math.h>
#include "analysis_patch_no_order.hpp"


#define PI 3.14159265

using namespace std;

ANALYSIS_PATCH_NO_ORDER::ANALYSIS_PATCH_NO_ORDER(PSF *system, GROUP *sel1, GROUP *sel2, vector<GROUP*> sels, int vector1d, int vector2d, int voidf, string input_cluster_name, string output_cluster_name, string filename, float dist_crit, float dr)
{
    this->system = system;
    this->sels = sels;

//    this->sel1 = sel1;
//    this->sel2 = sel2;

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
    this->system->crosslinking_flag.resize(system->NATOM,0);
    this->system->functionalizing_flag.resize(system->NATOM,0);
    //char fileSpec[filename.length()+1];
    //snprintf(fileSpec, sizeof(fileSpec),"%s",filename.c_str());
    this->file_temp = new ofstream(filename);
   // fstream *input_cluster = new fstream(input_cluster_name);

    input_cluster = new ifstream(input_cluster_name);
/*
    try
   {
      input_cluster.open(input_cluster_name.c_str());
      if (input_cluster.fail()) throw input_cluster_name;  // the exception being checked
      /*
      while (input_cluster.good())  // check next character
      {
         cout << descrip << ' ' << price << endl;
         inFile >> descrip >> price;
      }
      inFile.close();
      return 0;
   }
   catch(string e)
   {
      cout << e << " was not successfully opened.\n Please check that the file currently exists." << endl;
      exit(1);
   }
*/

/*
    input_cluster->exceptions(std::ifstream::failbit);
    try
    {
        input_cluster->open(input_cluster_name,ios::in);
    }
    catch(ifstream::failure e)
    {
        cerr << "Exception opening/reading file '" << input_cluster_name << "' : " << std::endl;
        cerr << "Please check the path of the file and if it exists." << endl;
    }

*/
    
    output_cluster = new ofstream(output_cluster_name);
   // clusters.resize(system->segments.size(), vector<CLUSTER*>(system->segments[0].size(), new CLUSTER(-1)));
    residue_cluster_ind.resize(system->segments.size(), vector<int>(system->segments[0].size() ));

}

void ANALYSIS_PATCH_NO_ORDER::init() {
}

vector<vector<vector<int>>> ANALYSIS_PATCH_NO_ORDER::head_cell(vector<vector<int>> segments_ind) {
    vector<vector<vector<int>>> head;
    head.resize(xcount,vector<vector<int>>(ycount,vector<int>(zcount,-1)));//,vector<int>(zcount)));

    int ixcell,iycell,izcell;
     for (auto &segment:segments_ind) {
	    for (int ind : segment) {
            ixcell = int ((system->x[ind] + system->pbc[0]*0.5)/cellsize);
            if (ixcell < 0) ixcell = 0;
            else if (ixcell > xcount-1) ixcell = xcount - 1;

            iycell = int ((system->y[ind] + system->pbc[2]*0.5)/cellsize);
            if (iycell < 0) iycell = 0;
            else if (iycell > ycount-1) iycell = ycount - 1;

            izcell = int ((system->z[ind] + system->pbc[5]*0.5)/cellsize);
            if (izcell < 0) izcell = 0;
            else if (izcell > zcount-1) izcell = zcount - 1;

            linkedlist[ind] = head[ixcell][iycell][izcell];
            head[ixcell][iycell][izcell] = ind;
        }
    }
    return head;
}

int ANALYSIS_PATCH_NO_ORDER::neighbor_cell_ind(int i, int i_incr, int n) {
    int i_adj;
    i_adj = i + i_incr;
    if (i_adj < 0) i_adj += n;
    if (i_adj >= n) i_adj -= n;
    return i_adj;
}


string ANALYSIS_PATCH_NO_ORDER::patchtype(string name1, string resname1, string name2, string resname2) {
    string patchtype;
    if (name1 == "C2" && resname1 == "STYR" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PST1";
    else if (name1 == "C2" && resname1 == "STYR" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PST2";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PST3";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PST4";
    else if (name1 == "C2" && resname1 == "DVB" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PST2";
    else if (name1 == "C2" && resname1 == "STYR" && name2 == "C1" && resname2 == "DVB" ) patchtype = "PST2";
    else if (name1 == "C1" && resname1 == "DVB" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PST3";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C2" && resname2 == "DVB" ) patchtype = "PST3";
    else if (name1 == "C3" && resname1 == "DVB" && name2 == "C2" && resname2 == "STYR" ) patchtype = "PSTVB7";
    else if (name1 == "C2" && resname1 == "STYR" && name2 == "C3" && resname2 == "DVB" ) patchtype = "PSTVB8";
    else if (name1 == "C4" && resname1 == "DVB" && name2 == "C1" && resname2 == "STYR" ) patchtype = "PSTVB3";
    else if (name1 == "C1" && resname1 == "STYR" && name2 == "C4" && resname2 == "DVB" ) patchtype = "PSTVB4";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO31";
    else if (name1 == "CE1"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO11";
    else if (name1 == "CD1"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO21";
    else if (name1 == "CD2"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO41";
    else if (name1 == "CE2"                      && name2 == "C1" && resname2 == "NC4" ) patchtype = "ORTHO51";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO32";
    else if (name1 == "CE1"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO12";
    else if (name1 == "CD1"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO22";
    else if (name1 == "CD2"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO42";
    else if (name1 == "CE2"                      && name2 == "C2" && resname2 == "NC4" ) patchtype = "ORTHO52";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO33";
    else if (name1 == "CE1"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO13";
    else if (name1 == "CD1"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO23";
    else if (name1 == "CD2"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO43";
    else if (name1 == "CE2"                      && name2 == "C3" && resname2 == "NC4" ) patchtype = "ORTHO53";
    else if (name1 == "CG" && resname1 == "STYR" && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO34";
    else if (name1 == "CE1"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO14";
    else if (name1 == "CD1"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO24";
    else if (name1 == "CD2"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO44";
    else if (name1 == "CE2"                      && name2 == "C4" && resname2 == "NC4" ) patchtype = "ORTHO54";
    else error1.error_exit("Cannot find the patch type!!"); 
    return patchtype;
}

void ANALYSIS_PATCH_NO_ORDER::flagallinres(int segid, int resid) {
    segid -= 1;
    resid -= 1;
    int natoms = system->segments[segid][resid].size();
    for (int i = 0; i < natoms; i++) {
        int flag_index = system->segments[segid][resid][i];
        system->crosslinking_flag[flag_index] = 1;
    }
}

void ANALYSIS_PATCH_NO_ORDER::flagallinresifcrosslinked(int segid, int resid) {
    segid -= 1;
    resid -= 1;
    int natoms = system->segments[segid][resid].size();
    int flag_temp = 0;
    for (int i = 0; i < natoms; i++) {
        int flag_index = system->segments[segid][resid][i];
        if (system->atomname[flag_index] == "C1" && system->charge[flag_index] > -0.15) flag_temp = 1;
        else if (system->atomname[flag_index] == "C2" && system->charge[flag_index] > -0.15) flag_temp = 1;
        else if (system->atomname[flag_index] == "C3" && system->charge[flag_index] > -0.15) flag_temp = 1;
        else if (system->atomname[flag_index] == "C4" && system->charge[flag_index] > -0.15) flag_temp = 1;
    }

    if (flag_temp == 1) {
        for (int i = 0; i < natoms; i++) {
            int flag_index = system->segments[segid][resid][i];
            system->crosslinking_flag[flag_index] = 1;
        }
    }
}

void ANALYSIS_PATCH_NO_ORDER::flagfunctionalization(int segid, int resid) {
    segid -= 1;
    resid -= 1;
    int natoms = system->segments[segid][resid].size();
    int flag_temp = 0;
    for (int i = 0; i < natoms; i++) {
        int flag_index = system->segments[segid][resid][i];
        if ((system->atomname[flag_index] == "CD1" || system->atomname[flag_index] == "CD2" || system->atomname[flag_index] == "CE1" || system->atomname[flag_index] == "CE2") && system->charge[flag_index] > -0.05) flag_temp = 1;
        else if (system->atomname[flag_index] == "CG" && system->resname[flag_index] == "STYR" && system->charge[flag_index] > -0.05 ) flag_temp = 1;
    }

    if (flag_temp == 1) {
        for (int i = 0; i < natoms; i++) {
            int flag_index = system->segments[segid][resid][i];
            system->functionalizing_flag[flag_index] = 1;
        }
    }
}

bool ANALYSIS_PATCH_NO_ORDER::is_empty(std::ifstream *pFile)
{
    return pFile->peek() == std::ifstream::traits_type::eof();
}


void ANALYSIS_PATCH_NO_ORDER::initialize_clusters() {
    int icluster = 0;
    nresidues = 0;
    for (int i = 0; i < system->segments.size(); i++) {
        for (int i1 = 0; i1 < system->segments[i].size(); i1++) {
            residue_cluster_ind[i][i1] = icluster;
            CLUSTER *cluster_temp = new CLUSTER(icluster);
            cluster_temp->residue_members.push_back({i,i1});
            clusters.push_back(cluster_temp);
            icluster++;
            nresidues++;
        }
    }
    nclusters = icluster;

}

void ANALYSIS_PATCH_NO_ORDER::organize_clusters(int cluster_id, int cluster_main_id) {
    if (clusters[cluster_id]->kid_clusters.size() > 0) {
        clusters[cluster_main_id]->residue_members.insert(clusters[cluster_main_id]->residue_members.end(),clusters[cluster_id]->residue_members.begin(),clusters[cluster_id]->residue_members.end());
        for (int i = 0; i < clusters[cluster_id]->kid_clusters.size(); i++) {
            organize_clusters(clusters[cluster_id]->kid_clusters[i]->cluster_ind, cluster_main_id);
        }
    }
}

void ANALYSIS_PATCH_NO_ORDER::reduce_clusters() {
    int icluster = 0;
    int segid_temp;
    int resid_temp;
    for (int i = 0; i < clusters.size(); i++) {
        if (clusters[i]->parent_cluster_ind == -1) {
            clusters[i]->cluster_ind = icluster;
            for (int i1 = 0; i1 < clusters[i]->residue_members.size(); i1++) {
                segid_temp = clusters[i]->residue_members[i1][0];
                resid_temp = clusters[i]->residue_members[i1][1];
                residue_cluster_ind[segid_temp][resid_temp] = icluster;
            }
            clusters[i]->kid_clusters.clear();
            icluster++;
        }
    }

    auto it = clusters.begin();
    while (it != clusters.end()) {
        if ( (*it)->parent_cluster_ind != -1) {
//            cout << "parent_cluster_ind: " << (*it)->parent_cluster_ind << endl;
            it = clusters.erase(it);
        } else {
            ++it;
        }
    }
    nclusters = clusters.size();

}



void ANALYSIS_PATCH_NO_ORDER::find_initial_clusters() {
    for (int i = 0; i < system->NBONDS; i++) {
        int iatom = system->ibond[i][0];
            int jatom = system->ibond[i][1];
            int segid_temp0 = system->segid[iatom];
            int resid_temp0 = system->resid[iatom];
            int segid_temp1 = system->segid[jatom];
            int resid_temp1 = system->resid[jatom];
            int cluster1;
            int cluster2;
            if (segid_temp0 != segid_temp1 || resid_temp0 != resid_temp1) {
                cluster1 = residue_cluster_ind[segid_temp0][resid_temp0];
                cluster2 = residue_cluster_ind[segid_temp1][resid_temp1]; 
                if (cluster1 < cluster2) {
                    clusters[cluster2]->parent_cluster_ind = cluster1;
                    clusters[cluster1]->kid_clusters.push_back(clusters[cluster2]);
                } else {
                    clusters[cluster1]->parent_cluster_ind = cluster2;
                    clusters[cluster2]->kid_clusters.push_back(clusters[cluster1]);
                }
            }
    }

    for (int i = 0; i < clusters.size(); i++) {
        int cluster_ind_main;
        if (clusters[i]->parent_cluster_ind == -1) {
            cluster_ind_main = clusters[i]->cluster_ind;
            organize_clusters(i,cluster_ind_main);
        }
    }
    reduce_clusters();
}



void ANALYSIS_PATCH_NO_ORDER::compute_void() {
    vector<float> r(3,0.0);
    vector<float> r1(3,0.0);
    vector<float> disp(3,0.0);
    vector<float> distances;
    vector<vector<int>> candidate_pairs;
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
    vector<vector<vector<vector<int>>>> heads;

    //cout <<"dist_crit: " << this->dist_crit<< "  dr: " << dr << " rdf bin: " << nbins << endl;
    //cout << "rdf nbins: " << nbins << endl;

    vector<float> rdf(nbins,0.0);
    int nsels = sels.size();
    if ( (nsels % 2) != 0) error1.error_exit("ERROR: Odd number of groups. Please select pairs of groups for patches!!");
    int npairs = nsels/2;

    for (int i = 0; i < nsels; i++) {
        if (sels[i]->NATOM == 0) error1.error_exit("ERROR: sels doesn't contain any atoms!");
        //cout << "segments_ind size: " << sels[i]->segments_ind.size()<< endl;
    }
    
    //if (sel1->NATOM == 0) error1.error_exit("ERROR: sel1 doesn't contain any atoms!");
    //if (sel2->NATOM == 0) error1.error_exit("ERROR: sel2 doesn't contain any atoms!");



    if (system->pbc[0] < 0.01 || system->pbc[2] < 0.01 || system->pbc[5] < 0.01 ) error1.error_exit("ERROR: Box size not specified!");
    //cout << "sel1->NATOM: " << sel1->NATOM << endl; //for debug purpose
    //cout << "sel2->NATOM: " << sel2->NATOM << endl; //for debug purpose
//    cout << "M_PI" << M_PI << endl;// for debug purpose

    xcount = int(system->pbc[0]/cellsize);
    ycount = int(system->pbc[2]/cellsize);
    zcount = int(system->pbc[5]/cellsize);

    linkedlist.resize(system->NATOM);

//    vector<vector<vector<int>>> head1 = head_cell(sel1->segments_ind);
//    vector<vector<vector<int>>> head2 = head_cell(sel2->segments_ind);

    for (int i = 0; i < nsels; i++) {
        vector<vector<vector<int>>> headtemp = head_cell(sels[i]->segments_ind);
        heads.push_back(headtemp);
    }

    int xcount = heads[0].size();
    int ycount = heads[0][0].size();
    int zcount = heads[0][0][0].size();

// Consider if residues are NC4 functionalized or not
    for (int i = 0; i < system->segments.size(); i++) {
        for (int i1 = 0; i1 < system->segments[i].size(); i1++) {
            int ind1 = system->segments[i][i1][0];
            int segid_temp0 = system->segid[ind1];
            int resid_temp0 = system->resid[ind1];
            if (system->resname[ind1] == "NC4") flagallinresifcrosslinked(segid_temp0, resid_temp0);
            if (system->resname[ind1] == "STYR" || system->resname[ind1] == "DVB"  ) flagfunctionalization(segid_temp0,resid_temp0);
        }
    }

// Initialize the residue clusters
    initialize_clusters();

// Merge the initial clusters
    find_initial_clusters(); 

    cout << "nresidues: " << nresidues << "; nclusters: " << nclusters << endl;

// Find the distance between possible croslinking pairs
    for (int iselpair = 0; iselpair < nsels/2; iselpair++) {
        int isel1 = iselpair*2;
        int isel2 = isel1 + 1;
        for (int i = 0; i < xcount; i++) {
            for (int j = 0; j < ycount; j++) {
                for (int k = 0; k < zcount; k++) {
                    int ind1 = heads[isel1][i][j][k];
                    while (1) {
                        if (ind1 < 0) break;
                        r[0] = system->x[ind1];
	                    r[1] = system->y[ind1];
	                    r[2] = system->z[ind1];
                        for (int iprime = -1; iprime <=1; iprime++) {
                            int i2 = neighbor_cell_ind(i,iprime,xcount);
                            for (int jprime = -1; jprime <=1; jprime++) {
                                int j2 = neighbor_cell_ind(j,jprime,ycount);
                                for (int kprime = -1; kprime <=1; kprime++) {
                                    int k2 = neighbor_cell_ind(k,kprime,zcount);
                                    int ind2 = heads[isel2][i2][j2][k2];
                                    while (1) {
                                        if (ind2 < 0) break;
                                        if (!(system->segid[ind1] == system->segid[ind2] && system->resid[ind1] == system->resid[ind2])) {

                                            r1[0] = system->x[ind2];
    	                                    r1[1] = system->y[ind2];
    	                                    r1[2] = system->z[ind2];

                                            disp = getDistPoints(r, r1);
                                            dist2 = disp[0]*disp[0] + disp[1]*disp[1] + disp[2]*disp[2];
                                            dist = sqrt(dist2);
                                    
                                            distances.push_back(dist);
                                            candidate_pairs.push_back({ind1,ind2});

                                        }
                                        ind2 = linkedlist[ind2];
                                    }
                                }
                            }
                        }
                        ind1 = linkedlist[ind1];
                    }
                }
            }
        }
    }

    //vector<float> vector_test = distances;

    vector<int> a_ind = heapSort(distances,distances.size());




    for (int i = 0; i < a_ind.size(); i++ ) {
        int i_sorted = a_ind[i];
        float local_dist = distances[i];
        int ind1 = candidate_pairs[i_sorted][0];
        int ind2 = candidate_pairs[i_sorted][1];
        int segid1 = system->segid[ind1];        
        int segid2 = system->segid[ind2]; 
        int resid1 = system->resid[ind1];       
        int resid2 = system->resid[ind2];       

        int segid_temp;
        int resid_temp;
        
        int functionalizing = 0;
        int functionalized = 0;

        if (local_dist < dist_crit && system->crosslinking_flag[ind1] == 0 && system->crosslinking_flag[ind2] == 0) {
        /*    if (system->resname[ind1] == "NC4") {
                functionalizing = 1;
                segid_temp = system->segid[ind1];
                resid_temp = system->resid[ind1];
                flagallinresifcrosslinked(segid_temp, resid_temp);
            }
            if (system->resname[ind2] == "NC4") {
                functionalizing = 1;
                segid_temp = system->segid[ind2];
                resid_temp = system->resid[ind2];
                flagallinresifcrosslinked(segid_temp, resid_temp);
            }
        */

            if (system->resname[ind1] == "NC4" || system->resname[ind2] == "NC4") {
                functionalizing = 1;
                if (system->functionalizing_flag[ind1] == 1 || system->functionalizing_flag[ind2] == 1) functionalized = 1;
            }

            //if (system->crosslinking_flag[ind1] == 0 && system->crosslinking_flag[ind2] == 0) {
                if (! (functionalizing == 1 && functionalized == 1) ) {
                    string patchtype = this->patchtype(system->atomname[ind1],system->resname[ind1],system->atomname[ind2],system->resname[ind2]);
                    *this->file_temp << "patch " << patchtype << " " << segid1 << ":" << resid1 << " " << segid2 << ":" << resid2 << endl;
                    system->crosslinking_flag[ind1] = 1;
                    system->crosslinking_flag[ind2] = 1;
                }
            //}

            if (system->resname[ind1] == "NC4") {
                segid_temp = system->segid[ind1] ;
                resid_temp = system->resid[ind1] ;
                flagallinres(segid_temp, resid_temp);
            }
            if (system->resname[ind2] == "NC4") {
                segid_temp = system->segid[ind2] ;
                resid_temp = system->resid[ind2] ;
                flagallinres(segid_temp, resid_temp);
            }
        }
    }

    cout << "N_pairs: " << a_ind.size() << endl; 



}


ANALYSIS_PATCH_NO_ORDER::~ANALYSIS_PATCH_NO_ORDER()
{
    system = NULL;
    //sel1 = NULL;
    //sel2 = NULL;
    this->rdf_count.clear();
//    fclose(this->outfile_box);
    this->file_temp->close();
}

