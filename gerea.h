/***************************************
*
* GEREA
* Gene Expression Reulator Enrhchment Analysis
* This program was designed to
* search the enriched gene expresion
* regulators
* Mod Data: 10/20/2018
* Author: Tinghua Huang
*
***************************************/

#ifndef CIRCLE_H
#define CIRCLE_H

using namespace std;

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <set>
#include <vector>
#include "math.h"
#include <stdlib.h>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
//#include <gmp.h>
//#include "crossprob.h"

//classes
class gene {
public:
    //gene(string id);
    string ref_id;
    string name;
    string desc;
    string cls;
    string expresion; //up, down, nodiff, cls
    vector<double> expressions;
    double fc;
    double qvalue;
    vector<double> fcs;
    vector<double> qvalues;
};

//class transf;

class target {
public:
    target(string id);
    //virtual ~target();
    string ref_id;
    string name;
    string desc;
    bool present;
    gene *gene_pt;
    //transf *parent;
    string rel_type;
    string dri_type;
    //string output;  //diff, up, down, stable
    //string glass; //class, not_assigned
    double fc;
    double qvalue;
};

//class regulon;

/*
class target {
public:
    rtarget(string id);
    //virtual ~target();
    string ref_id;
    string name;
    string desc;
    bool present;
    gene *gene_pt;
    //regulon *parent;
    string rel_type;
    string dri_type;
    //string output;  //diff, up, down, stable
    //string glass; //class, not_assigned
    double fc;
    double qvalue;
};
*/

class bionstat {
public:
    //bionstat();
    //map<string, long> target_n;
    int type;
    int an, bn, cn, dn, en, fn, gn, hn, in, jn, kn, ln, mn, qn, rn;
    double p0;
    double p1;
    double p2;
    vector<double> p1s;
    vector<double> p2s;
    double pvalue;
    double fdr_BH;
    int get_tar_n();
    bool p_nequal0();
    double fdr_BH_1;
    double fdr_BH_2;
    double fdr_BH_0;
    double pup;
    double pdown;
    void cal_pup_pdown();
};

class transf {
public:
    //transf(string id);
    //virtual ~transf();
    string ref_id;
    string name;
    string desc;
    string rel_type;
    vector<target> targets;
    //bionstat stat;
    void get_tids(set<string> &tids);
    void get_ptids(set<string> &tids);
    void find_targets(string target_id, vector<target*> &tars_pt);
    void encode_network1(vector<string> &ecd_strs);
    void encode_network2(vector<string> &ecd_strs);
    void encode_network3(vector<string> &ecd_strs);
    void encode_network4(vector<string> &ecd_strs);
    void encode_network5(vector<string> &ecd_strs, double q_threshold);
    void cal_dtype(string p_rtype);
    void cal_dtype2(string p_rtype);
    void cal_upregulation(vector<double> &qbackground, vector<double> &qdiff);
    void cal_downregulation(vector<double> &qbackground, vector<double> &qdiff);
    void cal_upregulation(vector<double> &qtarget);
    void cal_upregulation(vector<double> &qtarget, unsigned long perm_i);
    void cal_downregulation(vector<double> &qtarget);
    void cal_downregulation(vector<double> &qtarget, unsigned long perm_i);
    //void cal_abde(int &an, int &bn, int &dn, int &en, double qvalue_threshold);
    void cal_ac1(int &an, int &cn);
    void cal_ace(int &an, int &cn, int &en);
    void cal_beh(int &bn, int &en, int &hn);
    void cal_adg(int &an, int &dn, int &gn);
    void cal_be(int &bn, int &en);
    void cal_ad(int &an, int &dn);
    //
    void cal_adg2(int &an, int &dn, int &gn, int &jn, int &mn);
    void cal_beh2(int &bn, int &en, int &hn, int &kn, int &qn);
};

class regulon {
public:
    //regulon(string id);
    //virtual ~regulon();
    string ref_id;
    string name;
    string desc;
    vector<target> targets;
    vector<gene*> _genes;
    map<string, transf> transfs;
    bionstat stat;
    void loading_data(map<string, gene> &genes);
    void find_targets(string target_id, vector<target*> &tars_pt);
    void get_tids(set<string> &tids);
    void get_ptids(set<string> &tids);
    void encode_network1(vector<string> &ecd_strs);
    void encode_network2(vector<string> &ecd_strs);
    void encode_network3(vector<string> &ecd_strs);
    void encode_network4(vector<string> &ecd_strs);
    void encode_network5(vector<string> &ecd_strs, double q_threshold);
    void cal_dtype();
    void cal_dtype2();
    //void cal_abde(int &an, int &bn, int &dn, int &en, double qvalue_threshold);
    //void cal_cf(int &cn, int &fn, double qvalue_threshold);
    void cal_upregulation(vector<double> &qbackground, vector<double> &qdiff);
    void cal_downregulation(vector<double> &qbackground, vector<double> &qdiff);
    void cal_upregulation(vector<double> &qtarget);
    void cal_upregulation(vector<double> &qtarget, unsigned long perm_i);
    void cal_downregulation(vector<double> &qtarget);
    void cal_downregulation(vector<double> &qtarget, unsigned long perm_i);
    void cal_ac1(int &an, int &cn);
    void cal_bd1(int &bn, int &dn);
    void cal_ace(int &an, int &cn, int &en);
    void cal_bdf(int &bn, int &dn, int &fn);
    void cal_cfi(int &cn, int &fn, int &in);
    void cal_beh(int &bn, int &en, int &hn);
    void cal_adg(int &an, int &dn, int &gn);
    void cal_cf(int &cn, int &fn);
    void cal_be(int &bn, int &en);
    void cal_ad(int &an, int &dn);
    //
    void cal_adg2(int &an, int &dn, int &gn, int &jn, int &mn);
    void cal_beh2(int &bn, int &en, int &hn, int &kn, int &qn);
    void cal_cfi2(int &cn, int &fn, int &in, int &ln, int &rn);
};

class data {
public:
    data();
    string ref_id;
    string name;
    string desc;
    int data_seto;
    int data_type;
    double q_threshold;
    double fc_threshold;
    string class_str;
    map<string, gene> genes;
    void data_config(string &tmp_str);
    void save_data1(string &tmp_str);
    void save_data2(string &tmp_str);
    void save_data_title(string &tmp_str);
    void save_data_class(string &tmp_str);
    void save_data_exp(string &tmp_str);
    void cal_exp1();
    void cal_exp2();
    vector<double> get_qbackground();
    vector<double> get_qbackground(unsigned long perm_i);
    //gene* find_gene(string gene_id);
    //bool exi_gene(string gene_id);
    //void get_did(set<string> &tids);
    //void get_uid(set<string> &tids);
    vector<string> titles;
    vector<vector<int>> rdm_classes;
    vector<int> sample_classes;
    unsigned long sample_n;
    void gen_rdm_classes(unsigned long permutation_n);
    vector<int> class_set;
    map<int, unsigned long> class_set_n;
    void cal_fcs(unsigned long permutation_n);
    void cal_pvalues(unsigned long permutation_n);

};

class network {
public:
    network();
    //virtual ~network();
    int db_type;
    //int data_type;
    map<string, regulon> regulons;
    void db_config(string &tmp_str);
    void save_link1(string &tmp_str);
    void save_link2(string &tmp_str);
    long find_regulon(string regulon_id);
    //void data_config(string &tmp_str);
    void cal_dtype2();
    void cal_dtype();
    //void cal_abcdef(double qvalue_threshold);
    void cal_abcd1();
    void cal_abcdef2();
    void cal_abcdefghi();
    void cal_abcdefghi2();
    void cal_abcdef();
    void BH_correction();
    //void cal_rtype();
    //void run_stat();
    /*void run_stat1();
    void run_stat2();
    void run_stat3();
    void run_stat4();*/
    //void run_stat2();

    //void run_stat4_2();
    void run_stat_KS(int target_n);
    void run_stat_KS2(int target_n, vector<double> &qbackground);
    //void run_stat_KS3(int target_n);
    void BH_correction_2();
    void BH_correction_3();
};

class gerea {
public:
    gerea(string id);
    //virtual ~gerea();
    string session_id;
    string name;
    string desc;
    network network_it;
    data data_it;
    unsigned long target_n;
    unsigned long permutation_n;
    void reading_network(string file_name);
    void reading_data(string file_name);
    void reading_data_exp(string file_name);
    void data_permutation();
    void loading_data2network();
    void run_stat_dt2_dbt2();
    void run_stat_dt2_dbt3();
    void run_stat_KS3();
    void run_stat_KS3(unsigned long perm_i);
    //void run_stat();
    void sortby_fdr(map<string, double> &fdr_map, list<string> &sortby);
    void Print_result(string output_dir);
    void Print_result1(string output_dir);
    void Print_result3(string output_dir);
    void Print_result4(string output_dir);
    void Print_details(string output_dir);
    void Print_details1(string output_dir);
    void Print_details2(string output_dir);
    void Print_details3(string output_dir);
    void Print_details4(string output_dir);

    void Print_result4_2(string output_dir);

    void Print_details5(string output_dir);
    void Print_result5(string output_dir);

    void Print_result6(string output_dir);

};

#endif
