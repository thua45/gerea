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

using namespace std;
#include "gerea.h"

// utility
void trim_string(string &src);
vector<string> str_split(string& src, string delimit);
double BH_Jet(map<string, double> &qmap, string pid);

extern"C" {
	void fexact_(int *nrow, int *ncol, double *table, int *ldtabl, double *expect, double *percnt, double *emin, double *prt, double *pre);
}

vector<string> str_split(string& src, string delimit) {
	string null_subst = "";
	vector<string> strvec;
	if (src.empty() || delimit.empty()) {
		//cout<<"Split: Empty String\n";
		return strvec;
	}
	int deli_len = delimit.length();
	unsigned long index = string::npos, last_search_position = 0;
	while ((index = src.find(delimit, last_search_position)) != string::npos) {
		if (index == last_search_position) {
			// a blank block
			//strvec.push_back(null_subst);
		} else {
			strvec.push_back(src.substr(last_search_position, index - last_search_position));
			last_search_position = index + deli_len;
		}
	}
	string last_one = src.substr(last_search_position);
	if (!last_one.empty()) {
		strvec.push_back(last_one);
	}
	return strvec;
}

void trim_string(string &src) {
	unsigned long pr = string::npos;
	pr = src.find('\r');
	if(pr != string::npos) {
		src.erase(pr, src.size() - pr);
		//src = src.substr(0, pr);
	}

}

target::target(string id) {
	ref_id = id;
	present = false;
}

/*
target::rtarget(string id) {
	ref_id = id;
	present = false;
}
*/

//gene::gene(string id) {
//	ref_id = id;
//}

data::data() {
	//qvalue_threshold = 0.05;
}

gerea::gerea(string id) {
	session_id = id;
}

network::network()
{
	//db_type = 0;
}

int bionstat::get_tar_n() {
	int tn = 0;
	if (type == 1) {
		tn = an + cn;
	}
	else if (type == 2) {
		tn = an + cn + en;
	}
	else if (type == 3) {
		tn = an + dn + bn + en;
	}
	else if (type == 4) {
		tn = an + dn + gn + bn + en + hn;
	}

	return tn;
}

bool bionstat::p_nequal0() {
	bool nequal0 = false;
	if (type == 1) {
		if (p0 > 0 && p1 > 0) {
			nequal0 = true;
		}
	}
	else if (type == 2) {
		if (p0 > 0 && p1 > 0 && p2 > 0) {
			nequal0 = true;
		}
	}
	else if (type == 3) {
		if (p0 > 0 && p1 > 0 && p2 > 0) {
			nequal0 = true;
		}
	}
	else if (type == 4) {
		if (p0 > 0 && p1 > 0 && p2 > 0) {
			nequal0 = true;
		}
	}

	return nequal0;
}

void transf::find_targets(string target_id, vector<target*> &tars_pt) {
    //tars_pt.clear();
    //vector<*target> tars_it;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		/*
		if ((*tar_it).ref_id == target_id && (*tar_it).rel_type != "unknown") {
			tars_pt.push_back(&(*tar_it));
		}
		*/
		if ((*tar_it).ref_id == target_id) {
			tars_pt.push_back(&(*tar_it));
		}
	}
}

void regulon::find_targets(string target_id, vector<target*> &tars_pt) {
    tars_pt.clear();
    //vector<*target> tars_it;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		/*
		if ((*tar_it).ref_id == target_id && (*tar_it).rel_type != "unknown") {
			tars_pt.push_back(&(*tar_it));
		}
		*/
		if ((*tar_it).ref_id == target_id) {
			tars_pt.push_back(&(*tar_it));
		}
	}

	//vector<*target> tmp_tars_pt;
	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
		/*
        if ((*trs_it).second.rel_type == "unknown") {
            continue;
        }
        */
		(*trs_it).second.find_targets(target_id, tars_pt);
		//tars_pt.insert(tmp_tars_pt.begin(), tmp_tars_pt.end());
	}
}

void transf::get_tids(set<string> &tids) {
	tids.clear();
	vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		tids.insert(tar_it->ref_id);
	}
}

void regulon::get_tids(set<string> &tids) {
	tids.clear();
	vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		tids.insert(tar_it->ref_id);
	}

	set<string> tmp_tids;
	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
		(*trs_it).second.get_tids(tmp_tids);
		tids.insert(tmp_tids.begin(), tmp_tids.end());
	}

}

void transf::get_ptids(set<string> &tids) {
	tids.clear();
	vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if (tar_it->present == true) {
			tids.insert(tar_it->ref_id);
		}
	}
}

void regulon::get_ptids(set<string> &tids) {
	tids.clear();
	vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if (tar_it->present == true) {
			tids.insert(tar_it->ref_id);
		}
	}

	set<string> tmp_tids;
	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
		(*trs_it).second.get_ptids(tmp_tids);
		tids.insert(tmp_tids.begin(), tmp_tids.end());
	}
}

void network::db_config(string &tmp_str) {
	vector<string> tmp_vec;
	tmp_vec = str_split(tmp_str, "=");
	if (tmp_vec.size() != 2) {
		cout<<"reading db: invalid line \""<<tmp_str<<"\"\n";
		exit(1);
	}
	if (tmp_vec[0] == "#db_type") {
		db_type = atoi(tmp_vec[1].c_str());
	}

}

void network::save_link1(string &tmp_str) {
	vector<string> tmp_vec;
	tmp_vec = str_split(tmp_str, "\t");
	if (tmp_vec.size() == 0) {
		return;
	}
	if (tmp_vec.size() != 2 && tmp_vec.size() != 3) {
		cout<<"reading db: invalid line \""<<tmp_str<<"\"\n";
		return;
	}
	//long eid;
	if (tmp_vec.size() == 2) {
		//long eid = atoi(tmp_vec[0].c_str());
		string regulon_id = tmp_vec[0];
		string rel_type = "unknown";
		string target_id = tmp_vec[1];

		target new_target(target_id);
		new_target.rel_type = rel_type;
		regulons[regulon_id].targets.push_back(new_target);
	}
	else if (tmp_vec.size() == 3) {
		//cout<<"tmp_vec[0]:"<<tmp_vec[0]<<endl;

		//long eid = atoi(tmp_vec[0].c_str());
		string regulon_id = tmp_vec[0];
		string rel_type1 = "unknown";
		string transf_id = tmp_vec[1];
		string rel_type2 = "unknown";
		string target_id = tmp_vec[2];

		//cout<<"target_id:"<<target_id<<endl;

		target new_target(target_id);
		new_target.rel_type = rel_type2;
        regulons[regulon_id].transfs[transf_id].rel_type = rel_type1;
        regulons[regulon_id].transfs[transf_id].targets.push_back(new_target);
	}
}

void network::save_link2(string &tmp_str) {
	vector<string> tmp_vec;
	tmp_vec = str_split(tmp_str, "\t");
	if (tmp_vec.size() == 0) {
		return;
	}
	if (tmp_vec.size() != 3 && tmp_vec.size() != 5) {
		cout<<"reading db: invalid line \""<<tmp_str<<"\"\n";
		return;
	}
	//long eid;
	if (tmp_vec.size() == 3) {
		//long eid = atoi(tmp_vec[0].c_str());
		string regulon_id = tmp_vec[0];
		string rel_type = tmp_vec[1];
		string target_id = tmp_vec[2];

		target new_target(target_id);
		new_target.rel_type = rel_type;
		regulons[regulon_id].targets.push_back(new_target);
	}
	else if (tmp_vec.size() == 5) {
		//long eid = atoi(tmp_vec[0].c_str());
		string regulon_id = tmp_vec[0];
		string rel_type1 = tmp_vec[1];
		string transf_id = tmp_vec[2];
		string rel_type2 = tmp_vec[3];
		string target_id = tmp_vec[4];

		target new_target(target_id);
		new_target.rel_type = rel_type2;
        regulons[regulon_id].transfs[transf_id].rel_type = rel_type1;
        regulons[regulon_id].transfs[transf_id].targets.push_back(new_target);
	}
}

void gerea::reading_network(string file_name) {
	//string file_name = Data_Cfg.db_dir + Data_Cfg.database + ".Entity";
	//file_name = "network.txt";
	ifstream reader;
	reader.open(file_name.c_str(), ios::in);
	if (reader.fail()) {
		cout<<"cannot open file:\""<<file_name<<"\"\n";
		exit(1);
	}
	string tmp_str;
	vector<string> tmp_vec;

	if (!reader.eof()) {
		getline(reader, tmp_str, '\n');
		trim_string(tmp_str);
		if (tmp_str[0] == '#') {
			//cout<<tmp_str<<endl;

			network_it.db_config(tmp_str);
			//cout<<"db_type:"<<network_it.db_type<<endl;
		}
		else {
			cout<<"unknown database header:\""<<tmp_str<<"\"\n";
			exit(1);
		}
	}

	while (!reader.eof()) {
		getline(reader, tmp_str, '\n');
		trim_string(tmp_str);
		//cout<<tmp_str<<endl;

		if (network_it.db_type == 1) {
			network_it.save_link1(tmp_str);
			//cout<<"saved:"<<tmp_str<<endl;
		}
		else if (network_it.db_type == 2 || network_it.db_type == 3) {
			network_it.save_link2(tmp_str);
		}

	}
	reader.close();
}

void data::cal_exp1() {
	map<string, gene>::iterator gene_it;
	// for each element in the regulons
	for (gene_it = genes.begin(); gene_it != genes.end(); gene_it++) {
		if ((*gene_it).second.cls == class_str) {
			(*gene_it).second.cls = "cls";
		}
		else {
			(*gene_it).second.cls = "---";
		}
	}
}

void data::cal_exp2() {
	map<string, gene>::iterator gene_it;
	// for each element in the regulons
	for (gene_it = genes.begin(); gene_it != genes.end(); gene_it++) {
		if ((*gene_it).second.qvalue <= q_threshold && (*gene_it).second.fc >= fc_threshold) {
			(*gene_it).second.expresion = "up";
		}
		else if ((*gene_it).second.qvalue <= q_threshold && (*gene_it).second.fc <= (-1) * fc_threshold) {
			(*gene_it).second.expresion = "down";
		}
		else {
			(*gene_it).second.expresion = "nodiff";
		}
	}
}

void data::data_config(string &tmp_str) {
	vector<string> tmp_vec;
	tmp_vec = str_split(tmp_str, "=");
	if (tmp_vec.size() != 2) {
		cout<<"reading data: invalid line \""<<tmp_str<<"\"\n";
		exit(1);
	}
	if (tmp_vec[0] == "#data_type") {
		data_type = atoi(tmp_vec[1].c_str());
	}
	if (data_seto == 1 && data_type == 2) {
		cout<<"data_type conflict with -c option!"<<endl;
			exit(1);
	}
	else if (data_seto == 2 && data_type == 1) {
		cout<<"data_type conflict with -q and -f option!"<<endl;
			exit(1);
	}
}

void data::save_data1(string &tmp_str) {
	vector<string> tmp_vec;
	tmp_vec = str_split(tmp_str, "\t");
	if (tmp_vec.size() == 0) {
		return;
	}
	if (tmp_vec.size() != 2) {
		cout<<"invalid data line: \""<<tmp_str<<"\"\n";
		exit(1);
	}
	string gene_id = tmp_vec[0];
	string cls = tmp_vec[1];

	gene new_gene;
	new_gene.cls = cls;
	new_gene.ref_id = gene_id;
	genes[gene_id] = new_gene;
}

void data::save_data2(string &tmp_str) {
	vector<string> tmp_vec;
	tmp_vec = str_split(tmp_str, "\t");
	if (tmp_vec.size() == 0) {
		return;
	}
	if (tmp_vec.size() != 3) {
		cout<<"invalid data line: \""<<tmp_str<<"\"\n";
		exit(1);
	}
	string gene_id = tmp_vec[0];
	float fc = atof(tmp_vec[1].c_str());
	float qvalue = atof(tmp_vec[2].c_str());

	gene new_gene;
	new_gene.fc = fc;
	new_gene.qvalue = qvalue;
	new_gene.ref_id = gene_id;
	genes[gene_id] = new_gene;
}

void gerea::reading_data(string file_name) {
	//string file_name = Data_Cfg.work_dir + Data_Cfg.session_in + ".drp";
	//file_name = "data.txt";
	ifstream reader;
	reader.open(file_name.c_str(), ios::in);
	if (reader.fail()) {
		cout<<"cannot open file:\""<<file_name<<"\"\n";
		exit(1);
	}
	string tmp_str;
	vector<string> tmp_vec;

	if (!reader.eof()) {
		getline(reader, tmp_str, '\n');
		trim_string(tmp_str);
		if (tmp_str[0] == '#') {
			data_it.data_config(tmp_str);
		}
		else {
			cout<<"unknown data header:\""<<tmp_str<<"\"\n";
			exit(1);
		}
	}

	while (!reader.eof()) {
		getline(reader, tmp_str, '\n');
		trim_string(tmp_str);
		if (data_it.data_type == 1) {
			data_it.save_data1(tmp_str);
		}
		else if (data_it.data_type == 2) {
			data_it.save_data2(tmp_str);
		}

	}

	if (data_it.data_type == 1) {
		//cout<<"cal_exp"<<endl;

		data_it.cal_exp1();
	}
	else if (data_it.data_type == 2) {
		//cout<<"cal_exp"<<endl;

		data_it.cal_exp2();
	}
}

void regulon::loading_data(map<string, gene> &genes) {
    map<string, gene>::iterator gene_it;
    vector<target*> tars_pt;
    vector<target*>::iterator tars_pt_it;
	// for each element in the regulons
	for (gene_it = genes.begin(); gene_it != genes.end(); gene_it++) {
        tars_pt.clear();
	    find_targets((*gene_it).second.ref_id, tars_pt);
	    if (tars_pt.size() == 0) {
            _genes.push_back(&(*gene_it).second);
            continue;
	    }
	    else {
	    	//cout<<"found_tars:"<<tars_pt.size()<<endl;

            for (tars_pt_it = tars_pt.begin(); tars_pt_it != tars_pt.end(); tars_pt_it++) {
                (**tars_pt_it).present = true;
                (**tars_pt_it).gene_pt = &(*gene_it).second;
            }
	    }
	}
}

void gerea::loading_data2network() {
    map<string, regulon>::iterator reg_it;
	//map<string, transf>::iterator trs_it;
	// for each element in the regulons
	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		(*reg_it).second.loading_data(data_it.genes);
	}
}

double BH_Jet(map<string, double> &qmap, string pid) {
	double pvalue;
	double qBH;
	map<string, double>::iterator qrt;
	pvalue = qmap[pid];
	long qi = 1;
	for (qrt=qmap.begin(); qrt!=qmap.end(); qrt++) {
		if ((*qrt).second<pvalue) {
			qi++;
		}
	}
	qBH = pvalue*double(qmap.size())/double(qi);
	if (qBH>1) {
		qBH = 1;
	}
	else if (qBH<0) {
		qBH = 0;
	}
	//cout<<qBH<<endl;
	return qBH;
}

void network::BH_correction_2() {
	map<string, double> qdiff_map;
	map<string, regulon>::iterator reg_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		if ((*reg_it).second.stat.get_tar_n() == 0) {
			continue;
		}
		string reg_id = (*reg_it).first;
		//double pvalue = (*reg_it).second.stat.pvalue;
		qdiff_map[(*reg_it).first] = (*reg_it).second.stat.p0;
	}

	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		(*reg_it).second.stat.fdr_BH_0 = BH_Jet(qdiff_map, (*reg_it).first);
	}
	//p1
	qdiff_map.clear();
    for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		if ((*reg_it).second.stat.get_tar_n() == 0) {
			continue;
		}
		string reg_id = (*reg_it).first;
		//double pvalue = (*reg_it).second.stat.pvalue;
		qdiff_map[(*reg_it).first] = (*reg_it).second.stat.p1;
	}

	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		(*reg_it).second.stat.fdr_BH_1 = BH_Jet(qdiff_map, (*reg_it).first);
	}
	//p2
	qdiff_map.clear();
    for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		if ((*reg_it).second.stat.get_tar_n() == 0) {
			continue;
		}
		string reg_id = (*reg_it).first;
		//double pvalue = (*reg_it).second.stat.pvalue;
		qdiff_map[(*reg_it).first] = (*reg_it).second.stat.p2;
	}

	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		(*reg_it).second.stat.fdr_BH_2 = BH_Jet(qdiff_map, (*reg_it).first);
	}
}

void network::BH_correction() {
	map<string, double> qdiff_map;
	map<string, regulon>::iterator reg_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		if ((*reg_it).second.stat.get_tar_n() == 0) {
			continue;
		}
		string reg_id = (*reg_it).first;
		//double pvalue = (*reg_it).second.stat.pvalue;
		qdiff_map[(*reg_it).first] = (*reg_it).second.stat.pvalue;
	}

	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		(*reg_it).second.stat.fdr_BH = BH_Jet(qdiff_map, (*reg_it).first);
	}
}

void transf::cal_ac1(int &an, int &cn) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).cls == "cls") {
            an++;
		}
		else if ((*(*tar_it).gene_pt).cls == "---") {
			cn++;
		}
		else {
			cout<<"unknown class type!"<<endl;
			exit(1);
		}
	}
}

void regulon::cal_ac1(int &an, int &cn) {
	//begin
	an = 0; cn = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).cls == "cls") {
            an++;
		}
		else if ((*(*tar_it).gene_pt).cls == "---") {
			cn++;
		}
		else {
			cout<<"unknown class type!"<<endl;
			exit(1);
		}
	}

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_ac1(an, cn);
	}
}

void regulon::cal_bd1(int &bn, int &dn) {
    //map<string, gene>::iterator gene_it;
    bn = 0; dn = 0;
    vector<gene*>::iterator genes_pt_it;
    for (genes_pt_it = _genes.begin(); genes_pt_it != _genes.end(); genes_pt_it++) {
        if ((**genes_pt_it).cls == "cls") {
            bn++;
        }
        else if ((**genes_pt_it).cls == "---") {
            dn++;
        }
        else {
			cout<<"unknown class type!"<<endl;
			exit(1);
		}
    }
}

void network::cal_abcd1() {
	map<string, regulon>::iterator reg_it;
	bionstat *stat_pt;
	//map<string, transf>::iterator trs_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		stat_pt = &(*reg_it).second.stat;
		(*stat_pt).type = 1;
		(*reg_it).second.cal_ac1((*stat_pt).an, (*stat_pt).cn);
		(*reg_it).second.cal_bd1((*stat_pt).bn, (*stat_pt).dn);
	}
}

void transf::cal_ace(int &an, int &cn, int &en) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up") {
            an++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down") {
			cn++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff") {
			en++;
		}
	}
}

void regulon::cal_ace(int &an, int &cn, int &en) {
	//begin
	an = 0; cn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up") {
            an++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down") {
			cn++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff") {
			en++;
		}
	}

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_ace(an, cn, en);
	}
}

void regulon::cal_bdf(int &bn, int &dn, int &fn) {
    //map<string, gene>::iterator gene_it;
    bn = 0; dn = 0; fn = 0;
    vector<gene*>::iterator genes_pt_it;
    for (genes_pt_it = _genes.begin(); genes_pt_it != _genes.end(); genes_pt_it++) {
        if ((**genes_pt_it).expresion == "up") {
            bn++;
        }
        else if ((**genes_pt_it).expresion == "down") {
            dn++;
        }
        else if ((**genes_pt_it).expresion == "nodiff") {
            fn++;
        }
    }
}

void network::cal_abcdef2() {
	map<string, regulon>::iterator reg_it;
	bionstat *stat_pt;
	//map<string, transf>::iterator trs_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		stat_pt = &(*reg_it).second.stat;
		(*stat_pt).type = 2;
		(*reg_it).second.cal_ace((*stat_pt).an, (*stat_pt).cn, (*stat_pt).en);
		(*reg_it).second.cal_bdf((*stat_pt).bn, (*stat_pt).dn, (*stat_pt).fn);
	}
}

void network::run_stat1() {
	map<string, regulon>::iterator reg_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		int an, bn, cn, dn, en, fn;
		an = (*reg_it).second.stat.an;
		bn = (*reg_it).second.stat.bn;
		cn = (*reg_it).second.stat.cn;
		dn = (*reg_it).second.stat.dn;

		if (an + bn + cn + dn == 0) {
			(*reg_it).second.stat.p0 = 0;
		}
		else {
			(*reg_it).second.stat.p0 = double(an + cn) / double(an + bn + cn + dn);
		}

		if (an + bn == 0) {
			(*reg_it).second.stat.p1 = 0;
		}
		else {
			(*reg_it).second.stat.p1 = double(an) / double(an + bn);
		}

		int ncol, nrow;
    	double emin, expect, percnt, pre, prt;
    	double table[4];
    	ncol = 2;
    	nrow = 2;
    	expect = 0;
    	percnt = 1;
    	emin = 80;

		//double pvalue;
    	table[0] = an;
    	table[1] = cn;
    	table[2] = bn;
    	table[3] = dn;
    	//table[4] = cn;
    	//table[5] = fn;

    	fexact_(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre);
		//pvalue = fisher_ext(K, x, N, s);
		(*reg_it).second.stat.pvalue = pre;
	}
}

void network::run_stat2() {
	map<string, regulon>::iterator reg_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		int an, bn, cn, dn, en, fn, gn, hn, in;
		an = (*reg_it).second.stat.an;
		bn = (*reg_it).second.stat.bn;
		cn = (*reg_it).second.stat.cn;

		dn = (*reg_it).second.stat.dn;
		en = (*reg_it).second.stat.en;
		fn = (*reg_it).second.stat.fn;

		if (an + bn + cn + dn + en + fn == 0) {
			(*reg_it).second.stat.p0 = 0;
		}
		else {
			(*reg_it).second.stat.p0 = double(an + cn + en) / double(an + bn + cn + dn + en + fn);
		}

		if (an + cn + en == 0) {
			(*reg_it).second.stat.p1 = 0;
		}
		else {
			(*reg_it).second.stat.p1 = double(an) / double(an + +bn + cn + dn) * 2;
			(*reg_it).second.stat.p2 = double(cn) / double(an + +bn + cn + dn) * 2;
		}

		int ncol, nrow;
    	double emin, expect, percnt, pre, prt;
    	double table[6];
    	ncol = 2;
    	nrow = 3;
    	expect = 0;
    	percnt = 1;
    	emin = 80;

		//double pvalue;
    	table[0] = an;
    	table[1] = cn;
    	table[2] = en;

    	table[3] = bn;
    	table[4] = dn;
    	table[5] = fn;

    	fexact_(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre);
		//pvalue = fisher_ext(K, x, N, s);
		(*reg_it).second.stat.pvalue = pre;
	}
}

void network::run_stat4_2() {
	map<string, regulon>::iterator reg_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		int an, bn, cn, dn, en, fn, gn, hn, in, jn, kn, ln, mn, qn, rn;
		an = (*reg_it).second.stat.an;
		bn = (*reg_it).second.stat.bn;
		cn = (*reg_it).second.stat.cn;

		dn = (*reg_it).second.stat.dn;
		en = (*reg_it).second.stat.en;
		fn = (*reg_it).second.stat.fn;

		gn = (*reg_it).second.stat.gn;
		hn = (*reg_it).second.stat.hn;
		in = (*reg_it).second.stat.in;

		jn = (*reg_it).second.stat.jn;
		kn = (*reg_it).second.stat.kn;
		ln = (*reg_it).second.stat.ln;

		mn = (*reg_it).second.stat.mn;
		qn = (*reg_it).second.stat.qn;
		rn = (*reg_it).second.stat.rn;

		/*
		if (an + bn + cn + dn + en + fn + gn + hn + in == 0) {
			(*reg_it).second.stat.p0 = 0;
		}
		else {
			(*reg_it).second.stat.p0 = double(an + bn + dn + en + gn + hn) / double(an + bn + cn + dn + en + fn + gn + hn + in);
		}

		if (an + bn + cn + dn + en + fn == 0) {
			(*reg_it).second.stat.p1 = 0;
		}
		else {
			(*reg_it).second.stat.p1 = double(an + en) / double(an + bn + cn + dn + en + fn) * 2;
			(*reg_it).second.stat.p2 = double(bn + dn) / double(an + bn + cn + dn + en + fn) * 2;
		}
        */

		int ncol, nrow;
    	double emin, expect, percnt, pre, prt;
    	double table[9];
    	ncol = 3;
    	nrow = 3;
    	expect = 0;
    	percnt = 1;
    	emin = 80;

		//double pvalue;
    	table[0] = an;
    	table[1] = dn;
    	table[2] = gn;

    	table[3] = bn;
    	table[4] = en;
    	table[5] = hn;

    	table[6] = cn;
    	table[7] = fn;
    	table[8] = in;

    	fexact_(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre);
		//pvalue = fisher_ext(K, x, N, s);
		(*reg_it).second.stat.p0 = pre;

        double table_1[4];
    	ncol = 2;
    	nrow = 2;
    	expect = 0;
    	percnt = 1;
    	emin = 80;

		//double pvalue;
    	table_1[0] = an + en;
    	table_1[1] = jn + qn;

    	table_1[2] = bn + cn + dn + fn;
    	table_1[3] = kn + ln + mn + rn;

    	fexact_(&nrow, &ncol, table_1, &nrow, &expect, &percnt, &emin, &prt, &pre);
		//pvalue = fisher_ext(K, x, N, s);
		(*reg_it).second.stat.p1 = pre;

        double table_2[4];
    	ncol = 2;
    	nrow = 2;
    	expect = 0;
    	percnt = 1;
    	emin = 80;

		//double pvalue;
    	table_2[0] = bn + dn;
    	table_2[1] = kn + mn;

    	table_2[2] = an + cn + en + fn;
    	table_2[3] = jn + ln + qn + rn;

    	fexact_(&nrow, &ncol, table_2, &nrow, &expect, &percnt, &emin, &prt, &pre);
		//pvalue = fisher_ext(K, x, N, s);
		(*reg_it).second.stat.p2 = pre;
	}
}

void network::run_stat4() {
	map<string, regulon>::iterator reg_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		int an, bn, cn, dn, en, fn, gn, hn, in;
		an = (*reg_it).second.stat.an;
		bn = (*reg_it).second.stat.bn;
		cn = (*reg_it).second.stat.cn;

		dn = (*reg_it).second.stat.dn;
		en = (*reg_it).second.stat.en;
		fn = (*reg_it).second.stat.fn;

		gn = (*reg_it).second.stat.gn;
		hn = (*reg_it).second.stat.hn;
		in = (*reg_it).second.stat.in;

		if (an + bn + cn + dn + en + fn + gn + hn + in == 0) {
			(*reg_it).second.stat.p0 = 0;
		}
		else {
			(*reg_it).second.stat.p0 = double(an + bn + dn + en + gn + hn) / double(an + bn + cn + dn + en + fn + gn + hn + in);
		}

		if (an + bn + cn + dn + en + fn == 0) {
			(*reg_it).second.stat.p1 = 0;
		}
		else {
			(*reg_it).second.stat.p1 = double(an + en) / double(an + bn + cn + dn + en + fn) * 2;
			(*reg_it).second.stat.p2 = double(bn + dn) / double(an + bn + cn + dn + en + fn) * 2;
		}

		int ncol, nrow;
    	double emin, expect, percnt, pre, prt;
    	double table[9];
    	ncol = 3;
    	nrow = 3;
    	expect = 0;
    	percnt = 1;
    	emin = 80;

		//double pvalue;
    	table[0] = an;
    	table[1] = dn;
    	table[2] = gn;

    	table[3] = bn;
    	table[4] = en;
    	table[5] = hn;

    	table[6] = cn;
    	table[7] = fn;
    	table[8] = in;

    	fexact_(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre);
		//pvalue = fisher_ext(K, x, N, s);
		(*reg_it).second.stat.pvalue = pre;
	}
}

void regulon::cal_dtype() {
	//begin
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
        if ((*tar_it).rel_type == "unknown") {
            (*tar_it).dri_type = "unknown";
		}
		if ((*tar_it).rel_type == "positive") {
            (*tar_it).dri_type = "positive";
		}
		else if ((*tar_it).rel_type == "negative") {
			(*tar_it).dri_type = "negative";
		}
	}
	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_dtype((*trs_it).second.rel_type);
	}
}

void transf::cal_dtype(string p_rtype) {
	vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
    if ((*tar_it).present == false) {
        continue;
    }
    if (p_rtype == "positive") {
        if ((*tar_it).rel_type == "unknown") {
            (*tar_it).dri_type = "unknown";
		}
        else if ((*tar_it).rel_type == "positive") {
            (*tar_it).dri_type = "positive";
		}
        else if ((*tar_it).rel_type == "negative") {
            (*tar_it).dri_type = "negative";
		}
	}
	else if (p_rtype == "negative") {
        if ((*tar_it).rel_type == "unknown") {
            (*tar_it).dri_type = "unknown";
		}
        else if ((*tar_it).rel_type == "positive") {
            (*tar_it).dri_type = "negative";
		}
        else if ((*tar_it).rel_type == "negative") {
            (*tar_it).dri_type = "positive";
		}
	}
	else if (p_rtype == "unknown") {
        (*tar_it).dri_type = "unknown";
	}
	}
}

void network::cal_dtype() {
	map<string, regulon>::iterator reg_it;
	//map<string, transf>::iterator trs_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		(*reg_it).second.cal_dtype();
	}
}

void regulon::cal_cfi(int &cn, int &fn, int &in) {
    //map<string, gene>::iterator gene_it;
    cn = 0; fn = 0; in = 0;
    vector<gene*>::iterator genes_pt_it;
    for (genes_pt_it = _genes.begin(); genes_pt_it != _genes.end(); genes_pt_it++) {
        if ((**genes_pt_it).expresion == "up") {
            cn++;
        }
        if ((**genes_pt_it).expresion == "down") {
            fn++;
        }
        if ((**genes_pt_it).expresion == "nodiff") {
            in++;
        }
    }
}

void transf::cal_beh(int &bn, int &en, int &hn) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "negative") {
            bn++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "negative") {
			en++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "negative") {
			hn++;
		}
	}
}

void regulon::cal_beh(int &bn, int &en, int &hn) {
	//begin
	bn = 0; en = 0; hn = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "negative") {
            bn++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "negative") {
			en++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "negative") {
			hn++;
		}
	}

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_beh(bn, en, hn);
	}
}

void transf::cal_adg(int &an, int &dn, int &gn) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "positive") {
            an++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "positive") {
			dn++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "positive") {
			gn++;
		}
	}
}

void regulon::cal_adg(int &an, int &dn, int &gn) {
	//begin
	an = 0; dn = 0; gn = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "positive") {
            an++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "positive") {
			dn++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "positive") {
			gn++;
		}
	}

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_adg(an, dn, gn);
	}
}

void regulon::cal_cfi2(int &cn, int &fn, int &in, int &ln, int &rn) {
    //map<string, gene>::iterator gene_it;
    cn = 0; fn = 0; in = 0; ln = 0; rn = 0;
    vector<gene*>::iterator genes_pt_it;
    for (genes_pt_it = _genes.begin(); genes_pt_it != _genes.end(); genes_pt_it++) {
        if ((**genes_pt_it).expresion == "up") {
            cn++;
        }
        if ((**genes_pt_it).expresion == "down") {
            fn++;
        }
        if ((**genes_pt_it).expresion == "nodiff") {
            //in++;
            if ((**genes_pt_it).fc > 0.0) {
                ln++;
			}
			else {
                rn++;
			}
        }
    }
    in = ln + rn;
}

void transf::cal_beh2(int &bn, int &en, int &hn, int &kn, int &qn) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "negative") {
            bn++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "negative") {
			en++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "negative") {
			//hn++;
			if ((*(*tar_it).gene_pt).fc > 0.0) {
                kn++;
			}
			else {
                qn++;
			}
		}
	}
	hn = kn + qn;
}

void regulon::cal_beh2(int &bn, int &en, int &hn, int &kn, int &qn) {
	//begin
	bn = 0; en = 0; hn = 0; kn = 0; qn = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "negative") {
            bn++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "negative") {
			en++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "negative") {
			//hn++;
			if ((*(*tar_it).gene_pt).fc > 0.0) {
                kn++;
			}
			else {
                qn++;
			}
		}
	}
    hn = kn + qn;

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_beh2(bn, en, hn, kn, qn);
	}
}

void transf::cal_adg2(int &an, int &dn, int &gn, int &jn, int &mn) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "positive") {
            an++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "positive") {
			dn++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "positive") {
			//gn++;
			if ((*(*tar_it).gene_pt).fc > 0.0) {
                jn++;
			}
			else {
                mn++;
			}
		}
	}
	gn = jn + mn;
}

void regulon::cal_adg2(int &an, int &dn, int &gn, int &jn, int &mn) {
	//begin
	an = 0; dn = 0; gn = 0; jn = 0; mn = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).expresion == "up" && (*tar_it).dri_type == "positive") {
            an++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).expresion == "down" && (*tar_it).dri_type == "positive") {
			dn++;
		}
		else if ((*(*tar_it).gene_pt).expresion == "nodiff" && (*tar_it).dri_type == "positive") {
			if ((*(*tar_it).gene_pt).fc > 0.0) {
                jn++;
			}
			else {
                mn++;
			}
			//gn++;
		}
	}
    gn = jn + mn;

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_adg2(an, dn, gn, jn, mn);
	}
}

void network::cal_abcdefghi() {
	map<string, regulon>::iterator reg_it;
	bionstat *stat_pt;
	//map<string, transf>::iterator trs_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		stat_pt = &(*reg_it).second.stat;
		(*stat_pt).type = 4;
		(*reg_it).second.cal_adg((*stat_pt).an, (*stat_pt).dn, (*stat_pt).gn);
		(*reg_it).second.cal_beh((*stat_pt).bn, (*stat_pt).en, (*stat_pt).hn);
		(*reg_it).second.cal_cfi((*stat_pt).cn, (*stat_pt).fn, (*stat_pt).in);
	}
}

void network::cal_abcdefghi2() {
	map<string, regulon>::iterator reg_it;
	bionstat *stat_pt;
	//map<string, transf>::iterator trs_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		stat_pt = &(*reg_it).second.stat;
		(*stat_pt).type = 4;
		(*reg_it).second.cal_adg2((*stat_pt).an, (*stat_pt).dn, (*stat_pt).gn, (*stat_pt).jn, (*stat_pt).mn);
		(*reg_it).second.cal_beh2((*stat_pt).bn, (*stat_pt).en, (*stat_pt).hn, (*stat_pt).kn, (*stat_pt).qn);
		(*reg_it).second.cal_cfi2((*stat_pt).cn, (*stat_pt).fn, (*stat_pt).in, (*stat_pt).ln, (*stat_pt).rn);
	}
}

void network::run_stat3() {
	map<string, regulon>::iterator reg_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		int an, bn, cn, dn, en, fn, gn, hn, in;
		an = (*reg_it).second.stat.an;
		bn = (*reg_it).second.stat.bn;
		cn = (*reg_it).second.stat.cn;

		dn = (*reg_it).second.stat.dn;
		en = (*reg_it).second.stat.en;
		fn = (*reg_it).second.stat.fn;

		if (an + bn + cn + dn + en + fn == 0) {
			(*reg_it).second.stat.p0 = 0;
		}
		else {
			(*reg_it).second.stat.p0 = double(an + bn + dn + en) / double(an + bn + cn + dn + en + fn);
		}

		if (an + bn + cn == 0) {
			(*reg_it).second.stat.p1 = 0;
		}
		else {
			(*reg_it).second.stat.p1 = double(an) / double(an + bn + cn) * 2;
			(*reg_it).second.stat.p2 = double(bn) / double(an + bn + cn) * 2;
		}

		int ncol, nrow;
    	double emin, expect, percnt, pre, prt;
    	double table[6];
    	ncol = 3;
    	nrow = 2;
    	expect = 0;
    	percnt = 1;
    	emin = 80;

		//double pvalue;
    	table[0] = an;
    	table[1] = dn;

    	table[2] = bn;
    	table[3] = en;

    	table[4] = cn;
    	table[5] = fn;

    	fexact_(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre);
		//pvalue = fisher_ext(K, x, N, s);
		(*reg_it).second.stat.pvalue = pre;
	}
}

void regulon::cal_cf(int &cn, int &fn) {
    //map<string, gene>::iterator gene_it;
    cn = 0; fn = 0;
    vector<gene*>::iterator genes_pt_it;
    for (genes_pt_it = _genes.begin(); genes_pt_it != _genes.end(); genes_pt_it++) {
        if ((**genes_pt_it).cls == "cls") {
            cn++;
        }
        if ((**genes_pt_it).cls == "---") {
            fn++;
        }
    }
}

void transf::cal_be(int &bn, int &en) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).cls == "cls" && (*tar_it).dri_type == "negative") {
            bn++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).cls == "---" && (*tar_it).dri_type == "negative") {
			en++;
		}
	}
}

void regulon::cal_be(int &bn, int &en) {
	//begin
	bn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).cls == "cls" && (*tar_it).dri_type == "negative") {
            bn++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).cls == "---" && (*tar_it).dri_type == "negative") {
			en++;
		}
	}

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_be(bn, en);
	}
}

void transf::cal_ad(int &an, int &dn) {
	//begin
	//an = 0; bn = 0; dn = 0; en = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).cls == "cls" && (*tar_it).dri_type == "positive") {
            an++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).cls == "---" && (*tar_it).dri_type == "positive") {
			dn++;
		}
	}
}

void regulon::cal_ad(int &an, int &dn) {
	//begin
	an = 0; dn = 0;
    vector<target>::iterator tar_it;
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
        if ((*tar_it).present == false) {
            continue;
        }
		if ((*(*tar_it).gene_pt).cls == "cls" && (*tar_it).dri_type == "positive") {
            an++;
            //cout<<(*tar_it).ref_id<<endl;
		}
		else if ((*(*tar_it).gene_pt).cls == "---" && (*tar_it).dri_type == "positive") {
			dn++;
		}
	}

	map<string, transf>::iterator trs_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
        (*trs_it).second.cal_ad(an, dn);
	}
}

void network::cal_abcdef() {
	map<string, regulon>::iterator reg_it;
	bionstat *stat_pt;
	//map<string, transf>::iterator trs_it;
	for (reg_it = regulons.begin(); reg_it != regulons.end(); reg_it++) {
		stat_pt = &(*reg_it).second.stat;
		(*stat_pt).type = 3;
		(*reg_it).second.cal_ad((*stat_pt).an, (*stat_pt).dn);
		(*reg_it).second.cal_be((*stat_pt).bn, (*stat_pt).en);
		(*reg_it).second.cal_cf((*stat_pt).cn, (*stat_pt).fn);
	}
}

void gerea::run_stat() {
	if (data_it.data_type == 1 && network_it.db_type == 1) {
		network_it.cal_abcd1();
		network_it.run_stat1();
		network_it.BH_correction();
	}
	else if (data_it.data_type == 2 && network_it.db_type == 1) {
		network_it.cal_abcdef2();
		network_it.run_stat2();
		network_it.BH_correction();
	}
	else if (data_it.data_type == 1 && network_it.db_type == 2) {
		network_it.cal_dtype();
		network_it.cal_abcdef();
		network_it.run_stat3();
		network_it.BH_correction();
	}
	else if (data_it.data_type == 2 && network_it.db_type == 2) {
		network_it.cal_dtype();
		//network_it.cal_abcdefghi();
		network_it.cal_abcdefghi2();
		network_it.run_stat4_2();
		network_it.BH_correction_2();
	}
	else if (data_it.data_type == 2 && network_it.db_type == 3) {
		network_it.cal_dtype();
		//network_it.cal_abcdefghi();
		network_it.cal_abcdefghi2();
		network_it.run_stat4_2();
		network_it.BH_correction_2();
	}
	/*
    network_it.cal_dtype();
    //cout<<"qvalue_threshold: "<<data_it.qvalue_threshold<<endl;
    network_it.cal_abcdef(data_it.qvalue_threshold);
    network_it.run_stat2();
    network_it.BH_correction();
    */
}

void gerea::sortby_fdr(map<string, double> &fdr_map, list<string> &sortby) {
	map<string, double>::iterator fdr_rt;
	list<string>::iterator sort_rt;
	for (fdr_rt = fdr_map.begin(); fdr_rt != fdr_map.end(); fdr_rt++) {
		bool found = false;
		for (sort_rt = sortby.begin(); sort_rt != sortby.end(); sort_rt++) {
			if ((*fdr_rt).second < fdr_map[*sort_rt]) {
				found = true;
				sortby.insert(sort_rt, (*fdr_rt).first);
				break;
			}
		}
		if (found == false) {
			sortby.push_back((*fdr_rt).first);
		}
	}
}

void gerea::Print_result1(string output_dir) {
	string file_name = output_dir + "/" + session_id + ".ger.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}

	writer << "#result_type=1\n";
	//writer << "regulon\ta\tb\tc\td\tp0\tp1\tpvalue\tfdr\n";

	map<string, regulon>::iterator reg_it;
	bionstat *bstat;

	map<string, double> fdr_map;
	list<string> sortby;
	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		bstat = &((*reg_it).second.stat);
		if ((*bstat).get_tar_n() > 0 && (*bstat).p_nequal0()) {
			fdr_map[(*reg_it).first] = (*reg_it).second.stat.fdr_BH;
		}
	}
	sortby_fdr(fdr_map, sortby);
	list<string>::iterator sort_rt;
	for (sort_rt = sortby.begin(); sort_rt != sortby.end(); sort_rt++) {
		bstat = &(network_it.regulons[*sort_rt].stat);
	//for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		//bstat = &((*reg_it).second.stat);
		if ((*bstat).get_tar_n() > 0) {
			//writer << (*reg_it).first;
			writer << (*sort_rt);
            writer << "\t" << (*bstat).an;
            writer << "\t" << (*bstat).bn;
            writer << "\t" << (*bstat).cn;
            writer << "\t" << (*bstat).dn;
            writer << "\t" << (*bstat).p0;
            writer << "\t" << (*bstat).p1;
            writer << "\t" << (*bstat).pvalue;
            writer << "\t" << (*bstat).fdr_BH;
            writer << endl;
		}
	}
	writer.close();
}

void gerea::Print_result3(string output_dir) {
	string file_name = output_dir + "/" + session_id + ".ger.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}

	writer << "#result_type=3\n";
	//writer << "regulon\ta\tb\tc\td\te\tf\tp0\tp1\tp2\tpvalue\tfdr\n";

	map<string, regulon>::iterator reg_it;
	bionstat *bstat;

	map<string, double> fdr_map;
	list<string> sortby;
	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		bstat = &((*reg_it).second.stat);
		if ((*bstat).get_tar_n() > 0 && (*bstat).p_nequal0()) {
			fdr_map[(*reg_it).first] = (*reg_it).second.stat.fdr_BH;
		}
	}
	sortby_fdr(fdr_map, sortby);
	list<string>::iterator sort_rt;
	for (sort_rt = sortby.begin(); sort_rt != sortby.end(); sort_rt++) {
		bstat = &(network_it.regulons[*sort_rt].stat);

	//for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		//bstat = &((*reg_it).second.stat);
		if ((*bstat).get_tar_n() > 0) {
			writer << (*sort_rt);
            writer << "\t" << (*bstat).an;
            writer << "\t" << (*bstat).bn;
            writer << "\t" << (*bstat).cn;
            writer << "\t" << (*bstat).dn;
            writer << "\t" << (*bstat).en;
            writer << "\t" << (*bstat).fn;
            writer << "\t" << (*bstat).p0;
            writer << "\t" << (*bstat).p1;
            writer << "\t" << (*bstat).p2;
            writer << "\t" << (*bstat).pvalue;
            writer << "\t" << (*bstat).fdr_BH;
            writer << endl;
		}
	}
	writer.close();
}

void gerea::Print_result4_2(string output_dir) {
	string file_name = output_dir + "/" + session_id + ".ger.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}

	writer << "#result_type=3\n";
	//writer << "regulon\ta\tb\tc\td\te\tf\tg\th\ti\tp0\tp1\tp2\n";

	map<string, regulon>::iterator reg_it;
	bionstat *bstat;

	map<string, double> fdr_map;
	list<string> sortby;
	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		bstat = &((*reg_it).second.stat);
		//if ((*bstat).get_tar_n() > 0 && (*bstat).p_nequal0()) {
		if ((*bstat).get_tar_n() >= target_n) {
            if ((*reg_it).second.stat.fdr_BH_1 <= (*reg_it).second.stat.fdr_BH_2) {
                fdr_map[(*reg_it).first] = (*reg_it).second.stat.fdr_BH_1;
            }
            else {
                fdr_map[(*reg_it).first] = (*reg_it).second.stat.fdr_BH_2;
            }
		}
	}
	sortby_fdr(fdr_map, sortby);
	list<string>::iterator sort_rt;
	for (sort_rt = sortby.begin(); sort_rt != sortby.end(); sort_rt++) {
		bstat = &(network_it.regulons[*sort_rt].stat);

	//for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		//bstat = &((*reg_it).second.stat);
		if ((*bstat).get_tar_n() > target_n) {
			writer << (*sort_rt);
            writer << "\t" << (*bstat).an;
            writer << "\t" << (*bstat).bn;
            writer << "\t" << (*bstat).cn;
            writer << "\t" << (*bstat).dn;
            writer << "\t" << (*bstat).en;
            writer << "\t" << (*bstat).fn;
            writer << "\t" << (*bstat).gn;
            writer << "\t" << (*bstat).hn;
            writer << "\t" << (*bstat).in;

            /*
            writer << "\t" << (*bstat).jn;
            writer << "\t" << (*bstat).kn;
            writer << "\t" << (*bstat).ln;
            writer << "\t" << (*bstat).mn;
            writer << "\t" << (*bstat).qn;
            writer << "\t" << (*bstat).rn;
            */

            writer << "\t" << (*bstat).fdr_BH_0;
            writer << "\t" << (*bstat).fdr_BH_1;
            writer << "\t" << (*bstat).fdr_BH_2;
            writer << endl;
		}
	}
	writer.close();
}

void gerea::Print_result4(string output_dir) {
	string file_name = output_dir + "/" + session_id + ".ger.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}

	writer << "#result_type=4\n";
	//writer << "regulon\ta\tb\tc\td\te\tf\tg\th\ti\tp0\tp1\tp2\tpvalue\tfdr\n";

	map<string, regulon>::iterator reg_it;
	bionstat *bstat;

	map<string, double> fdr_map;
	list<string> sortby;
	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		bstat = &((*reg_it).second.stat);
		if ((*bstat).get_tar_n() > 0 && (*bstat).p_nequal0()) {
			fdr_map[(*reg_it).first] = (*reg_it).second.stat.fdr_BH;
		}
	}
	sortby_fdr(fdr_map, sortby);
	list<string>::iterator sort_rt;
	for (sort_rt = sortby.begin(); sort_rt != sortby.end(); sort_rt++) {
		bstat = &(network_it.regulons[*sort_rt].stat);

	//for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		//bstat = &((*reg_it).second.stat);
		if ((*bstat).get_tar_n() > 0) {
			writer << (*sort_rt);
            writer << "\t" << (*bstat).an;
            writer << "\t" << (*bstat).bn;
            writer << "\t" << (*bstat).cn;
            writer << "\t" << (*bstat).dn;
            writer << "\t" << (*bstat).en;
            writer << "\t" << (*bstat).fn;
            writer << "\t" << (*bstat).gn;
            writer << "\t" << (*bstat).hn;
            writer << "\t" << (*bstat).in;
            writer << "\t" << (*bstat).p0;
            writer << "\t" << (*bstat).p1;
            writer << "\t" << (*bstat).p2;
            writer << "\t" << (*bstat).pvalue;
            writer << "\t" << (*bstat).fdr_BH;
            writer << endl;
		}
	}
	writer.close();
}

void gerea::Print_result(string output_dir) {
    //cout<<"data_type: "<<data_it.data_type<<endl;
    //cout<<"network_type: "<<network_it.db_type<<endl;
	if (data_it.data_type == 1 && network_it.db_type == 1) {
		Print_result1(output_dir);
	}
	else if (data_it.data_type == 2 && network_it.db_type == 1) {
		Print_result3(output_dir);
	}
	else if (data_it.data_type == 1 && network_it.db_type == 2) {
		Print_result3(output_dir);
	}
	else if (data_it.data_type == 2 && network_it.db_type == 2) {
        //cout<<"OK6"<<endl;
		Print_result4_2(output_dir);
	}
	else if (data_it.data_type == 2 && network_it.db_type == 3) {
        //cout<<"OK6"<<endl;
		Print_result4_2(output_dir);
	}
}

/*
void gerea::Print_result() {
	string file_name = session_id + ".ger.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}
	writer << "regulon\ta\tb\tc\td\te\tf\tp0\tp1\tp2\tpvalue\tfdr\n";

	map<string, regulon>::iterator reg_it;
	bionstat *bstat;
	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		bstat = &((*reg_it).second.stat);
		//if ((*bstat).get_tar_n() > 0) {
			writer << (*reg_it).first;
            writer << "\t" << (*bstat).an;
            writer << "\t" << (*bstat).bn;
            writer << "\t" << (*bstat).cn;
            writer << "\t" << (*bstat).dn;
            writer << "\t" << (*bstat).en;
            writer << "\t" << (*bstat).fn;
            writer << "\t" << (*bstat).p0;
            writer << "\t" << (*bstat).p1;
            writer << "\t" << (*bstat).p2;
            writer << "\t" << (*bstat).pvalue;
            writer << "\t" << (*bstat).fdr_BH;
            writer << endl;
		//}
	}
	writer.close();
}
*/

void transf::encode_network2(vector<string> &tmp_strs) {
	//tids.clear();
	tmp_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->gene_pt->expresion;
			tmp_strs.push_back(ecd_str);
		}
	}
}

void regulon::encode_network2(vector<string> &ecd_strs) {
	//tids.clear();
	ecd_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->gene_pt->expresion;
			ecd_strs.push_back(ecd_str);
		}
	}

	//set<string> tmp_tids;
	vector<string> tmp_strs;
	map<string, transf>::iterator trs_it;
	vector<string>::iterator str_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
		(*trs_it).second.encode_network2(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			ecd_str =  (*trs_it).second.rel_type + "\t" + (*trs_it).first + "\t" + (*str_it);
			ecd_strs.push_back(ecd_str);
		}
		//tids.insert(tmp_tids.begin(), tmp_tids.end());
	}
}

void gerea::Print_details(string output_dir) {
	if (data_it.data_type == 1 && network_it.db_type == 1) {
		Print_details1(output_dir);
	}
	else if (data_it.data_type == 2 && network_it.db_type == 1) {
		Print_details2(output_dir);
	}
	else if (data_it.data_type == 1 && network_it.db_type == 2) {
		Print_details3(output_dir);
	}
	else if (data_it.data_type == 2 && network_it.db_type == 2) {
		Print_details4(output_dir);
	}
	else if (data_it.data_type == 2 && network_it.db_type == 3) {
		Print_details4(output_dir);
	}
}

void gerea::Print_details2(string output_dir) {
	string file_name = session_id + ".links.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}
	//writer << "pathw_name\tdata_N\tgene_s\tdata_K\tgene_x\tp0\tp1\tpvalue\tfdr\tgene_name\n";

	map<string, regulon>::iterator reg_it;
	//bionstat *bstat;
	vector<string> tmp_strs;
	vector<string>::iterator str_it;

	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		tmp_strs.clear();
		(*reg_it).second.encode_network2(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			writer << (*reg_it).first;
			writer << "\t" << (*str_it);
			writer << endl;
		}
//		bstat = &((*reg_it).second.stat);
//		if ((*bstat).N_diff == 0) {
//			continue;
//		}
//		writer << (*reg_it).first;
//		writer << "\t" << (*bstat).data_N;
//		writer << "\t" << (*bstat).gene_s;
//		writer << "\t" << (*bstat).data_K;
//		writer << "\t" << (*bstat).gene_x;
//		writer << "\t" << (*bstat).p0;
//		writer << "\t" << (*bstat).p1;
//		writer << "\t" << (*bstat).pvalue;
//		writer << "\t" << (*bstat).fdr_BH;
//		writer << "\t" << "Tar_Str_Diff(pid)";
//		writer << endl;

	}
	writer.close();
}

void transf::encode_network1(vector<string> &tmp_strs) {
	//tids.clear();
	tmp_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->gene_pt->cls;
			tmp_strs.push_back(ecd_str);
		}
	}
}

void regulon::encode_network1(vector<string> &ecd_strs) {
	//tids.clear();
	ecd_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->gene_pt->cls;
			ecd_strs.push_back(ecd_str);
		}
	}

	//set<string> tmp_tids;
	vector<string> tmp_strs;
	map<string, transf>::iterator trs_it;
	vector<string>::iterator str_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
		(*trs_it).second.encode_network1(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			ecd_str =  (*trs_it).second.rel_type + "\t" + (*trs_it).first + "\t" + (*str_it);
			ecd_strs.push_back(ecd_str);
		}
		//tids.insert(tmp_tids.begin(), tmp_tids.end());
	}
}

void gerea::Print_details1(string output_dir) {
	string file_name = output_dir + "/" + session_id + ".links.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}
	//writer << "pathw_name\tdata_N\tgene_s\tdata_K\tgene_x\tp0\tp1\tpvalue\tfdr\tgene_name\n";

	map<string, regulon>::iterator reg_it;
	//bionstat *bstat;
	vector<string> tmp_strs;
	vector<string>::iterator str_it;

	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		tmp_strs.clear();
		(*reg_it).second.encode_network1(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			writer << (*reg_it).first;
			writer << "\t" << (*str_it);
			writer << endl;
		}
//		bstat = &((*reg_it).second.stat);
//		if ((*bstat).N_diff == 0) {
//			continue;
//		}
//		writer << (*reg_it).first;
//		writer << "\t" << (*bstat).data_N;
//		writer << "\t" << (*bstat).gene_s;
//		writer << "\t" << (*bstat).data_K;
//		writer << "\t" << (*bstat).gene_x;
//		writer << "\t" << (*bstat).p0;
//		writer << "\t" << (*bstat).p1;
//		writer << "\t" << (*bstat).pvalue;
//		writer << "\t" << (*bstat).fdr_BH;
//		writer << "\t" << "Tar_Str_Diff(pid)";
//		writer << endl;

	}
	writer.close();
}

void transf::encode_network3(vector<string> &tmp_strs) {
	//tids.clear();
	tmp_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->dri_type + "\t" + tar_it->gene_pt->cls;
			tmp_strs.push_back(ecd_str);
		}
	}
}

void regulon::encode_network3(vector<string> &ecd_strs) {
	//tids.clear();
	ecd_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->dri_type + "\t" + tar_it->gene_pt->cls;
			ecd_strs.push_back(ecd_str);
		}
	}

	//set<string> tmp_tids;
	vector<string> tmp_strs;
	map<string, transf>::iterator trs_it;
	vector<string>::iterator str_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
		(*trs_it).second.encode_network3(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			ecd_str =  (*trs_it).second.rel_type + "\t" + (*trs_it).first + "\t" + (*str_it);
			ecd_strs.push_back(ecd_str);
		}
		//tids.insert(tmp_tids.begin(), tmp_tids.end());
	}
}

void gerea::Print_details3(string output_dir) {
	string file_name = output_dir + "/" + session_id + ".links.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}
	//writer << "pathw_name\tdata_N\tgene_s\tdata_K\tgene_x\tp0\tp1\tpvalue\tfdr\tgene_name\n";

	map<string, regulon>::iterator reg_it;
	//bionstat *bstat;
	vector<string> tmp_strs;
	vector<string>::iterator str_it;

	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		tmp_strs.clear();
		(*reg_it).second.encode_network3(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			writer << (*reg_it).first;
			writer << "\t" << (*str_it);
			writer << endl;
		}
//		bstat = &((*reg_it).second.stat);
//		if ((*bstat).N_diff == 0) {
//			continue;
//		}
//		writer << (*reg_it).first;
//		writer << "\t" << (*bstat).data_N;
//		writer << "\t" << (*bstat).gene_s;
//		writer << "\t" << (*bstat).data_K;
//		writer << "\t" << (*bstat).gene_x;
//		writer << "\t" << (*bstat).p0;
//		writer << "\t" << (*bstat).p1;
//		writer << "\t" << (*bstat).pvalue;
//		writer << "\t" << (*bstat).fdr_BH;
//		writer << "\t" << "Tar_Str_Diff(pid)";
//		writer << endl;

	}
	writer.close();
}

void transf::encode_network4(vector<string> &tmp_strs) {
	//tids.clear();
	tmp_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->dri_type + "\t" + tar_it->gene_pt->expresion;
			tmp_strs.push_back(ecd_str);
		}
	}
}

void regulon::encode_network4(vector<string> &ecd_strs) {
	//tids.clear();
	ecd_strs.clear();
	vector<target>::iterator tar_it;
	string ecd_str = "";
	for (tar_it = targets.begin(); tar_it != targets.end(); tar_it++) {
		if ((*tar_it).present == true) {
			ecd_str = tar_it->rel_type + "\t" + tar_it->ref_id + "\t" + tar_it->dri_type + "\t" + tar_it->gene_pt->expresion;
			ecd_strs.push_back(ecd_str);
		}
	}

	//set<string> tmp_tids;
	vector<string> tmp_strs;
	map<string, transf>::iterator trs_it;
	vector<string>::iterator str_it;
	for (trs_it = transfs.begin(); trs_it != transfs.end(); trs_it++) {
		(*trs_it).second.encode_network4(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			ecd_str =  (*trs_it).second.rel_type + "\t" + (*trs_it).first + "\t" + (*str_it);
			ecd_strs.push_back(ecd_str);
		}
		//tids.insert(tmp_tids.begin(), tmp_tids.end());
	}
}

void gerea::Print_details4(string output_dir) {
	string file_name = output_dir + "/" + session_id + ".links.txt";
	ofstream writer;
	writer.open(file_name.c_str(), ios::out);
	if (writer.fail()) {
		cout << "cannot open file " << file_name << "\n";
		exit(1);
	}
	//writer << "pathw_name\tdata_N\tgene_s\tdata_K\tgene_x\tp0\tp1\tpvalue\tfdr\tgene_name\n";

	map<string, regulon>::iterator reg_it;
	//bionstat *bstat;
	vector<string> tmp_strs;
	vector<string>::iterator str_it;

	for (reg_it = network_it.regulons.begin(); reg_it != network_it.regulons.end(); reg_it++) {
		tmp_strs.clear();
		(*reg_it).second.encode_network4(tmp_strs);
		for (str_it = tmp_strs.begin(); str_it != tmp_strs.end(); str_it++) {
			writer << (*reg_it).first;
			writer << "\t" << (*str_it);
			writer << endl;
		}
//		bstat = &((*reg_it).second.stat);
//		if ((*bstat).N_diff == 0) {
//			continue;
//		}
//		writer << (*reg_it).first;
//		writer << "\t" << (*bstat).data_N;
//		writer << "\t" << (*bstat).gene_s;
//		writer << "\t" << (*bstat).data_K;
//		writer << "\t" << (*bstat).gene_x;
//		writer << "\t" << (*bstat).p0;
//		writer << "\t" << (*bstat).p1;
//		writer << "\t" << (*bstat).pvalue;
//		writer << "\t" << (*bstat).fdr_BH;
//		writer << "\t" << "Tar_Str_Diff(pid)";
//		writer << endl;

	}
	writer.close();
}
