//============================================================================
// Name        : gerea-new.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <getopt.h>
#include "gerea.h"
#include <map>
#include <string>

#include <iostream>
using namespace std;

char const short_options[] = "hn:e:q:f:c:s:t:d:";
struct option long_options[] =
{
    {"help", 0, NULL, 'h'},
    {"network", 1, NULL, 'n'},
	{"expression", 1, NULL, 'e'},
	{"qvalue", 1, NULL, 'q'},
	{"foldchange", 1, NULL, 'f'},
	{"class", 1, NULL, 'c'},
    {"session", 1, NULL, 's'},
    {"target_n", 1, NULL, 't'},
    {"dir", 1, NULL, 'd'},
    {0, 0, 0, 0}

};

int main(int argc, char *argv[]) {
	if (argc == 1)
	{
		printf("gerea -n network -e expression -s session -d dir [-q qvalue] [-f foldchange] [-c class]\n");
		exit(1);

	}

	map<char, int> options;
	options['n'] = 0;
	options['e'] = 0;
	options['s'] = 0;
	options['d'] = 0;

	string net_file, exp_file, session_name, output_dir;
	float q_threshold, fc_threshold;
	int target_n;
	string class_str;
	q_threshold = 0.05;
	fc_threshold = 1.5;
	class_str = "cls";
	target_n = 10;

	int data_seto = 0;

	int c;
	while((c=getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
	{
		switch (c)
		{
		case 'h':
			printf("gerea -n network -e expression -s session -d dir [-q qvalue] [-f foldchange] [-t target_n] [-c class]\n");
			exit(0);
			break;
		case 'n':
			net_file = optarg;
			options['n'] = 1;
			break;
		case 'e':
			exp_file = optarg;
			options['e'] = 1;
			break;
		case 'q':
			//cout<<optarg<<endl;
			q_threshold = atof(optarg);
			if (data_seto == 0) {
				data_seto = 2;
			}
			else if (data_seto == 1) {
				cout<<"conflict with -c option!"<<endl;
				exit(1);
			}
			break;
		case 'f':
			//cout<<optarg<<endl;
			fc_threshold = atof(optarg);
			if (data_seto == 0) {
				data_seto = 2;
			}
			else if (data_seto == 1) {
				cout<<"conflict with -c option!"<<endl;
				exit(1);
			}
			break;
        case 't':
			//cout<<optarg<<endl;
			target_n = atoi(optarg);
			break;
		case 'c':
			class_str = optarg;
			if (data_seto == 0) {
				data_seto = 1;
			}
			else if (data_seto == 2) {
				cout<<"conflict with -q and -f option!"<<endl;
				exit(1);
			}
			break;
		case 's':
			session_name = optarg;
			options['s'] = 1;
			break;
		case 'd':
			output_dir = optarg;
			options['d'] = 1;
			break;
		}
	}

	if (options['n'] == 0) {
		cout<<"-n is required!"<<endl;
		exit(1);
	}
	if (options['e'] == 0) {
		cout<<"-e is required!"<<endl;
		exit(1);
	}
	if (options['s'] == 0) {
		cout<<"-s is required!"<<endl;
		exit(1);
	}
	if (options['d'] == 0) {
		cout<<"-d is required!"<<endl;
		exit(1);
	}

	gerea gerea_it(session_name);
	//cout<<q_threshold<<endl;
	//cout<<fc_threshold<<endl;
    gerea_it.target_n = target_n;

	gerea_it.data_it.data_seto = data_seto;
	gerea_it.data_it.q_threshold = q_threshold;
	gerea_it.data_it.fc_threshold = fc_threshold;
	//gerea_it.data_it.fc_threshold = fc_threshold;
	gerea_it.data_it.class_str = class_str;

	//cout<<gerea_it.data_it.q_threshold<<"\t"<<"OK0"<<endl;

	gerea_it.reading_network(net_file);
	//cout<<"OK1"<<endl;

	gerea_it.reading_data(exp_file);
	//cout<<"OK2"<<endl;

	gerea_it.loading_data2network();
	//cout<<"OK3"<<endl;

	gerea_it.run_stat();
	//cout<<"OK4"<<endl;

	gerea_it.Print_result(output_dir);
	//cout<<"OK5"<<endl;

	gerea_it.Print_details(output_dir);

	return 0;

}
