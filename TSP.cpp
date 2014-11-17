#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <getopt.h>
#include <cstring>
#include <vector>
#include <list>
#include <cmath>
#include <stack>
#include <limits>
#include <sstream>
#include <deque>
#include "algorithms.h"
using namespace std;

void printHelp();

int main(int argc, char * argv[]) {
	std:: ios_base::sync_with_stdio(false);

	cout << setprecision(2);
	cout << fixed;

	struct option longOpts[] = {
		{"help", 		no_argument, 		NULL, 'h'},
		{NULL,			0,					NULL,  0 }
	};

	opterr = 1;
	int opt = 0, index = 0;
	while ((opt = getopt_long(argc,argv,"h:",longOpts,&index)) != -1) {
		switch (opt) {
			case 'h':
				printHelp();
				return 0;
			case '?':
			default:
				printHelp();
				exit(1);
		}
	}
		
	int num_locations, x, y;
	double total = 0;
	cin >> num_locations;
	deque<loc> optpath;
	list<loc> path;
	
	//Reads in points and create an upper bound and the deque for the opt calculation
	for(int i = 0; i < num_locations; ++i){
		cin >> x >> y;
		loc in(x, y, i);	
		UBound(path, in, i, total);
		optpath.push_back(in);
	}
	//Runs program
	Optimal(optpath, total);	
	return 0;
	
	
}

static const char helpText[] = "Input number of locations, then coordiantes on an XY plane, will output most efficient hamiltonian cycle through all of the nodes";
void printHelp() {
	cout << helpText << flush;
}
