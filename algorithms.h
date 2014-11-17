#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include <list>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <getopt.h>
#include <cstring>
#include <vector>
#include <cmath>
#include <limits>
#include <deque>
#include <stack>
#include <sstream>
using namespace std;

//Location struct
struct loc{
	int x;
	int y;
	int input;
	double distance;
	int parent;
	loc(int x_in, int y_in, int in) : x(x_in), y(y_in), input(in), distance(numeric_limits<double>::infinity()), parent(-1) {}
};


//Distance function
inline double distance(const loc &a, const loc &b){
	double dx = pow(a.x - b.x, 2);
	double dy = pow(a.y - b.y, 2);
	return sqrt(dx+dy);
}


//Finds the next closest node to add to the MST
loc least_finder(vector<loc> &in, int index, double &totalWeight){
	double least = in[index].distance;	
	int out = index;
	for(unsigned int i = index; i < in.size(); i++){
		if(in[i].distance < least){
			least = in[i].distance;
			out = i;
		}
	}
	if(index != 0) totalWeight += least;
	swap(in[index], in[out]);
	return in[index];	
}



//Creates an upper bound
void UBound(list<loc> &path, loc &in, int i, double &total){
	if(i == 0){
		path.push_back(in);
		path.push_back(in);
	}else if(i == 1){
		path.insert(++path.begin(), in);
		total += 2*distance(*path.begin(), in);
	}else{
		//Insertion method, finds best fit and adds it to the corresponding location
		double min = numeric_limits<double>::infinity();
		auto start = path.begin();
		auto end = ++path.begin();
		auto placer = end;
		double temp = 0;
		for(unsigned int i = 0; i < path.size() - 1; ++i){
			temp = -distance(*start, *end) + distance(*start, in) + distance(in, *end);
			if(temp < min){
				min = temp;
				placer = end;
			}
			++start;
			++end;
		}
		total += min;
		path.insert(placer, in);
	}
	
}

//Next two functions are the MST modified for the OPTTSP, simpler and no printing but overall same logic
inline void finder(vector<loc> &in, loc next, int index, vector<vector<double>> &adj){
	for(unsigned int i = index; i < in.size(); ++i){
		
		double dist = adj[next.input][in[i].input];
		
		if(dist < in[i].distance){
				in[i].distance = dist;
				in[i].parent = next.input;
		}
	}
}

inline double estimate(deque<loc> &lin, vector<vector<double>> &adj){
		double totalWeight = 0;
		vector<loc> in(lin.begin(), lin.end());	
		for(unsigned int a = 0; a < in.size(); a++){
			loc next = least_finder(in, a, totalWeight);
			finder(in, next, a, adj);
		}
		return totalWeight;
}




//This function connects the calculated MST to the already calculated path (makes for tighter lower bound)
double pruner(deque<loc> &current, deque<loc> &rest, vector<vector<double>>& adj_matrix){
	double s_close = numeric_limits<double>::infinity();
	int s_index = -1;
	for(unsigned int i = 0; i < rest.size(); i++){
		if(adj_matrix[current.front().input][rest[i].input] < s_close){
			s_close = adj_matrix[current.front().input][rest[i].input];
			s_index = i;
		}
	}
	double b_close = numeric_limits<double>::infinity();
	int b_index = -1;
	for(unsigned int i = 0; i < rest.size(); i++){
		if(adj_matrix[current.back().input][rest[i].input] < b_close){
			b_close = adj_matrix[current.back().input][rest[i].input];
			b_index = rest[i].input;
		}
	}
	double out = s_close + b_close;
	
	//Makes sure that the connections dont connect to the same path, finds next nearest if it does
	if(s_index == b_index){
		double bn_close = numeric_limits<double>::infinity();
		double sn_close = numeric_limits<double>::infinity();
		int bn_index = -1;
		int sn_index = -1;
		for(unsigned int i = 0; i < rest.size(); i++){
			if(adj_matrix[current.back().input][rest[i].input] < bn_close && (bn_index != b_index)){
				bn_close = adj_matrix[current.back().input][rest[i].input];
				bn_index = rest[i].input;
			}
		}
		for(unsigned int i = 0; i < rest.size(); i++){
			if(adj_matrix[current.front().input][rest[i].input] < bn_close && (bn_index != b_index)){
				sn_close = adj_matrix[current.front().input][rest[i].input];
				sn_index = rest[i].input;
			}
		}
		
		if(sn_close < bn_close){
			s_index = sn_index;
			out = b_close + sn_close;
		}else{
			b_index = bn_index;
			out = s_close + bn_close; 
		}
	}
	
	//a single 2-opt for the connections
	double n_dist = adj_matrix[current.front().input][b_index] + adj_matrix[current.back().input][s_index];
	if(n_dist < out) out = n_dist;
	return out;
}

//The main genPerms function that generates all permutations of a path
void genPerms(deque<loc> &in, deque<loc> &perm, deque<loc> &best, double &u_bound, 
			double &current, double &currentBest, vector<vector<double>>& adj_matrix){
		
	unsigned int k, size = in.size();
	//Logic for when a complete permutation has been calculated
	if(in.empty()){
		double last = adj_matrix[perm.front().input][perm.back().input];
		current += last;
		
		if(current < currentBest){
			currentBest = current;
			
			if(currentBest < u_bound) u_bound = currentBest;
			
			best = perm;
		}
		current -= last;
		return;
	}
	
	//Loop that greates each permutation
	for(k = 0; k < size; ++k){
		//Add nodes to permutation
		perm.push_back(in.front());
		in.pop_front();
		double add = 0;
		//Update the running total distance
		if(perm.size() >= 2){
			add = adj_matrix[perm.back().input][perm[perm.size() - 2].input];
			current += add;
		}
		//Prunes if the lower bound is greater than the upper bound and if the permutation doesnt start with a 0
		if((current + estimate(in, adj_matrix) + pruner(perm, in, adj_matrix)) <= u_bound && perm.front().input == 0){
			genPerms(in, perm, best, u_bound, current, currentBest, adj_matrix);
		}
		//Update the running total distance and add nodes back to the original deque
		current -= add;
		in.push_back(perm.back());
		perm.pop_back();
	}
}

//Main function for the optimal tsp
void Optimal(deque<loc> &in, double &u_bound){
	//Creates adjacency matrix for efficiency, intialized to infinity
	vector<vector<double>> adj_matrix(in.size(), vector<double>(in.size(), numeric_limits<double>::infinity()));
	for(unsigned int i = 0; i < in.size(); ++i){
		for(unsigned int k = 0; k < in.size(); ++k){
			if(i != k) adj_matrix[i][k] = distance(in[i], in[k]);
		}
	}
	//Initialize current best to infinity
	double currentBest = numeric_limits<double>::infinity(), current = 0;
	deque<loc> s, best;
	//Main call of genPerms
	genPerms(in, s, best, u_bound, current, currentBest, adj_matrix);
	stringstream out;
		
	//Prints the best path
	cout << currentBest << "\n";
	while(!best.empty()){
		out << best.front().input << ' ';
		best.pop_front();
	}
	cout << out.str() << endl;	

}
	

#endif
