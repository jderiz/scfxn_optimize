#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sstream>

using namespace std;

const bool VERBOSE=false;

int
aa2idx (char AA) {
	static std::map<char,int> aamap;
	if (aamap.empty()) {
		aamap['A']=0;  aamap['C']=1;  aamap['D']=2;  aamap['E']=3;  aamap['F']=4;
		aamap['G']=5;  aamap['H']=6;  aamap['I']=7;  aamap['K']=8;  aamap['L']=9;
		aamap['M']=10; aamap['N']=11; aamap['P']=12; aamap['Q']=13; aamap['R']=14;
		aamap['S']=15; aamap['T']=16; aamap['V']=17; aamap['W']=18; aamap['Y']=19;
	}
	return aamap[AA];
}

char
idx2aa (int ii) {
	static std::map<char,int> aamap;
	if (aamap.empty()) {
		aamap[0]='A';  aamap[1]='C';  aamap[2]='D';  aamap[3]='E';  aamap[4]='F';
		aamap[5]='G';  aamap[6]='H';  aamap[7]='I';  aamap[8]='K';  aamap[9]='L';
		aamap[10]='M'; aamap[11]='N'; aamap[12]='P'; aamap[13]='Q'; aamap[14]='R';
		aamap[15]='S'; aamap[16]='T'; aamap[17]='V'; aamap[18]='W'; aamap[19]='Y';
	}
	return aamap[ii];
}

class NelderMead {
private:
	double ALPHA, BETA, GAMMA;

	std::vector< double > wts_best;
	double score_best;

public:
	NelderMead() {
		ALPHA=1;
		BETA=0.5;
		GAMMA=2.0;

		score_best = 1e6;
		wts_best.resize(20,0.0);
	}

	void
	checkpoint(
		std::vector< double > const &wts_i,
		double score_i
	) {
		if (score_i<score_best) {
			score_best = score_i;
			wts_best = wts_i;
		}
	}

	void
	recover_best(
		std::vector< double > &wts_i,
		double &score_i
	) {
		score_i = score_best;
		wts_i = wts_best;
	}

	// Helper function - find the lowest, 2nd highest and highest position
	void
	FindMarkers (
		std::vector<double> const &Y,
		int &idx_lo,
		int &idx_nhi,
		int &idx_hi
	) {
		int N = Y.size()-1;

		idx_lo=0;
		if (Y[0]>Y[1]) {
			idx_hi=0; idx_nhi=1;
		} else {
			idx_hi=1; idx_nhi=0;
		}

		for(int i=0; i<=N; ++i) {
			if (Y[i] < Y[idx_lo]) idx_lo= i;
			if (Y[i] > Y[idx_hi]) {
				idx_nhi = idx_hi; idx_hi = i;
			} else if (Y[i] > Y[idx_nhi] && idx_hi != i) {
				idx_nhi = i;
			}
		}
	}

	void CalcCentroid (
		std::vector< std::vector<double> > const&P,
		int idx_hi,
		std::vector<double> &C
	) {
		int N = P.size()-1;

		C.resize(N);
		for (int j=0; j<N; ++j) {
			C[j]=0.0;
			for(int i=0; i<=N; ++i) {
				if (i!=idx_hi) C[j] += P[i][j];
			}
			C[j] /= N;
		}
	}

	void
	AdjustCentroid (
		std::vector<double> &C,
		std::vector< std::vector<double> > const&P,
		int idx_hi_o,
		int idx_hi
	) {
		int N = P.size()-1;

		if (idx_hi_o != idx_hi) {
			for (int j=0; j<N; ++j) {
				C[j] += (P[idx_hi_o][j] - P[idx_hi][j]) / N;
			}
		}
	}

	void
	CalcReflection(
		std::vector<double> const &p1,
		std::vector<double> const &p2,
		double scale,
		std::vector<double> &p_refl
	) {
		int N = p1.size();

		p_refl.resize(N);
		for (int j=0; j<N; ++j) {
			p_refl[j] = p1[j] + scale*(p1[j]-p2[j]);
		}
	}

	void
	run(
		std::vector<double> const& guesses,
		std::vector<double> const& scales,
		double (*evaluatorFn)(std::vector<double> const&),
		double tol = 5e-4,
		int itmax = 2000,
		bool verbose = false
	) {
		int N = guesses.size();
		std::vector< std::vector<double> > vertices(N+1);
		std::vector<double> scores(N+1,0);

		// construct vertices
		vertices[0] = guesses;
		for (int i=0; i<N; ++i) {
			vertices[i+1] = guesses;
			vertices[i+1][i] += scales[i];
		}

		// evaluate vertices
		for (int i=0; i<=N; ++i) {
			scores[i] = (*evaluatorFn)(vertices[i]);
			checkpoint( vertices[i], scores[i] );
		}

		// run
		bool recalc=true;
		std::vector<double> C(N,0); // centroid
		std::vector<double> PR(N,0); // reflection point
		double YR = 0; // reflection point value
		int ihi_o = 0;

		int iter = 1;

		while (iter++<itmax) {
			int ilo, inhi, ihi;
			FindMarkers(scores, ilo, inhi, ihi);

			// Stopping conditions
			double rtol = 2*abs(scores[ihi]-scores[ilo])/(abs(scores[ihi])+abs(scores[ilo])+1e-12);
			if (rtol<tol) break;

			if (recalc) {
				CalcCentroid(vertices, ihi, C);
			} else {
				AdjustCentroid(C, vertices, ihi_o, ihi);
			}

			recalc = false;

			// Determine the reflection point, evaluate its value
			CalcReflection(C, vertices[ihi], ALPHA, PR);
			YR = (*evaluatorFn)(PR);
			checkpoint( PR, YR );

			if (YR < scores[ilo]) {
				std::vector<double> PE(N,0);
				CalcReflection(C, PR, -GAMMA, PE);
				double YE = (*evaluatorFn)(PE);
				checkpoint( PE, YE );
				if (YE < scores[ilo]) {
					vertices[ihi] = PE; scores[ihi] = YE;
				} else {
					vertices[ihi] = PR; scores[ihi] = YR;
				}
			} else if (YR >= scores[inhi]) {
				if (YR < scores[ihi] ) {
					vertices[ihi] = PR; scores[ihi] = YR;
				}

				std::vector<double> PC(N,0);
				CalcReflection(C, vertices[ihi], -BETA, PC);
				double YC=(*evaluatorFn)(PC);
				checkpoint( PC, YC );

				if (YC < scores[ihi]) {
					vertices[ihi] = PC; scores[ihi] = YC;
				} else {
					for(int i=0; i<=N; ++i) {
						if (i!=ilo) {
							std::vector<double> vi_old = vertices[i];
							CalcReflection(vertices[ilo], vi_old, -BETA, vertices[i]);
							scores[i] = (*evaluatorFn)(vertices[i]);
							checkpoint( vertices[i], scores[i] );
						}
					}
					recalc = true;
				}
			} else {
				vertices[ihi] = PR; scores[ihi] = YR;
			}

			ihi_o = ihi;
		}
	}
};

//
class PointSSMData {
public:
	PointSSMData() {
		native_ = 0;
		data_.resize(20,0.0);
	}

	double &
	operator[](int i) {
		return data_[i];
	}

	double const &
	operator[](int i) const {
		return data_[i];
	}

	void native( int i) { native_=i; }
	int native( ) { return native_; }

	void normalize() {
		double sum=0.0;
		for (int i=0; i<20; ++i) sum += data_[i];
		for (int i=0; i<20; ++i) data_[i] /= sum;
	}

private:
	int native_;
	std::vector<double> data_;
};

// the global SSM data
std::vector <PointSSMData> ENERGIES, COUNTS;

// point evaluator
double
evaluator( std::vector<double> const &values ) {
	static int ITER=1;
	double kT = 0.5, kldiv_wt = 0.5; // params

	double profile_recov=0.0, recov=0.0, score_recov=0.0, score_diverge=0.0, score=0.0;

	std::vector<double> native_distr(20,0.0);
	std::vector<double> decoy_distr(20,0.0);

	int N = ENERGIES.size();

	for (int i=0; i<N; ++i) {
		double e_min=1e6, prof_wt=0.0, probsum=0.0;
		int aa_min = -1;

		for (int j=0; j<20; ++j) {
			double e_j = ENERGIES[i][j]	+ values[j];
			if (e_min > e_j) {
				e_min = e_j;
				aa_min = j;
			}
		}

		for (int j=0; j<20; ++j) {
			double e_j = ENERGIES[i][j]	+ values[j];
			double p_j = std::exp( - (e_j - e_min)/kT );
			probsum += p_j;
			prof_wt += COUNTS[i][j] * p_j;
		}

		score_recov -= prof_wt/probsum;
		if (aa_min == COUNTS[i].native()) {
			recov+=1.0;
		}

		native_distr[ COUNTS[i].native() ]++;
		decoy_distr[ aa_min ]++;
		if (COUNTS[i][aa_min] > 0) {
			profile_recov+=1.0;
		}
	}

	for (int j=0; j<20; ++j) {
		double p_i = native_distr[j] / N;
		double q_i = decoy_distr[j] / N;
		if (p_i == 0) { p_i = 0.01/N; }
		if (q_i == 0) { q_i = 0.01/N; }
		score_diverge += p_i * log( p_i/q_i );
	}
	recov /= N;
	profile_recov /= N;
	score_recov /= N;
	score = score_recov + kldiv_wt*score_diverge;

	if (VERBOSE) {
		std::cout << "[" << ITER++ << "] recov/profile/score/kl = " <<
			recov << " / " << profile_recov << " / " << score_recov << " / " << score_diverge << " / " << score << "\n";
	}

	return (score);
}

//
int main(int ARGV, char **ARGC) {
	if (ARGV < 2) {
		std::cerr << "usage: " << ARGC[0] << " <file>\n";
		return 1;
	}

	double start[] = {2.3, 3.4, -1.8, -2.1, 1.5, 1.3, 0.3, 1.8, -1.1, 1.0, 1.4, -1.0, -2.9, -1.1, -1.3, -0.1, 1.2, 2.6, 2.9, 0.9 };
	double scale[] = { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };

	std::vector<double> startV(start, start + sizeof(start) / sizeof(start[0]) );
	std::vector<double> scaleV(scale, scale + sizeof(scale) / sizeof(scale[0]) );

	ifstream ifs (ARGC[1], ifstream::in);

	for (int arg=1;arg<ARGV; ++arg) {
		std::string tag(ARGC[arg]);

		// temp storage for 1 struct worth of data
		std::vector <PointSSMData> energies_i, counts_i;
		std::vector <int> resids;
		std::string line;

		std::string scorefile = tag+".ENERGIES";
		ifstream ifener (scorefile.c_str(), ifstream::in);
		std::string buf;

		while ( getline (ifener,line) ) {
			std::stringstream ss(line);
			std::vector<string> tokens; // Create vector to hold our words
			while (ss >> buf)
				tokens.push_back(buf);

			if (tokens.size() < 22) continue;

			PointSSMData energies_ij;
			energies_ij.native( atoi(tokens[1].c_str()) );
			for (int i=0; i<20; ++i) {
				energies_ij[i] = atof(tokens[i+2].c_str());
			}
			energies_i.push_back( energies_ij );
		}

		std::string countfile = tag+".COUNTS";
		ifstream ifcount (countfile.c_str(), ifstream::in);

		while ( getline (ifcount,line) ) {
			stringstream ss(line);
			std::vector<string> tokens; // Create vector to hold our words
			while (ss >> buf)
			tokens.push_back(buf);

			if (tokens.size() < 22) continue;

			PointSSMData counts_ij;
			counts_ij.native( aa2idx(tokens[1][0]) );
			for (int i=0; i<20; ++i) {
				counts_ij[i] = atof(tokens[i+2].c_str());
			}
			counts_ij.normalize();
			counts_i.push_back( counts_ij );
		}

		for (int i=0; i<energies_i.size(); ++i) {
			if (counts_i.size() > energies_i[i].native()-1) {
				ENERGIES.push_back( energies_i[i] );
				COUNTS.push_back( counts_i[energies_i[i].native()-1] );
			}
		}
	}


	int N = ENERGIES.size();
	if (VERBOSE) std::cout << "Read " << N << " residues.\n";

	NelderMead optimizer;
	optimizer.run( startV, scaleV, &evaluator );

	std::vector< double > final_ref;
	double final_score;
	optimizer.recover_best(final_ref,final_score);

	std::cout << final_score << " " << N << std::endl << std::endl;

	std::cout <<  "METHOD_WEIGHTS ref " <<std::setprecision(5);
	for (int i=0; i<20; ++i) {
		std::cout << final_ref[i] << " ";
	}
	std::cout << "\n";

}

